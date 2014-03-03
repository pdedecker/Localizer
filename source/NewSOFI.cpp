/*
 Copyright 2008-2014 Peter Dedecker.
 
 This file is part of Localizer.
 
 Localizer is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 Localizer is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with Localizer.  If not, see <http://www.gnu.org/licenses/>.
 
 
 Additional permission under GNU GPL version 3 section 7
 
 If you modify this Program, or any covered work, by
 linking or combining it with libraries required for interaction
 with analysis programs such as Igor Pro or Matlab, or to acquire
 data from or control hardware related to an experimental measurement,
 the licensors of this Program grant you additional permission
 to convey the resulting work.
 */

#include "NewSOFI.h"

#include <map>
#include <sstream>

#include <eigen3/Eigen/Eigen>
#include "tbb/tbb.h"
#include "tbb/spin_mutex.h"
#include "XOPStandardHeaders.h"

#include "NewSOFIKernels.h"

int RawSOFIWorker(std::shared_ptr<ImageLoader> imageLoader, const std::vector<std::shared_ptr<SOFIFrameVerifier> >& frameVerifiers, const int firstImageToProcess, const int lastImageToProcess, int &imagesProcessedSoFar, const int& totalNumberOfImagesToProcess, std::shared_ptr<ProgressReporter> progressReporter, const std::vector<std::pair<int, std::vector<SOFIKernel> > >& orders, std::map<PixelCombination,ImagePtr,ComparePixelCombinations>& pixelMap, std::vector<ImagePtr>& sofiImages, bool wantAverageImage, ImagePtr& averageImage);

double Prefactor(int nPartitions);
Eigen::MatrixXd EvaluatePartition(const Partition& partition, const std::map<PixelCombination,ImagePtr,ComparePixelCombinations>& pixelMap, const int nOutputRows, const int nOutputCols);
Eigen::MatrixXd EvaluatePartitionsSet(const GroupOfPartitions& groupOfPartitions, const std::map<PixelCombination,ImagePtr,ComparePixelCombinations>& pixelMap, const int nOutputRows, const int nOutputCols);

void PerformPixelationCorrection(ImagePtr imageToCorrect, const int order);

void DoNewSOFI(std::shared_ptr<ImageLoader> imageLoader, SOFIOptions& options, std::shared_ptr<ProgressReporter> progressReporter, std::vector<ImagePtr>& sofiOutputImages) {
    int nImages = imageLoader->getNImages();
    if (nImages == 0)
        throw std::runtime_error("SOFI without input images");
    
    int order = options.order;
    const std::vector<std::shared_ptr<SOFIFrameVerifier> >& frameVerifiers = options.frameVerifiers;
    bool wantAverageImage = options.wantAverageImage;
    ImagePtr& averageImage = options.averageImage;
    
    std::vector<std::pair<int, std::vector<SOFIKernel> > > orders;
    orders.push_back(std::pair<int, std::vector<SOFIKernel> >(order, KernelsForOrder(order)));
    size_t nOrders = orders.size();
    int firstImageToProcess = 0;
    int lastImageToProcess = 0;
    int imagesProcessedSoFar = 0;
    int totalNumberOfImagesToProcess = nImages;
    std::map<PixelCombination,ImagePtr,ComparePixelCombinations> pixelMap;
    std::vector<std::vector<ImagePtr> > allSubImages(nOrders);
    std::vector<int> nImagesInSubCalculation;
    int totalNumberOfImagesIncluded = 0;
    if (wantAverageImage) {
        averageImage = ImagePtr(new Image(imageLoader->getXSize(), imageLoader->getYSize()));
        averageImage->setConstant(0.0);
    }
    
    progressReporter->CalculationStarted();
    
    // calculate a SOFI image for every batch of images
    int batchSize = 100;
    int nBatches;
    if (nImages % batchSize == 0) {
        nBatches = nImages / batchSize;
    } else {
        nBatches = nImages / batchSize + 1;
    }
    for (int n = 0; n < nBatches; ++n) {
        firstImageToProcess = lastImageToProcess;
        lastImageToProcess = std::min(firstImageToProcess + batchSize - 1, nImages - 1);
        std::vector<ImagePtr> subImages;
        int nImagesIncluded = RawSOFIWorker(imageLoader, frameVerifiers, firstImageToProcess, lastImageToProcess, imagesProcessedSoFar, totalNumberOfImagesToProcess, progressReporter, orders, pixelMap, subImages, wantAverageImage, averageImage);
        nImagesInSubCalculation.push_back(nImagesIncluded);
        totalNumberOfImagesIncluded += nImagesIncluded;
        for (size_t j = 0; j < nOrders; ++j) {
            allSubImages[j].push_back(subImages[j]);
        }
    }
    
    // combine the batches into a single average image per order and perform pixelation correction
    std::vector<ImagePtr> mergedImages(nOrders);
    tbb::parallel_for<size_t>(0, nOrders, [&](size_t j) {
        int order = orders.at(j).first;
        const std::vector<ImagePtr>& theseSubImages = allSubImages.at(j);
        ImagePtr mergedImage(new Image(theseSubImages.at(0)->rows(), theseSubImages.at(0)->cols()));
        mergedImage->setConstant(0.0);
        double summedContribution = 0.0;
        for (int batch = 0; batch < nBatches; ++batch) {
            double contribution = static_cast<double>(nImagesInSubCalculation.at(batch)) / static_cast<double>(batchSize);
            *mergedImage += contribution * *(theseSubImages.at(batch));
            summedContribution += contribution;
        }
        *mergedImage /= summedContribution;
        
        if (options.doPixelationCorrection)
            PerformPixelationCorrection(mergedImage, order);
        
        mergedImages.at(j) = mergedImage;
    });
    
    sofiOutputImages = mergedImages;
    
    if (wantAverageImage)
        *averageImage /= static_cast<double>(totalNumberOfImagesIncluded);
    
    progressReporter->CalculationDone();
}

void MinMaxDeltas(const PixelCombination& pixelCombination, int& minRowDelta, int& maxRowDelta, int& minColDelta, int& maxColDelta) {
    minRowDelta = 1000;
    minColDelta = 1000;
    maxRowDelta = -1000;
    maxColDelta = -1000;
    
    for (size_t i = 0; i < pixelCombination.size(); ++i) {
        const std::pair<int,int>& pixel = pixelCombination[i];
        if (pixel.first < minRowDelta)
            minRowDelta = pixel.first;
        if (pixel.first > maxRowDelta)
            maxRowDelta = pixel.first;
        if (pixel.second < minColDelta)
            minColDelta = pixel.second;
        if (pixel.second > maxColDelta)
            maxColDelta = pixel.second;
    }
}

PixelCombination OffsetPixelCombination(const PixelCombination& pixelCombination, const int rowDelta, const int colDelta) {
    PixelCombination offsetCombination(pixelCombination);
    for (size_t i = 0; i < offsetCombination.size(); ++i) {
        offsetCombination[i].first += rowDelta;
        offsetCombination[i].second += colDelta;
    }
    return offsetCombination;
}

int RawSOFIWorker(std::shared_ptr<ImageLoader> imageLoader, const std::vector<std::shared_ptr<SOFIFrameVerifier> >& frameVerifiers, const int firstImageToProcess, const int lastImageToProcess, int &imagesProcessedSoFar, const int& totalNumberOfImagesToProcess, std::shared_ptr<ProgressReporter> progressReporter, const std::vector<std::pair<int, std::vector<SOFIKernel> > >& orders, std::map<PixelCombination,ImagePtr,ComparePixelCombinations>& pixelMap, std::vector<ImagePtr>& sofiImages, bool wantAverageImage, ImagePtr& averageImage) {
    int nRows = imageLoader->getXSize();
    int nCols = imageLoader->getYSize();
    int nImagesToInclude = lastImageToProcess - firstImageToProcess + 1;
    int nImagesIncluded = 0;
    
    // store all needed pixel combinations in the map
    // if the map is non-empty, we assume that it was already setup by a previous call to RawSOFIWorker, so we can prevent unnecessary work.
    if (pixelMap.empty()) {
        int allCombinations = 0;
        for (size_t i = 0; i < orders.size(); ++i) {
            const std::vector<SOFIKernel>& kernels = orders[i].second;
            for (auto kernelIt = kernels.cbegin(); kernelIt != kernels.cend(); ++kernelIt) {
                const std::vector<GroupOfPartitions>& GroupOfPartitions = kernelIt->combinations;
                for (auto partitionsSetIt = GroupOfPartitions.cbegin(); partitionsSetIt != GroupOfPartitions.cend(); ++partitionsSetIt) {
                    for (auto partitionsIt = partitionsSetIt->cbegin(); partitionsIt != partitionsSetIt->cend(); ++partitionsIt) {
                        for (auto pixelCombinationIt = partitionsIt->cbegin(); pixelCombinationIt != partitionsIt->cend(); ++pixelCombinationIt) {
                            allCombinations += 1;
                            int minRowDelta, maxRowDelta, minColDelta, maxColDelta;
                            MinMaxDeltas(*pixelCombinationIt, minRowDelta, maxRowDelta, minColDelta, maxColDelta);
                            PixelCombination offsetCombination = OffsetPixelCombination(*pixelCombinationIt, -1 * (2 + minRowDelta), -1 * (2 + minColDelta));
                            if (!pixelMap.count(offsetCombination)) {
                                maxRowDelta -= 2 + minRowDelta;
                                maxColDelta -= 2 + minColDelta;
                                int nRowsToCalculate = nRows - maxRowDelta;
                                int nColsToCalculate = nCols - maxColDelta;
                                ImagePtr matrix(new Image(nRowsToCalculate, nColsToCalculate));
                                pixelMap.insert(std::pair<PixelCombination, ImagePtr>(offsetCombination, matrix));
                            }
                        }
                    }
                }
            }
        }
        //std::ostringstream ostream;
        //ostream << "Map has total of " << pixelMap.size() << " entries (" << static_cast<int>(static_cast<double>(pixelMap.size()) / static_cast<double>(allCombinations) * 100.0) << "% retained)\r";
        //XOPNotice(ostream.str().c_str());
    }
    
    // clear all accumulated pixel combinations
    tbb::parallel_do(pixelMap.begin(), pixelMap.end(), [](std::pair<PixelCombination, ImagePtr> item) {
        ImagePtr matrix = item.second;
        matrix->setConstant(0.0);
    });
    
    // calculate all products over the images
    for (int n = firstImageToProcess; n <= lastImageToProcess; ++n) {
        // update progress and check for abort
        int abortStatus = progressReporter->UpdateCalculationProgress(imagesProcessedSoFar, totalNumberOfImagesToProcess);
        if (abortStatus)
            throw USER_ABORTED("user abort");
        
        imagesProcessedSoFar += 1;
        ImagePtr currentImage = imageLoader->readImage(n);
        
        // check that the frame is valid
        bool isValidFrame = true;
        for (auto verifierIt = frameVerifiers.cbegin(); verifierIt != frameVerifiers.cend(); ++verifierIt) {
            if (!(*verifierIt)->isValidFrame(currentImage)) {
                isValidFrame = false;
                break;
            }
        }
        if (!isValidFrame)
            continue;
        nImagesIncluded += 1;
        if (wantAverageImage)
            *averageImage += *currentImage;
        
        // do the actual calculation
        tbb::parallel_do(pixelMap.begin(), pixelMap.end(), [=](std::pair<PixelCombination,ImagePtr> item) {
            const PixelCombination& currentCombination = item.first;
            ImagePtr matrix = item.second;
            int nRowsToCalculate = matrix->rows();
            int nColsToCalculate = matrix->cols();
            for (int col = 2; col < nColsToCalculate; ++col) {
                for (int row = 2; row < nRowsToCalculate; ++row) {
                    double product = 1.0;
                    for (size_t i = 0; i < currentCombination.size(); ++i) {
                        product *= (*currentImage)(row + currentCombination[i].first, col + currentCombination[i].second);
                    }
                    (*matrix)(row - 2, col - 2) += product;
                }
            }
        });
    }
    
    // normalize by number of images
    tbb::parallel_do(pixelMap.begin(), pixelMap.end(), [=](std::pair<PixelCombination, ImagePtr> item) {
        ImagePtr matrix = item.second;
        (*matrix) /= static_cast<double>(nImagesToInclude);
    });
    
    // and make the SOFI images
    sofiImages.clear();
    for (size_t i = 0; i < orders.size(); ++i) {
        const std::pair<int, std::vector<SOFIKernel> >& calculation = orders[i];
        int order = calculation.first;
        const std::vector<SOFIKernel>& kernels = calculation.second;
        ImagePtr sofiImage(new Image(order * (nRows - 4), order * (nCols - 4)));
        tbb::parallel_do(kernels.cbegin(), kernels.cend(), [=,&pixelMap,&sofiImage](const SOFIKernel& kernel) {
            Eigen::MatrixXd evaluated(nRows - 4, nCols - 4);
            evaluated.setConstant(0.0);
            for (auto partitionsSetIt = kernel.combinations.cbegin(); partitionsSetIt != kernel.combinations.cend(); ++partitionsSetIt) {
                evaluated += EvaluatePartitionsSet(*partitionsSetIt, pixelMap, nRows - 4, nCols - 4);
            }
            if (kernel.combinations.size() > 1)
                evaluated /= static_cast<double>(kernel.combinations.size());
            
            for (int col = 0; col < nCols - 4; ++col) {
                for (int row = 0; row < nRows - 4; ++row) {
                    int baseOutputRow = row * order;
                    int baseOutputCol = col * order;
                    (*sofiImage)(baseOutputRow + kernel.outputDeltaX, baseOutputCol + kernel.outputDeltaY) = evaluated(row, col);
                }
            }
        });
        sofiImages.push_back(sofiImage);
    }
    
    return nImagesIncluded;
}

Eigen::MatrixXd EvaluatePartitionsSet(const GroupOfPartitions& groupOfPartitions, const std::map<PixelCombination,ImagePtr,ComparePixelCombinations>& pixelMap, const int nOutputRows, const int nOutputCols) {
    Eigen::MatrixXd accumulated(nOutputRows, nOutputCols);
    accumulated.setConstant(0.0);
    for (auto it = groupOfPartitions.cbegin(); it != groupOfPartitions.cend(); ++it) {
        int nSetsInPartition = it->size();
        accumulated += Prefactor(nSetsInPartition) * EvaluatePartition(*it, pixelMap, nOutputRows, nOutputCols);
    }
    
    return accumulated;
}

Eigen::MatrixXd EvaluatePartition(const Partition& partition, const std::map<PixelCombination,ImagePtr,ComparePixelCombinations>& pixelMap, const int nOutputRows, const int nOutputCols) {
    Eigen::MatrixXd result(nOutputRows, nOutputCols);
    result.setConstant(1.0);
    int minRowDelta, maxRowDelta, minColDelta, maxColDelta;
    for (size_t i = 0; i < partition.size(); ++i) {
        const PixelCombination& subset = partition[i];
        MinMaxDeltas(subset, minRowDelta, maxRowDelta, minColDelta, maxColDelta);
        PixelCombination offsetCombination = OffsetPixelCombination(subset, -1 * (2 + minRowDelta), -1 * (2 + minColDelta));
        const Image& productMatrix = *(pixelMap.at(offsetCombination));
        for (int col = 0; col < nOutputCols; ++col) {
            for (int row = 0; row < nOutputRows; ++row) {
                result(row, col) *= productMatrix(row + minRowDelta + 2, col + minColDelta + 2);
            }
        }
    }
    
    return result;
}

double Prefactor(int nPartitions) {
    switch (nPartitions) {
        case 1:
            return 1.0;
        case 2:
            return -1.0;
        case 3:
            return 2.0;
        case 4:
            return -6.0;
        case 5:
            return 24.0;
        case 6:
            return -120.0;
        case 7:
            return 720.0;
        case 8:
            return -5040.0;
        case 9:
            return 40320.0;
        default:
            throw std::runtime_error("no prefactor");
            return 1.0;
    }
}

void PerformPixelationCorrection(ImagePtr imageToCorrect, const int order) {
    int nRows = imageToCorrect->rows();
    int nCols = imageToCorrect->cols();
    
    Eigen::ArrayXXd aFactor(order, order), bTerm(order, order);
    int nKernelRows = order, nKernelCols = order;
    int nPixelsOfEachKind = nRows * nCols / (order * order);
    int kindOfRow, kindOfCol;
    
    // calculate averages
    Eigen::ArrayXXd pixelAverages(order, order);
    pixelAverages.setConstant(0.0);
    for (int col = 0; col < nCols; ++col) {
        for (int row = 0; row < nRows; ++row) {
            kindOfRow = row % nKernelRows;
            kindOfCol = col % nKernelCols;
            pixelAverages(kindOfRow, kindOfCol) += (*imageToCorrect)(row, col);
        }
    }
    pixelAverages /= static_cast<double>(nPixelsOfEachKind);
    
    // calculate variances
    Eigen::ArrayXXd pixelVariances(order, order);
    pixelVariances.setConstant(0.0);
    for (int col = 0; col < nCols; ++col) {
        for (int row = 0; row < nRows; ++row) {
            kindOfRow = row % nKernelRows;
            kindOfCol = col % nKernelCols;
            pixelVariances(kindOfRow, kindOfCol) += square<double>((*imageToCorrect)(row, col) - pixelAverages(kindOfRow, kindOfCol));
        }
    }
    pixelVariances /= static_cast<double>(nPixelsOfEachKind - 1);
    
    // now determine correction factors, arbitrary to the first element
    for (int col = 0; col < nKernelCols; ++col) {
        for (int row = 0; row < nKernelRows; ++row) {
            if (row == 0 && col == 0) {
                aFactor(row, col) = 1.0;
                bTerm(row, col) = 0.0;
                continue;
            }
            aFactor(row, col) = std::sqrt(pixelVariances(0, 0) / pixelVariances(row, col));
            bTerm(row, col) = pixelAverages(0, 0) - aFactor(row, col) * pixelAverages(row, col);
        }
    }
    
    // perform the correction
    for (int col = 0; col < nCols; ++col) {
        for (int row = 0; row < nRows; ++row) {
            kindOfRow = row % nKernelRows;
            kindOfCol = col % nKernelCols;
            (*imageToCorrect)(row, col) = (*imageToCorrect)(row, col) * aFactor(kindOfRow, kindOfCol) + bTerm(kindOfRow, kindOfCol);
        }
    }
}
