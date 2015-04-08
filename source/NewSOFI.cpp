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
#include <boost/circular_buffer.hpp>
#include "tbb/tbb.h"
#include "tbb/spin_mutex.h"
#ifdef WITH_IGOR
#include "XOPStandardHeaders.h"
#include "IgorUtilities.h"
#endif

#include "NewSOFIKernels.h"

int RawSOFIWorker(std::shared_ptr<ImageLoader> imageLoader, const std::vector<std::shared_ptr<SOFIFrameVerifier> >& frameVerifiers, const int firstImageToProcess, const int lastImageToProcess, int &imagesProcessedSoFar, const int& totalNumberOfImagesToProcess, std::shared_ptr<ProgressReporter> progressReporter, const std::vector<std::pair<int, std::vector<SOFIKernel> > >& orders, const std::vector<int>& timeLags, bool isAuto, std::map<PixelCombination,ImagePtr,ComparePixelCombinations>& pixelMap, std::vector<ImagePtr>& sofiImages, bool wantAverageImage, ImagePtr& averageImage, bool wantJackKnife, std::vector<std::vector<ImagePtr> >& jackKnifeImages, const std::vector<double>& pixelCombinationWeights, bool wantDebugMessages);

ImagePtr AssembleSOFIImage(const int nInputRows, const int nInputCols, const int order, const bool isAuto, const std::vector<SOFIKernel>& kernels, const std::map<PixelCombination,ImagePtr,ComparePixelCombinations>& pixelMap, const std::vector<double>& pixelCombinationWeights, std::vector<std::vector<double>>& usedCombinationWeights);
double Prefactor(int nPartitions);
Eigen::MatrixXd EvaluatePartition(const Partition& partition, const std::map<PixelCombination,ImagePtr,ComparePixelCombinations>& pixelMap, const int nOutputRows, const int nOutputCols);
Eigen::MatrixXd EvaluatePartitionsSet(const GroupOfPartitions& groupOfPartitions, const std::map<PixelCombination,ImagePtr,ComparePixelCombinations>& pixelMap, const int nOutputRows, const int nOutputCols);

void PerformPixelationCorrection(ImagePtr imageToCorrect, bool alsoCorrectVariance, const int order);
void AssembleJackKnifeImages(std::vector<std::vector<std::vector<ImagePtr> > >& jackKnifeBatchSOFIImages, const std::vector<int>& imagesIncludedInBatch, const int batchSize, const std::vector<std::vector<ImagePtr> >& batchSOFIImages, std::vector<std::vector<ImagePtr> >& jackKnifeOutputImages);

void JackKnife(std::shared_ptr<ImageLoader> imageLoader, const int firstImageToInclude, const int lastImageToInclude, const int order, const std::vector<int>& timeLags, const bool isAuto, const std::vector<SOFIKernel>& kernels, std::map<PixelCombination,ImagePtr,ComparePixelCombinations>& pixelMap, const std::vector<double>& pixelCombinationWeights, std::vector<ImagePtr>& jackKnifeImages);

int NumberOfPixelCombinationsInKernels(const std::vector<SOFIKernel>& kernels);
void PrintVirtualPixelInfo(const std::vector<SOFIKernel>& kernels, const std::vector<std::vector<double>>& combinationWeights);

void DoNewSOFI(std::shared_ptr<ImageLoader> imageLoader, SOFIOptions& options, std::shared_ptr<ProgressReporter> progressReporter, std::vector<ImagePtr>& sofiOutputImages) {
    int nImages = imageLoader->getNImages();
    if (nImages <= 2)
        throw std::runtime_error("too few input images");
    if ((imageLoader->getXSize() < 5) || (imageLoader->getYSize() < 5))
        throw std::runtime_error("images too small");
    
    const std::vector<int> orders = options.orders;
    const bool isAuto = !options.wantCrossCumulant;
    std::vector<int> lagTimes = options.lagTimes;
    int nOrders = orders.size();
    int requestedBatchSize = options.batchSize;
    const std::vector<std::shared_ptr<SOFIFrameVerifier> >& frameVerifiers = options.frameVerifiers;
    bool wantAverageImage = options.wantAverageImage;
    ImagePtr& averageImage = options.averageImage;
    ImagePtr batchAverageImage;
    bool wantJackKnife = options.wantJackKnife;
    bool wantDebugMessages = options.wantDebugMessages;
    std::vector<std::vector<ImagePtr> >& jackKnifeOutputImages = options.jackKnifeImages;
    const std::vector<double>& pixelCombinationWeights = options.pixelCombinationWeights;
    bool haveExplicitPixelCombinationWeights = !pixelCombinationWeights.empty();
    jackKnifeOutputImages.clear();
    std::vector<std::vector<ImagePtr> > batchSOFIImages;    // SOFI image of each batch, 1 vector per order. Only used for jackknife.
    batchSOFIImages.resize(nOrders);
    std::vector<std::vector<std::vector<ImagePtr> > > jackKnifeBatchSOFIImages;  // jackknife SOFI images. batch / order / images
    std::vector<int> imagesIncludedInBatch;
    
    if (!lagTimes.empty()) {
        // have explicit lag times
        if (nOrders > 1)
            throw std::runtime_error("only a single order can be calculated if time lags are specified");
        if (orders.at(0) == 1)
            throw std::runtime_error("lag times do not make sense for a first order calculation");
        if (lagTimes.size() != (orders.at(0) - 1))
            throw std::runtime_error("SOFI calculation of order n requires specification of exactly (n - 1) lag times");
        // prepend zero lag time
        lagTimes.insert(lagTimes.begin(), 0);
    } else {
        // take all zero lag times
        int largestOrder = *std::max_element(orders.cbegin(), orders.cend());
        lagTimes = std::vector<int>(largestOrder, 0);
    }
    
    std::vector<std::pair<int, std::vector<SOFIKernel> > > kernelPairs;
    for (int i = 0; i < nOrders; ++i) {
        if (isAuto) {
            kernelPairs.push_back(std::pair<int, std::vector<SOFIKernel> >(orders[i], AutoKernelsForOrder(orders[i], lagTimes)));
        } else {
            kernelPairs.push_back(std::pair<int, std::vector<SOFIKernel> >(orders[i], KernelsForOrder(orders[i], lagTimes, options.pixelCombinationCutoff)));
        }
    }
    
    // if weights are specified then we can only calculate a single order
    // and the number of points needs to match.
    if (haveExplicitPixelCombinationWeights) {
        if (nOrders > 1)
            throw std::runtime_error("only a single order can be calculated if weights are specified");
        if (pixelCombinationWeights.size() != NumberOfPixelCombinationsInKernels(kernelPairs.at(0).second))
            throw std::runtime_error("number of weights does not match number of pixel combinations");
    }
    
    int firstImageToProcess = 0;
    int lastImageToProcess = 0;
    int imagesProcessedSoFar = 0;
    int totalNumberOfImagesToProcess = nImages;
    std::map<PixelCombination,ImagePtr,ComparePixelCombinations> pixelMap;
    std::vector<std::vector<ImagePtr> > allSubImages(nOrders);
    double summedContributions = 0.0;
    int totalNumberOfImagesIncluded = 0;
    if (wantAverageImage) {
        averageImage = ImagePtr(new Image(imageLoader->getXSize(), imageLoader->getYSize()));
        averageImage->setConstant(0.0);
        batchAverageImage = ImagePtr(new Image(*averageImage));
    }
    
    progressReporter->CalculationStarted();
    // calculate a SOFI image for every batch of images
    std::vector<ImagePtr> mergedImages(nOrders);
    int batchSize = (requestedBatchSize > 0) ? requestedBatchSize : 100;
    int nBatches;
    if (nImages % batchSize <= 2) {
        nBatches = nImages / batchSize;
    } else {
        nBatches = nImages / batchSize + 1;
    }
    jackKnifeBatchSOFIImages.resize(nBatches);
    for (unsigned int n = 0; n < nBatches; ++n) {
        firstImageToProcess = lastImageToProcess;
        lastImageToProcess = std::min(firstImageToProcess + batchSize - 1, nImages - 1);
        std::vector<ImagePtr> subImages;
        int nImagesIncluded = RawSOFIWorker(imageLoader, frameVerifiers, firstImageToProcess, lastImageToProcess, imagesProcessedSoFar, totalNumberOfImagesToProcess, progressReporter, kernelPairs, lagTimes, isAuto, pixelMap, subImages, wantAverageImage, batchAverageImage, wantJackKnife, jackKnifeBatchSOFIImages.at(n), pixelCombinationWeights, wantDebugMessages);
        totalNumberOfImagesIncluded += nImagesIncluded;
        double contribution = static_cast<double>(nImagesIncluded) / static_cast<double>(batchSize);
        summedContributions += contribution;
        for (size_t j = 0; j < nOrders; ++j) {
            if (mergedImages[j].get() == NULL) {
                mergedImages[j] = ImagePtr(new Image(subImages[j]->rows(), subImages[j]->cols()));
                mergedImages[j]->setConstant(0.0);
            }
            *mergedImages[j] += *subImages[j] * contribution;
        }
        if (wantAverageImage) {
            *averageImage += *batchAverageImage * contribution;
        }
        if (wantJackKnife) {
            for (int j = 0; j < nOrders; ++j) {
                batchSOFIImages.at(j).push_back(subImages[j]);
            }
            imagesIncludedInBatch.push_back(nImagesIncluded);
        }
    }
    
    // create output images and perform pixelation correction
    tbb::parallel_for<int>(0, nOrders, [&](int j) {
        *mergedImages[j] /= summedContributions;
        
        if (options.doPixelationCorrection)
            PerformPixelationCorrection(mergedImages[j], options.alsoCorrectVariance, orders[j]);
    });
    
    sofiOutputImages = mergedImages;
    
    if (wantAverageImage)
        *averageImage /= summedContributions;
    
    if (wantJackKnife) {
        AssembleJackKnifeImages(jackKnifeBatchSOFIImages, imagesIncludedInBatch, batchSize, batchSOFIImages, jackKnifeOutputImages);
        if (options.doPixelationCorrection) {
            for (int orderIndex = 0; orderIndex < nOrders; ++orderIndex) {
                std::vector<ImagePtr>& imagesForThisOrder = jackKnifeOutputImages.at(orderIndex);
                tbb::parallel_for<size_t>(0, imagesForThisOrder.size(), [&](size_t i) {
                    PerformPixelationCorrection(imagesForThisOrder[i], options.alsoCorrectVariance, orders.at(orderIndex));
                });
            }
        }
    }
    
    progressReporter->CalculationDone();
}

void MinMaxDeltas(const PixelCombination& pixelCombination, int& minRowDelta, int& maxRowDelta, int& minColDelta, int& maxColDelta) {
    minRowDelta = 1000;
    minColDelta = 1000;
    maxRowDelta = -1000;
    maxColDelta = -1000;
    
    for (size_t i = 0; i < pixelCombination.size(); ++i) {
        const std::tuple<int,int, int>& pixel = pixelCombination[i];
        if (std::get<0>(pixel) < minRowDelta)
            minRowDelta = std::get<0>(pixel);
        if (std::get<0>(pixel) > maxRowDelta)
            maxRowDelta = std::get<0>(pixel);
        if (std::get<1>(pixel) < minColDelta)
            minColDelta = std::get<1>(pixel);
        if (std::get<1>(pixel) > maxColDelta)
            maxColDelta = std::get<1>(pixel);
    }
}

PixelCombination OffsetPixelCombination(const PixelCombination& pixelCombination, const int rowDelta, const int colDelta) {
    PixelCombination offsetCombination(pixelCombination);
    for (size_t i = 0; i < offsetCombination.size(); ++i) {
        std::get<0>(offsetCombination[i]) += rowDelta;
        std::get<1>(offsetCombination[i]) += colDelta;
    }
    return offsetCombination;
}

void AssembleJackKnifeImages(std::vector<std::vector<std::vector<ImagePtr> > >& jackKnifeBatchSOFIImages, const std::vector<int>& imagesIncludedInBatch, const int batchSize, const std::vector<std::vector<ImagePtr> >& batchSOFIImages, std::vector<std::vector<ImagePtr> >& jackKnifeOutputImages) {
    int nBatches = imagesIncludedInBatch.size();
    int nOrders = batchSOFIImages.size();
    
    jackKnifeOutputImages.resize(nOrders);
    
    tbb::parallel_for<int>(0, nOrders, [&](int orderIndex) {
        for (int thisBatch = 0; thisBatch < nBatches; ++thisBatch) {
            std::vector<ImagePtr>& jackKnifeImagesInBatch = jackKnifeBatchSOFIImages.at(thisBatch).at(orderIndex);
            double contribution = static_cast<double>(imagesIncludedInBatch[thisBatch] - 1) / static_cast<double>(batchSize);
            for (int thisImage = 0; thisImage < jackKnifeImagesInBatch.size(); ++ thisImage) {
                *(jackKnifeImagesInBatch.at(thisImage)) *= contribution;
                double summedContribution = contribution;
                for (int otherBatch = 0; otherBatch < nBatches; ++otherBatch) {
                    if (otherBatch == thisBatch)
                        continue;
                    double otherContribution = static_cast<double>(imagesIncludedInBatch[otherBatch]) / static_cast<double>(batchSize);
                    *(jackKnifeImagesInBatch.at(thisImage)) += (*(batchSOFIImages.at(orderIndex).at(otherBatch))) * otherContribution;
                    summedContribution += otherContribution;
                }
                *(jackKnifeImagesInBatch.at(thisImage)) /= summedContribution;
                jackKnifeOutputImages.at(orderIndex).push_back(jackKnifeImagesInBatch.at(thisImage));
            }
        }
    });
}

int RawSOFIWorker(std::shared_ptr<ImageLoader> imageLoader, const std::vector<std::shared_ptr<SOFIFrameVerifier> >& frameVerifiers, const int firstImageToProcess, const int lastImageToProcess, int &imagesProcessedSoFar, const int& totalNumberOfImagesToProcess, std::shared_ptr<ProgressReporter> progressReporter, const std::vector<std::pair<int, std::vector<SOFIKernel> > >& orders, const std::vector<int>& timeLags, bool isAuto, std::map<PixelCombination,ImagePtr,ComparePixelCombinations>& pixelMap, std::vector<ImagePtr>& sofiImages, bool wantAverageImage, ImagePtr& averageImage, bool wantJackKnife, std::vector<std::vector<ImagePtr> >& jackKnifeImages, const std::vector<double>& pixelCombinationWeights, bool wantDebugMessages) {
    int nRows = imageLoader->getXSize();
    int nCols = imageLoader->getYSize();
    int nImagesIncluded = 0;
    sofiImages.clear();
    jackKnifeImages.clear();
    
    // only one order for explicit combination weights allowed
    if (!pixelCombinationWeights.empty() && (orders.size() != 1))
        throw std::logic_error("explicit weights but more than one order");
    
    // find min and max time lag
    int minTimeLag = 0, maxTimeLag = 0;
    for (size_t i = 0; i < timeLags.size(); i++) {
        minTimeLag = std::min(minTimeLag, timeLags[i]);
        maxTimeLag = std::max(maxTimeLag, timeLags[i]);
    }
    // required: lowest time lag is zero, all others are zero or postive
    if (minTimeLag != 0)
        throw std::runtime_error("at least one time lag must be zero, and all time lags must be zero or positive");
    int nImagesInBuffer = maxTimeLag - minTimeLag + 1;
    bool haveNonZeroTimeLag = ((minTimeLag != 0) || (maxTimeLag != 0));
    
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
    
    if (wantAverageImage) {
        if ((averageImage.get() == NULL) || (averageImage->rows() != nRows) || (averageImage->cols() != nCols)) {
            averageImage = ImagePtr(new Image(nRows, nCols));
        }
        averageImage->setConstant(0.0);
    }
    
    // initialize the buffer with the images.
    imageLoader->spoolTo(firstImageToProcess);
    boost::circular_buffer<ImagePtr> imageBuffer(nImagesInBuffer);
    for (int i = 0; (i < nImagesInBuffer - 1) && (i <= lastImageToProcess); i++) {
        imageBuffer.push_back(imageLoader->readNextImage());
    }
    
    // calculate all products over the images
    int firstImageWithOutput = firstImageToProcess;
    int lastImageWithOutput = lastImageToProcess - maxTimeLag;
    imagesProcessedSoFar += firstImageWithOutput - firstImageToProcess;
    for (int n = firstImageWithOutput; n <= lastImageWithOutput; ++n) {
        // update progress and check for abort
        int abortStatus = progressReporter->UpdateCalculationProgress(imagesProcessedSoFar, totalNumberOfImagesToProcess);
        if (abortStatus)
            throw USER_ABORTED("user abort");
        
        ImagePtr currentImage = imageLoader->readNextImage();
        imagesProcessedSoFar += 1;
        
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
        
        imageBuffer.push_back(currentImage);
        
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
                    if (haveNonZeroTimeLag) {
                        for (size_t i = 0; i < currentCombination.size(); ++i) {
                            const ImagePtr& imageAtTimeLag = imageBuffer.at(std::get<2>(currentCombination[i]));
                            product *= (*imageAtTimeLag)(row + std::get<0>(currentCombination[i]), col + std::get<1>(currentCombination[i]));
                        }
                    } else {
                        for (size_t i = 0; i < currentCombination.size(); ++i) {
                            product *= (*currentImage)(row + std::get<0>(currentCombination[i]), col + std::get<1>(currentCombination[i]));
                        }
                    }
                    (*matrix)(row - 2, col - 2) += product;
                }
            }
        });
    }
    imagesProcessedSoFar += lastImageToProcess - lastImageWithOutput;
    
    // normalize by number of images
    tbb::parallel_do(pixelMap.begin(), pixelMap.end(), [=](std::pair<PixelCombination, ImagePtr> item) {
        ImagePtr matrix = item.second;
        (*matrix) /= static_cast<double>(nImagesIncluded);
    });
    if (wantAverageImage) {
        *averageImage /= static_cast<double>(nImagesIncluded);
    }
    
    // and make the SOFI images
    sofiImages.resize(orders.size());
    tbb::parallel_for<int>(0, orders.size(), [&](int i) {
        std::vector<std::vector<double>> combinationWeights;
        const std::pair<int, std::vector<SOFIKernel> >& calculation = orders[i];
        int order = calculation.first;
        const std::vector<SOFIKernel>& kernels = calculation.second;
        ImagePtr sofiImage = AssembleSOFIImage(nRows, nCols, order, isAuto, kernels, pixelMap, pixelCombinationWeights, combinationWeights);
        sofiImages[i] = sofiImage;
        if (wantDebugMessages) {
            PrintVirtualPixelInfo(kernels, combinationWeights);
        }
    });
    
    // make JackKnife images if needed
    if (wantJackKnife) {
        jackKnifeImages.resize(orders.size());
        for (size_t i = 0; i < orders.size(); ++i) {
            const std::pair<int, std::vector<SOFIKernel> >& calculation = orders[i];
            int order = calculation.first;
            const std::vector<SOFIKernel>& kernels = calculation.second;
            std::vector<ImagePtr> jackKnifeImagesForThisOrder;
            JackKnife(imageLoader, firstImageToProcess, lastImageToProcess, order, timeLags, isAuto, kernels, pixelMap, pixelCombinationWeights, jackKnifeImagesForThisOrder);
            jackKnifeImages[i] = jackKnifeImagesForThisOrder;
        }
    }
    
    return nImagesIncluded;
}

ImagePtr AssembleSOFIImage(const int nInputRows, const int nInputCols, const int order, const bool isAuto, const std::vector<SOFIKernel>& kernels, const std::map<PixelCombination,ImagePtr,ComparePixelCombinations>& pixelMap, const std::vector<double>& pixelCombinationWeights, std::vector<std::vector<double>>& usedCombinationWeights) {
    int nVirtualPixelsPerPoint1D = (isAuto) ? 1 : order;
    ImagePtr sofiImage(new Image(nVirtualPixelsPerPoint1D * (nInputRows - 4), nVirtualPixelsPerPoint1D * (nInputCols - 4)));
    usedCombinationWeights.clear();
    usedCombinationWeights.resize(kernels.size());
    
    // reformat the explicit pixel combination weights to make it easier to apply them in parallel
    std::vector<std::vector<double> > explicitCombinationWeights;
    if (!pixelCombinationWeights.empty()) {
        explicitCombinationWeights.resize(kernels.size());
        int offset = 0;
        for (size_t kernelIndex = 0; kernelIndex < kernels.size(); kernelIndex++) {
            const std::vector<PixelCombination>& pixelCombinations = kernels.at(kernelIndex).pixelCombinations;
            std::vector<double>& kernelCombinationWeights = explicitCombinationWeights.at(kernelIndex);
            for (size_t combinationIndex = 0; combinationIndex < pixelCombinations.size(); combinationIndex++) {
                kernelCombinationWeights.push_back(pixelCombinationWeights.at(offset));
                offset++;
            }
        }
    }
    
    tbb::parallel_for(size_t(0), kernels.size(), [&](const size_t index) {
        const SOFIKernel& kernel = kernels.at(index);
        Eigen::MatrixXd evaluated(nInputRows - 4, nInputCols - 4);
        evaluated.setConstant(0.0);
        double accumulatedWeights = 0.0;
        for (size_t combinationIndex = 0; combinationIndex < kernel.combinations.size(); ++combinationIndex) {
            Eigen::MatrixXd evaluatedPartitionsSet = EvaluatePartitionsSet(kernel.combinations.at(combinationIndex), pixelMap, nInputRows - 4, nInputCols - 4);
            double weight;
            if (!explicitCombinationWeights.empty()) {
                weight = explicitCombinationWeights.at(index).at(combinationIndex);
            } else {
                weight = evaluatedPartitionsSet.mean();
            }
            accumulatedWeights += std::abs(weight);
            evaluated += evaluatedPartitionsSet * weight;
            usedCombinationWeights.at(index).push_back(weight);
        }
        if (accumulatedWeights != 0.0) {
            evaluated /= accumulatedWeights;
        }
        
        for (int col = 0; col < nInputCols - 4; ++col) {
            for (int row = 0; row < nInputRows - 4; ++row) {
                int baseOutputRow = row * nVirtualPixelsPerPoint1D;
                int baseOutputCol = col * nVirtualPixelsPerPoint1D;
                (*sofiImage)(baseOutputRow + kernel.outputDeltaX, baseOutputCol + kernel.outputDeltaY) = evaluated(row, col);
            }
        }
    });
    
    return sofiImage;
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

void PerformPixelationCorrection(ImagePtr imageToCorrect, bool alsoCorrectVariance, const int order) {
    if (order == 1)
        return; // no virtual pixels for 1st order cumulant
    
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
    
    Eigen::ArrayXXd pixelVariances(order, order);
    if (alsoCorrectVariance) {
        // calculate variances
        pixelVariances.setConstant(0.0);
        for (int col = 0; col < nCols; ++col) {
            for (int row = 0; row < nRows; ++row) {
                kindOfRow = row % nKernelRows;
                kindOfCol = col % nKernelCols;
                pixelVariances(kindOfRow, kindOfCol) += square<double>((*imageToCorrect)(row, col) - pixelAverages(kindOfRow, kindOfCol));
            }
        }
        pixelVariances /= static_cast<double>(nPixelsOfEachKind - 1);
    }
    
    // now determine correction factors, arbitrary to the first element
    for (int col = 0; col < nKernelCols; ++col) {
        for (int row = 0; row < nKernelRows; ++row) {
            if (row == 0 && col == 0) {
                aFactor(row, col) = 1.0;
                bTerm(row, col) = 0.0;
                continue;
            }
            if (alsoCorrectVariance) {
                aFactor(row, col) = std::sqrt(pixelVariances(0, 0) / pixelVariances(row, col));
                bTerm(row, col) = pixelAverages(0, 0) - aFactor(row, col) * pixelAverages(row, col);
            } else {
                aFactor(row, col) = pixelAverages(0, 0) / pixelAverages(row, col);
                bTerm(row, col) = 0.0;
            }
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

void JackKnife(std::shared_ptr<ImageLoader> imageLoader, const int firstImageToInclude, const int lastImageToInclude, const int order, const std::vector<int>& timeLags, const bool isAuto, const std::vector<SOFIKernel>& kernels, std::map<PixelCombination,ImagePtr,ComparePixelCombinations>& pixelMap, const std::vector<double>& pixelCombinationWeights, std::vector<ImagePtr>& jackKnifeImages) {
    int nInputRows = imageLoader->getXSize();
    int nInputCols = imageLoader->getYSize();
    
    // Mathematically pixelMap will not be changed by this function. However, its values will be manipulated and due to round-off
    // these will inevitably be different when this function returns.
    
    // find min and max time lag
    int minTimeLag = 0, maxTimeLag = 0;
    for (size_t i = 0; i < timeLags.size(); i++) {
        minTimeLag = std::min(minTimeLag, timeLags[i]);
        maxTimeLag = std::max(maxTimeLag, timeLags[i]);
    }
    // required: lowest time lag is zero, all others are zero or postive
    if (minTimeLag != 0)
        throw std::runtime_error("at least one time lag must be zero, and all time lags must be zero or positive");
    bool haveNonZeroTimeLag = ((minTimeLag != 0) || (maxTimeLag != 0));
    
    int nImagesToInclude = lastImageToInclude - maxTimeLag - firstImageToInclude + 1;
    jackKnifeImages.clear();
    jackKnifeImages.reserve(nImagesToInclude);
    
    // pixelMap will have been normalized by the number of images
    tbb::parallel_do(pixelMap.begin(), pixelMap.end(), [=](std::pair<PixelCombination, ImagePtr> item) {
        ImagePtr matrix = item.second;
        (*matrix) *= static_cast<double>(nImagesToInclude);
    });
    
    int nImagesInBuffer = maxTimeLag - minTimeLag + 1;
    int firstImageToProcess = firstImageToInclude;
    int lastImageToProcess = lastImageToInclude - maxTimeLag;
    
    imageLoader->spoolTo(firstImageToInclude);
    boost::circular_buffer<ImagePtr> imageBuffer(nImagesInBuffer);
    for (int i = firstImageToProcess; (i < nImagesInBuffer - 1) && (i <= lastImageToProcess); i++) {
        imageBuffer.push_back(imageLoader->readNextImage());
    }
    
    for (int i = firstImageToProcess; i <= lastImageToProcess; ++i) {
        ImagePtr currentImage = imageLoader->readNextImage();
        imageBuffer.push_back(currentImage);
        
        // subtract the contribution of the current image from the pixelMap and normalize
        tbb::parallel_do(pixelMap.begin(), pixelMap.end(), [=](std::pair<PixelCombination,ImagePtr> item) {
            const PixelCombination& currentCombination = item.first;
            ImagePtr matrix = item.second;
            int nRowsToCalculate = matrix->rows();
            int nColsToCalculate = matrix->cols();
            for (int col = 2; col < nColsToCalculate; ++col) {
                for (int row = 2; row < nRowsToCalculate; ++row) {
                    double product = 1.0;
                    if (haveNonZeroTimeLag) {
                        for (size_t i = 0; i < currentCombination.size(); ++i) {
                            const ImagePtr& imageAtTimeLag = imageBuffer.at(std::get<2>(currentCombination[i]));
                            product *= (*imageAtTimeLag)(row + std::get<0>(currentCombination[i]), col + std::get<1>(currentCombination[i]));
                        }
                    } else {
                        for (size_t i = 0; i < currentCombination.size(); ++i) {
                            product *= (*currentImage)(row + std::get<0>(currentCombination[i]), col + std::get<1>(currentCombination[i]));
                        }
                    }
                    (*matrix)(row - 2, col - 2) -= product;
                }
            }
            
            // normalize the pixelMap
            (*matrix) /= static_cast<double>(nImagesToInclude - 1);
        });
        
        // calculate a SOFI image without the current image's contribution
        std::vector<std::vector<double>> usedCombinationWeights;
        ImagePtr partialSOFI = AssembleSOFIImage(nInputRows, nInputCols, order, isAuto, kernels, pixelMap, pixelCombinationWeights, usedCombinationWeights);
        
        // undo the normalization and add the contribution of the current image back in
        tbb::parallel_do(pixelMap.begin(), pixelMap.end(), [=](std::pair<PixelCombination,ImagePtr> item) {
            const PixelCombination& currentCombination = item.first;
            ImagePtr matrix = item.second;
            
            // undo normalization
            (*matrix) *= static_cast<double>(nImagesToInclude - 1);
            
            // add the contribution of the current image back in
            int nRowsToCalculate = matrix->rows();
            int nColsToCalculate = matrix->cols();
            for (int col = 2; col < nColsToCalculate; ++col) {
                for (int row = 2; row < nRowsToCalculate; ++row) {
                    double product = 1.0;
                    if (haveNonZeroTimeLag) {
                        for (size_t i = 0; i < currentCombination.size(); ++i) {
                            const ImagePtr& imageAtTimeLag = imageBuffer.at(std::get<2>(currentCombination[i]));
                            product *= (*imageAtTimeLag)(row + std::get<0>(currentCombination[i]), col + std::get<1>(currentCombination[i]));
                        }
                    } else {
                        for (size_t i = 0; i < currentCombination.size(); ++i) {
                            product *= (*currentImage)(row + std::get<0>(currentCombination[i]), col + std::get<1>(currentCombination[i]));
                        }
                    }
                    (*matrix)(row - 2, col - 2) += product;
                }
            }
        });
        
        // store the partial SOFI image
        jackKnifeImages.push_back(partialSOFI);
    }
    
    // re-normalize the pixel map
    tbb::parallel_do(pixelMap.begin(), pixelMap.end(), [=](std::pair<PixelCombination, ImagePtr> item) {
        ImagePtr matrix = item.second;
        (*matrix) /= static_cast<double>(nImagesToInclude);
    });
}

int NumberOfPixelCombinationsInKernels(const std::vector<SOFIKernel>& kernels) {
    int nCombinations = 0;
    for (auto it = kernels.cbegin(); it != kernels.cend(); it++) {
        nCombinations += it->pixelCombinations.size();
    }
    
    return nCombinations;
}

void PrintVirtualPixelInfo(const std::vector<SOFIKernel>& kernels, const std::vector<std::vector<double>>& combinationWeights) {
    std::ostringstream ss;
    
    for (size_t i = 0; i < kernels.size(); i++) {
        const SOFIKernel& singleKernel = kernels[i];
        ss << "(" << singleKernel.outputDeltaX << "," << singleKernel.outputDeltaY << ")\r";
        for (size_t pixelIndex = 0; pixelIndex < singleKernel.pixelCombinations.size(); pixelIndex++) {
            ss << "\t" << singleKernel.scores[pixelIndex] << "\t" << combinationWeights.at(i).at(pixelIndex) << "\t";
            ss << "{";
            for (auto pixelPairIt = singleKernel.pixelCombinations[pixelIndex].cbegin(); pixelPairIt != singleKernel.pixelCombinations[pixelIndex].cend(); pixelPairIt++) {
                ss << "(" << std::get<0>(*pixelPairIt) << "," << std::get<1>(*pixelPairIt) << "," << std::get<2>(*pixelPairIt) << ")";
            }
            ss << "}\r";
        }
    }
    
#ifdef WITH_IGOR
    PrintToHistory(ss.str());
#endif
}

Eigen::MatrixXd PixelCombinationsForOrderAsMatrix(int order, double pixelCombinationCutoff) {
    std::vector<int> dummyLagTimes(order, 0);
    std::vector<SOFIKernel> kernels = KernelsForOrder(order, dummyLagTimes, pixelCombinationCutoff);
    int nVirtualPixels = kernels.size();
    
    int nCombinations = NumberOfPixelCombinationsInKernels(kernels);
    
    Eigen::MatrixXd combinations(nCombinations, 2 * order + 2);
    combinations.setConstant(99);
    int offset = 0;
    for (int kernelIndex = 0; kernelIndex < nVirtualPixels; kernelIndex++) {
        const SOFIKernel& kernel = kernels.at(kernelIndex);
        for (size_t combinationIndex = 0; combinationIndex < kernel.pixelCombinations.size(); combinationIndex++) {
            const PixelCombination& pixelCombination = kernel.pixelCombinations.at(combinationIndex);
            combinations(offset, 0) = kernel.outputDeltaX;
            combinations(offset, 1) = kernel.outputDeltaY;
            for (size_t pixelsInCombinationIndex = 0; pixelsInCombinationIndex < pixelCombination.size(); pixelsInCombinationIndex++) {
                combinations(offset, 2 + 2 * pixelsInCombinationIndex) = std::get<0>(pixelCombination[pixelsInCombinationIndex]);
                combinations(offset, 2 + 2 * pixelsInCombinationIndex + 1) = std::get<1>(pixelCombination[pixelsInCombinationIndex]);
            }
            offset++;
        }
    }
    
    return combinations;
}
