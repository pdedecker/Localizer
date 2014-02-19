
#include "NewSOFI.h"

#include <map>
#include <sstream>

#include <eigen3/Eigen/Eigen>
#include "tbb/tbb.h"
#include "tbb/spin_mutex.h"
#include "XOPStandardHeaders.h"

double Prefactor(int nPartitions);
Eigen::MatrixXd EvaluatePartition(const Partition& partition, const std::map<PixelCombination,ImagePtr,ComparePixelCombinations>& pixelMap);
Eigen::MatrixXd EvaluatePartitionsSet(const GroupOfPartitions& groupOfPartitions, const std::map<PixelCombination,ImagePtr,ComparePixelCombinations>& pixelMap);

void PerformPixelationCorrection(ImagePtr imageToCorrect, const int order);

void DoNewSOFI(std::shared_ptr<ImageLoader> imageLoader, std::shared_ptr<ProgressReporter> progressReporter, const int order, std::vector<ImagePtr>& sofiOutputImages) {
    int nRows = imageLoader->getXSize();
    int nCols = imageLoader->getYSize();
    int nImages = imageLoader->getNImages();
    
    progressReporter->CalculationStarted();
    
    std::vector<SOFIKernel> kernels = KernelsForOrder(order);
    std::map<PixelCombination, ImagePtr , ComparePixelCombinations> pixelMap;
    
    // store all needed pixel combinations in the map
    int allCombinations = 0;
    for (auto kernelIt = kernels.cbegin(); kernelIt != kernels.cend(); ++kernelIt) {
        const std::vector<GroupOfPartitions>& GroupOfPartitions = kernelIt->combinations;
        for (auto partitionsSetIt = GroupOfPartitions.cbegin(); partitionsSetIt != GroupOfPartitions.cend(); ++partitionsSetIt) {
            for (auto partitionsIt = partitionsSetIt->cbegin(); partitionsIt != partitionsSetIt->cend(); ++partitionsIt) {
                for (auto subsetIt = partitionsIt->cbegin(); subsetIt != partitionsIt->cend(); ++subsetIt) {
                    allCombinations += 1;
                    if (!pixelMap.count(*subsetIt)) {
                        ImagePtr matrix(new Image(nRows - 4, nCols - 4));
                        matrix->setConstant(0.0);
                        pixelMap.insert(std::pair<PixelCombination, ImagePtr>(*subsetIt, matrix));
                    }
                }
            }
        }
    }
    
    std::ostringstream ostream;
    ostream << "Map has total of " << pixelMap.size() << " entries (" << static_cast<int>(static_cast<double>(pixelMap.size()) / static_cast<double>(allCombinations) * 100.0) << "% retained)\r";
    XOPNotice(ostream.str().c_str());
    
    // calculate all products over the images
    for (int n = 0; n < nImages; ++n) {
        //tbb::parallel_for(0, nImages, [&](const int n) {
        int abortStatus = progressReporter->UpdateCalculationProgress(n, nImages);
        if (abortStatus)
            throw USER_ABORTED("user abort");
        ImagePtr currentImage = imageLoader->readImage(n);
        tbb::parallel_do(pixelMap.begin(), pixelMap.end(), [=](std::pair<PixelCombination,ImagePtr> item) {
            const PixelCombination& currentCombination = item.first;
            ImagePtr matrix = item.second;
            for (int col = 2; col < nCols - 2; ++col) {
                for (int row = 2; row < nRows - 2; ++row) {
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
    for (auto it = pixelMap.begin(); it != pixelMap.end(); ++it) {
        *(it->second) /= static_cast<double>(nImages);
    }
    
    // and make the SOFI image
    ImagePtr sofiImage(new Image(order * (nRows - 4), order * (nCols - 4)));
    tbb::parallel_do(kernels.cbegin(), kernels.cend(), [=,&pixelMap,&sofiImage](const SOFIKernel& kernel) {
        Eigen::MatrixXd evaluated(nRows - 4, nCols - 4);
        evaluated.setConstant(0.0);
        for (auto partitionsSetIt = kernel.combinations.cbegin(); partitionsSetIt != kernel.combinations.cend(); ++partitionsSetIt) {
            evaluated += EvaluatePartitionsSet(*partitionsSetIt, pixelMap);
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

    PerformPixelationCorrection(sofiImage, order);
    sofiOutputImages.clear();
    sofiOutputImages.push_back(sofiImage);
    
    progressReporter->CalculationDone();
}

void RawSOFIWorker(std::shared_ptr<ImageLoader> imageLoader, const int firstImageToProcess, const int lastImageToProcess, int &imagesProcessedSoFar, const int& totalNumberOfImagesToProcess, std::shared_ptr<ProgressReporter> progressReporter, const std::vector<std::pair<int, std::vector<SOFIKernel> > >& orders, std::map<PixelCombination,ImagePtr,ComparePixelCombinations>& pixelMap, std::vector<ImagePtr>& sofiImages) {
    int nRows = imageLoader->getXSize();
    int nCols = imageLoader->getYSize();
    
    // store all needed pixel combinations in the map
    // if the map is non-empty, we assume that it was already setup by a previous call to RawSOFIWorker, so we can prevent unnecessary work.
    if (pixelMap.empty()) {
        int allCombinations = 0;
        for (size_t i = 0; i < orders.size(); ++i) {
            const std::pair<int, std::vector<SOFIKernel> >& calculation = orders[i];
            const std::vector<SOFIKernel>& kernels = calculation.second;
            for (auto kernelIt = kernels.cbegin(); kernelIt != kernels.cend(); ++kernelIt) {
                const std::vector<GroupOfPartitions>& GroupOfPartitions = kernelIt->combinations;
                for (auto partitionsSetIt = GroupOfPartitions.cbegin(); partitionsSetIt != GroupOfPartitions.cend(); ++partitionsSetIt) {
                    for (auto partitionsIt = partitionsSetIt->cbegin(); partitionsIt != partitionsSetIt->cend(); ++partitionsIt) {
                        for (auto subsetIt = partitionsIt->cbegin(); subsetIt != partitionsIt->cend(); ++subsetIt) {
                            allCombinations += 1;
                            if (!pixelMap.count(*subsetIt)) {
                                ImagePtr matrix(new Image(nRows - 4, nCols - 4));
                                pixelMap.insert(std::pair<PixelCombination, ImagePtr>(*subsetIt, matrix));
                            }
                        }
                    }
                }
            }
        }
        std::ostringstream ostream;
        ostream << "Map has total of " << pixelMap.size() << " entries (" << static_cast<int>(static_cast<double>(pixelMap.size()) / static_cast<double>(allCombinations) * 100.0) << "% retained)\r";
        XOPNotice(ostream.str().c_str());
    }
    
    // clear all accumulated pixel combinations
    tbb::parallel_do(pixelMap.begin(), pixelMap.end(), [](std::pair<PixelCombination, ImagePtr> item) {
        ImagePtr matrix = item.second;
        matrix->setConstant(0.0);
    });
    
    // calculate all products over the images
    for (int n = firstImageToProcess; n <= lastImageToProcess; ++n) {
        int abortStatus = progressReporter->UpdateCalculationProgress(imagesProcessedSoFar, totalNumberOfImagesToProcess);
        if (abortStatus)
            throw USER_ABORTED("user abort");
        ImagePtr currentImage = imageLoader->readImage(n);
        tbb::parallel_do(pixelMap.begin(), pixelMap.end(), [=](std::pair<PixelCombination,ImagePtr> item) {
            const PixelCombination& currentCombination = item.first;
            ImagePtr matrix = item.second;
            for (int col = 2; col < nCols - 2; ++col) {
                for (int row = 2; row < nRows - 2; ++row) {
                    double product = 1.0;
                    for (size_t i = 0; i < currentCombination.size(); ++i) {
                        product *= (*currentImage)(row + currentCombination[i].first, col + currentCombination[i].second);
                    }
                    (*matrix)(row - 2, col - 2) += product;
                }
            }
        });
        imagesProcessedSoFar += 1;
    }
    
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
                evaluated += EvaluatePartitionsSet(*partitionsSetIt, pixelMap);
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
}

Eigen::MatrixXd EvaluatePartitionsSet(const GroupOfPartitions& groupOfPartitions, const std::map<PixelCombination,ImagePtr,ComparePixelCombinations>& pixelMap) {
    int nRows = pixelMap.cbegin()->second->rows();
    int nCols = pixelMap.cbegin()->second->cols();
    Eigen::MatrixXd accumulated(nRows, nCols);
    accumulated.setConstant(0.0);
    for (auto it = groupOfPartitions.cbegin(); it != groupOfPartitions.cend(); ++it) {
        int nSetsInPartition = it->size();
        accumulated += Prefactor(nSetsInPartition) * EvaluatePartition(*it, pixelMap);
    }
    
    return accumulated;
}

Eigen::MatrixXd EvaluatePartition(const Partition& partition, const std::map<PixelCombination,ImagePtr,ComparePixelCombinations>& pixelMap) {
    int nRows = pixelMap.cbegin()->second->rows();
    int nCols = pixelMap.cbegin()->second->cols();
    Eigen::MatrixXd result(nRows, nCols);
    result.setConstant(1.0);
    for (size_t i = 0; i < partition.size(); ++i) {
        PixelCombination subset = partition[i];
        result = result.cwiseProduct(*(pixelMap.at(subset)));
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
