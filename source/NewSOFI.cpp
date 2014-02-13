
#include "NewSOFI.h"

#include <map>
#include <sstream>

#include <eigen3/Eigen/Eigen>
#include "tbb/tbb.h"
#include "tbb/spin_mutex.h"
#include "XOPStandardHeaders.h"

#include "NewSOFIKernels.h"

double Prefactor(int nPartitions);
Eigen::MatrixXd EvaluatePartition(const Partition& partition, const std::map<PixelCombination,ImagePtr,ComparePixelCombinations>& pixelMap);
Eigen::MatrixXd EvaluatePartitionsSet(const GroupOfPartitions& groupOfPartitions, const std::map<PixelCombination,ImagePtr,ComparePixelCombinations>& pixelMap);

void DoNewSOFI(std::shared_ptr<ImageLoader> imageLoader, std::shared_ptr<ProgressReporter> progressReporter, const int order, std::vector<ImagePtr>& sofiOutputImages) {
    int nRows = imageLoader->getXSize();
    int nCols = imageLoader->getYSize();
    int nImages = imageLoader->getNImages();
    
    progressReporter->CalculationStarted();
    
    std::vector<SOFIKernel> kernels;
    if (order == 2) {
        kernels = PixelCombinationsToKernels(sofiPixelCombinations2());
    } else if (order == 3) {
        kernels = PixelCombinationsToKernels(sofiPixelCombinations3());
    } else if (order == 4) {
        kernels = PixelCombinationsToKernels(sofiPixelCombinations4());
    } else if (order == 5) {
        kernels = PixelCombinationsToKernels(sofiPixelCombinations5());
    } else if (order == 6) {
        kernels = PixelCombinationsToKernels(sofiPixelCombinations6());
    } else {
        throw std::runtime_error("unsupported order");
    }
    
    std::map<PixelCombination, ImagePtr , ComparePixelCombinations> pixelMap;
    
    // store all needed pixel combinations in the map
    for (auto kernelIt = kernels.cbegin(); kernelIt != kernels.cend(); ++kernelIt) {
        const std::vector<GroupOfPartitions>& GroupOfPartitions = kernelIt->combinations;
        for (auto partitionsSetIt = GroupOfPartitions.cbegin(); partitionsSetIt != GroupOfPartitions.cend(); ++partitionsSetIt) {
            for (auto partitionsIt = partitionsSetIt->cbegin(); partitionsIt != partitionsSetIt->cend(); ++partitionsIt) {
                for (auto subsetIt = partitionsIt->cbegin(); subsetIt != partitionsIt->cend(); ++subsetIt) {
                    ImagePtr matrix(new Image(nRows - 4, nCols - 4));
                    matrix->setConstant(0.0);
                    pixelMap.insert(std::pair<PixelCombination, ImagePtr>(*subsetIt, matrix));
                }
            }
        }
    }
    
    std::ostringstream ostream;
    ostream << "Map has total of " << pixelMap.size() << " entries\r";
    XOPNotice(ostream.str().c_str());
    
    // calculate all products over the images
    for (int n = 0; n < nImages; ++n) {
        //tbb::parallel_for(0, nImages, [&](const int n) {
        progressReporter->UpdateCalculationProgress(n, nImages);
        ImagePtr currentImage = imageLoader->readImage(n);
        tbb::parallel_do(pixelMap.begin(), pixelMap.end(), [=](std::pair<PixelCombination,ImagePtr> item) {
            const PixelCombination& currentCombination = item.first;
            ImagePtr matrix = item.second;
            for (int col = 2; col < nCols - 2; ++col) {
                for (int row = 2; row < nRows - 2; ++row) {
                    double product = 1.0;
                    for (int i = 0; i < currentCombination.size(); ++i) {
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
    
    sofiOutputImages.clear();
    sofiOutputImages.push_back(sofiImage);
    
    progressReporter->CalculationDone();
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
    for (int i = 0; i < partition.size(); ++i) {
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
