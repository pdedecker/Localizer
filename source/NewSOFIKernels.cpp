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

#include "NewSOFIKernels.h"

#include <algorithm>
#include <sstream>

#ifdef WITH_IGOR
#include "XOPStandardHeaders.h"
#endif

#include "Defines.h"

std::vector<SOFIKernel> SOFIVirtualPixelsToKernels(const std::vector<SOFIVirtualPixel>& virtualPixels);
void SortPixelCombination(PixelCombination& combination);
std::vector<SOFIVirtualPixel> sofiPixelCombinations1();
std::vector<SOFIVirtualPixel> sofiPixelCombinations2();
std::vector<SOFIVirtualPixel> sofiPixelCombinations3();
std::vector<SOFIVirtualPixel> sofiPixelCombinations4();
std::vector<SOFIVirtualPixel> sofiPixelCombinations5();
std::vector<SOFIVirtualPixel> sofiPixelCombinations6();

std::vector<SOFIVirtualPixel> SOFIVirtualPixelsForOrder(int order);
std::string PrintSOFIVirtualPixels(const std::vector<SOFIVirtualPixel>& virtualPixels);

std::vector<SOFIKernel> KernelsForOrder(const int order) {
    std::vector<SOFIKernel> kernels;
    
#ifdef WITH_IGOR
    XOPNotice(PrintSOFIVirtualPixels(SOFIVirtualPixelsForOrder(2)).c_str());
    XOPNotice("\r");
#endif
    
    switch (order) {
        case 1:
            kernels = SOFIVirtualPixelsToKernels(sofiPixelCombinations1());
            break;
        case 2:
            kernels = SOFIVirtualPixelsToKernels(sofiPixelCombinations2());
            break;
        case 3:
            kernels = SOFIVirtualPixelsToKernels(sofiPixelCombinations3());
            break;
        case 4:
            kernels = SOFIVirtualPixelsToKernels(sofiPixelCombinations4());
            break;
        case 5:
            kernels = SOFIVirtualPixelsToKernels(sofiPixelCombinations5());
            break;
        case 6:
            kernels = SOFIVirtualPixelsToKernels(sofiPixelCombinations6());
            break;
        default:
            throw std::runtime_error("unsupported order");
            break;
    }
    
    for (auto kernelIt = kernels.begin(); kernelIt != kernels.end(); ++kernelIt) {
        std::vector<GroupOfPartitions>& GroupOfPartitions = kernelIt->combinations;
        for (auto partitionsSetIt = GroupOfPartitions.begin(); partitionsSetIt != GroupOfPartitions.end(); ++partitionsSetIt) {
            for (auto partitionsIt = partitionsSetIt->begin(); partitionsIt != partitionsSetIt->end(); ++partitionsIt) {
                for (auto pixelCombinationIt = partitionsIt->begin(); pixelCombinationIt != partitionsIt->end(); ++pixelCombinationIt) {
                    SortPixelCombination(*pixelCombinationIt);
                }
            }
        }
    }
    
    return kernels;
}

SOFIKernel SOFIVirtualPixelToKernel(const SOFIVirtualPixel& virtualPixel);

std::vector<SOFIKernel> SOFIVirtualPixelsToKernels(const std::vector<SOFIVirtualPixel>& virtualPixels) {
    std::vector<SOFIKernel> kernels;
    kernels.reserve(virtualPixels.size());
    
    for (auto it = virtualPixels.cbegin(); it != virtualPixels.cend(); ++it) {
        kernels.push_back(SOFIVirtualPixelToKernel(*it));
    }
    
    return kernels;
}

SOFIKernel SOFIVirtualPixelToKernel(const SOFIVirtualPixel& virtualPixel) {
    SOFIKernel kernel;
    kernel.outputDeltaX = virtualPixel.outputDeltaX;
    kernel.outputDeltaY = virtualPixel.outputDeltaY;
    for (auto it = virtualPixel.combinations.cbegin(); it != virtualPixel.combinations.cend(); ++it) {
        kernel.combinations.push_back(AllPartitions(*it));
    }
    return kernel;
}

void SortPixelCombination(PixelCombination& combination) {
    std::sort(combination.begin(), combination.end(), [](const std::pair<int,int>& pair1, const std::pair<int,int>& pair2) -> bool {
        if (pair1.first < pair2.first) {
            return true;
        } else if (pair1.first == pair2.first) {
            return (pair1.second < pair2.second);
        }
        return false;
    });
}

bool ComparePixelCombinations::operator() (const PixelCombination& comb1, const PixelCombination& comb2) const {
    if (comb1.size() < comb2.size())
        return true;
    if (comb1.size() > comb2.size())
        return false;
    
    for (size_t i = 0; i < comb1.size(); ++i) {
        if (comb1[i] == comb2[i])
            continue;
        return comb1[i] < comb2[i];
    }
    return false;
}

// partitioning code based on http://www.cs.bgu.ac.il/~orlovm/papers/partitions.pdf
std::pair<std::vector<int>, std::vector<int> > FirstPartition(const int nElements);
std::pair<std::vector<int>, std::vector<int> > LastPartition(const int nElements);
void NextPartition(std::vector<int>& kk, std::vector<int>& MM);

GroupOfPartitions AllPartitions(PixelCombination pixelCombination) {
    int nElements = pixelCombination.size();
    GroupOfPartitions groupOfPartitions;
    
    std::vector<int> kk = FirstPartition(nElements).first;
    std::vector<int> MM = FirstPartition(nElements).second;
    std::vector<int> kkLast = LastPartition(nElements).first;
    std::vector<int> MMLast = LastPartition(nElements).second;
    
    for ( ; ; ) {
        Partition partition;
        int nSubsets = *(std::max_element(kk.begin(), kk.end())) + 1;
        for (int i = 0; i < nSubsets; ++i) {
            PixelCombination subset;
            for (int j = 0; j < nElements; ++j) {
                if (kk[j] == i)
                    subset.push_back(pixelCombination[j]);
            }
            partition.push_back(subset);
        }
        groupOfPartitions.push_back(partition);
        
        if (kk == kkLast)
            break;
        NextPartition(kk, MM);
    }
    
    return groupOfPartitions;
}

std::pair<std::vector<int>, std::vector<int> > FirstPartition(const int nElements) {
    std::vector<int> kk(nElements, 0);
    std::vector<int> MM(nElements, 0);
    return std::pair<std::vector<int>, std::vector<int> >(kk, MM);
}

std::pair<std::vector<int>, std::vector<int> > LastPartition(const int nElements) {
    std::vector<int> kk(nElements);
    std::vector<int> MM(nElements);
    for (int i = 0; i < nElements; ++i) {
        kk[i] = i;
        MM[i] = i;
    }
    return std::pair<std::vector<int>, std::vector<int> >(kk, MM);
}

void NextPartition(std::vector<int>& kk, std::vector<int>& MM) {
    int nElements = kk.size();
    
    for (int i = nElements - 1; i > 0; --i) {
        if (kk[i] <= MM[i - 1]) {
            kk[i] = kk[i] + 1;
            MM[i] = std::max(MM[i], kk[i]);
            for (int j = i + 1; j < nElements; ++j) {
                kk[j] = kk[0];
                MM[j] = MM[i];
            }
            return;
        }
    }
}

std::string ParititionTest() {
    PixelCombination pixelCombination;
    for (int i = 0; i < 3; ++i) {
        pixelCombination.push_back(std::pair<int, int>(i, i));
    }
    
    std::ostringstream ss;
    GroupOfPartitions groupOfPartitions = AllPartitions(pixelCombination);
    int nPartitions = groupOfPartitions.size();
    ss << "found " << nPartitions << " partitions\r";
    for (int i = 0; i < nPartitions; ++i) {
        ss << "Partition " << i << ":\r";
        Partition partition = groupOfPartitions[i];
        ss << PrintPartition(partition);
    }
    
    return ss.str();
}

std::string PrintPartition(const Partition& partition) {
    std::ostringstream ss;
    int nSubsets = partition.size();
    for (int j = 0; j < nSubsets; ++j) {
        ss << "{";
        for (auto pixelIt = partition[j].cbegin(); pixelIt != partition[j].cend(); ++pixelIt) {
            ss << "(" << pixelIt->first << "," << pixelIt->second << ") ";
        }
        ss << "} ";
    }
    ss << "\r";
    return ss.str();
}

std::string PrintSOFIVirtualPixels(const std::vector<SOFIVirtualPixel>& virtualPixels) {
    std::ostringstream ss;
    
    for (int i = 0; i < virtualPixels.size(); i++) {
        const SOFIVirtualPixel& singleVirtualPixel = virtualPixels[i];
        ss << "(" << singleVirtualPixel.outputDeltaX << "," << singleVirtualPixel.outputDeltaY << ")\r";
        for (int pixelIndex = 0; pixelIndex < singleVirtualPixel.combinations.size(); pixelIndex++) {
            ss << "\t" << singleVirtualPixel.scores[pixelIndex] << "\t";
            ss << "{";
            for (auto pixelPairIt = singleVirtualPixel.combinations[pixelIndex].cbegin(); pixelPairIt != singleVirtualPixel.combinations[pixelIndex].cend(); pixelPairIt++) {
                ss << "(" << pixelPairIt->first << "," << pixelPairIt->second << ")";
            }
            ss << "}\r";
        }
    }
    
    return ss.str();
}

class PixelCombinationAccumulator {
public:
    PixelCombinationAccumulator(int maxNCombinations, int outputDeltaX, int outputDeltaY) :
        _maxNCombinations(maxNCombinations),
        _outputDeltaX(outputDeltaX),
        _outputDeltaY(outputDeltaY),
        _worstScore(std::numeric_limits<double>::infinity()),
        _worstScoreIndex(-10)
    {
    }
    
    void addCombination(const std::vector<std::pair<int, int>>& combination);
    SOFIVirtualPixel getCombination() const;
    int getOutputDeltaX() const {return _outputDeltaX;}
    int getOutputDeltaY() const {return _outputDeltaY;}
    
private:
    int _outputDeltaX, _outputDeltaY;
    int _maxNCombinations;
    
    double _worstScore;
    int _worstScoreIndex;
    
    std::vector<std::vector<std::pair<int, int>>> _pixelCombinations;
    std::vector<double> _scores;
    
};

void PixelCombinationAccumulator::addCombination(const std::vector<std::pair<int, int>>& combination) {
    // calculate product of all pairwise distances in the combination. This will be the metric for determining
    // which pixel combinations are good. The lower the distance, the better.
    double score = 1.0;
    for (size_t j = 0; j < combination.size(); j++) {
        for (size_t i = j + 1; i < combination.size(); i++) {
            double sqDistance = square((double)combination[i].first - (double)combination[j].first) + square((double)combination[i].second - (double)combination[j].second);
            score *= sqDistance;
        }
    }
    
    if (_pixelCombinations.size() < _maxNCombinations) {
        // still building up the vector
        _pixelCombinations.push_back(combination);
        _scores.push_back(score);
        if (_pixelCombinations.size() == _maxNCombinations) {
            auto it = std::max_element(_scores.cbegin(), _scores.cend());
            _worstScore = *it;
            _worstScoreIndex = it - _scores.cbegin();
        }
    } else {
        if (score < _worstScore) {
            _pixelCombinations[_worstScoreIndex] = combination;
            auto it = std::max_element(_scores.cbegin(), _scores.cend());
            _worstScore = *it;
            _worstScoreIndex = it - _scores.cbegin();
        }
    }
}

SOFIVirtualPixel PixelCombinationAccumulator::getCombination() const {
    std::vector<int> sortIndex(_pixelCombinations.size());
    for (int i = 0; i < sortIndex.size(); i++)
        sortIndex[i] = i;
    
    std::sort(sortIndex.begin(), sortIndex.end(), [=](int i, int j) -> bool {return _scores[i] < _scores[j];});
    
    SOFIVirtualPixel result;
    result.outputDeltaX = _outputDeltaX;
    result.outputDeltaY = _outputDeltaY;
    for (int i = 0; i < _pixelCombinations.size(); i++) {
        result.combinations.push_back(_pixelCombinations.at(sortIndex.at(i)));
        result.scores.push_back(_scores.at(sortIndex.at(i)));
    }
    
    return result;
}

std::vector<SOFIVirtualPixel> SOFIVirtualPixelsForOrder(int order) {
    if (Clip(order, 1, 6) != order)
        throw std::logic_error("unsupported order");
    
    int maxNCombinations = 10;
    int neighborhoodSize = 5;
    
    int nSubPixels = order * order;
    std::vector<PixelCombinationAccumulator> accumulators;
    for (int i = 0; i < nSubPixels; i++) {
        accumulators.push_back(PixelCombinationAccumulator(maxNCombinations, i / order, i % order));
    }
    
    std::vector<std::pair<int, int>> inputPixels;
    for (int i = -neighborhoodSize / 2; i <= neighborhoodSize / 2; i++) {
        for (int j = -neighborhoodSize / 2; j <= neighborhoodSize / 2; j++) {
            inputPixels.push_back(std::pair<int, int>(i, j));
        }
    }
    
    std::vector<std::pair<int, int>> pixelCombination(order);
    // http://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
    unsigned int v = (1 << order) - 1; // starting permutation of bits
    unsigned int lastPermutation = 1 << inputPixels.size();
    for ( ; ; ) {
        unsigned int t = (v | (v - 1)) + 1;
        unsigned int w = t | ((((t & -t) / (v & -v)) >> 1) - 1);  // next permutation of bits
        v = w;
        if (w > lastPermutation)
            break;
        
        int offset = 0;
        for (int j = 0; j < inputPixels.size(); j++) {
            if ((w >> j) & 1) {
                pixelCombination.at(offset) = inputPixels.at(j);
                offset++;
            }
        }
        
        int outputDeltaX = 0, outputDeltaY = 0;
        for (int j = 0; j < pixelCombination.size(); j++) {
            outputDeltaX += pixelCombination[j].first;
            outputDeltaY += pixelCombination[j].second;
        }
        
        if (!Within(outputDeltaX, 0, order - 1) || !Within(outputDeltaY, 0, order - 1))
            continue;
        
        assert(accumulators.at(outputDeltaX * order + outputDeltaY).getOutputDeltaX() == outputDeltaX);
        assert(accumulators.at(outputDeltaX * order + outputDeltaY).getOutputDeltaY() == outputDeltaY);
        accumulators.at(outputDeltaX * order + outputDeltaY).addCombination(pixelCombination);
    }
    
    std::vector<SOFIVirtualPixel> sofiPixelCombinations(accumulators.size());
    for (int i = 0; i < accumulators.size(); i++) {
        sofiPixelCombinations[i] = accumulators[i].getCombination();
    }
    
    return sofiPixelCombinations;
}

std::vector<SOFIVirtualPixel> sofiPixelCombinations1()
{
    std::vector<SOFIVirtualPixel> pixelCombinations;
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 0;
        pixelCombination.outputDeltaY = 0;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    return pixelCombinations;
}

std::vector<SOFIVirtualPixel> sofiPixelCombinations2()
{
    std::vector<SOFIVirtualPixel> pixelCombinations;
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 0;
        pixelCombination.outputDeltaY = 0;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 0;
        pixelCombination.outputDeltaY = 1;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 1;
        pixelCombination.outputDeltaY = 0;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 1;
        pixelCombination.outputDeltaY = 1;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    return pixelCombinations;
}
std::vector<SOFIVirtualPixel> sofiPixelCombinations3()
{
    std::vector<SOFIVirtualPixel> pixelCombinations;
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 0;
        pixelCombination.outputDeltaY = 0;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 0;
        pixelCombination.outputDeltaY = 1;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 0;
        pixelCombination.outputDeltaY = 2;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 1;
        pixelCombination.outputDeltaY = 0;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 1;
        pixelCombination.outputDeltaY = 1;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 1;
        pixelCombination.outputDeltaY = 2;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 2;
        pixelCombination.outputDeltaY = 0;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 2;
        pixelCombination.outputDeltaY = 1;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 2;
        pixelCombination.outputDeltaY = 2;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    return pixelCombinations;
}
std::vector<SOFIVirtualPixel> sofiPixelCombinations4()
{
    std::vector<SOFIVirtualPixel> pixelCombinations;
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 0;
        pixelCombination.outputDeltaY = 0;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 0;
        pixelCombination.outputDeltaY = 1;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 0;
        pixelCombination.outputDeltaY = 2;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 0;
        pixelCombination.outputDeltaY = 3;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 1;
        pixelCombination.outputDeltaY = 0;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 1;
        pixelCombination.outputDeltaY = 1;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 1;
        pixelCombination.outputDeltaY = 2;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 1;
        pixelCombination.outputDeltaY = 3;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 2;
        pixelCombination.outputDeltaY = 0;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 2;
        pixelCombination.outputDeltaY = 1;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 2;
        pixelCombination.outputDeltaY = 2;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 2;
        pixelCombination.outputDeltaY = 3;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 3;
        pixelCombination.outputDeltaY = 0;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 3;
        pixelCombination.outputDeltaY = 1;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 3;
        pixelCombination.outputDeltaY = 2;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(2, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 3;
        pixelCombination.outputDeltaY = 3;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    return pixelCombinations;
}

std::vector<SOFIVirtualPixel> sofiPixelCombinations5()
{
    std::vector<SOFIVirtualPixel> pixelCombinations;
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 0;
        pixelCombination.outputDeltaY = 0;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 0;
        pixelCombination.outputDeltaY = 1;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 0;
        pixelCombination.outputDeltaY = 2;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 0;
        pixelCombination.outputDeltaY = 3;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 0;
        pixelCombination.outputDeltaY = 4;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 1;
        pixelCombination.outputDeltaY = 0;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 1;
        pixelCombination.outputDeltaY = 1;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 1;
        pixelCombination.outputDeltaY = 2;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 1;
        pixelCombination.outputDeltaY = 3;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 1;
        pixelCombination.outputDeltaY = 4;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 2;
        pixelCombination.outputDeltaY = 0;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 2;
        pixelCombination.outputDeltaY = 1;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 2;
        pixelCombination.outputDeltaY = 2;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 2;
        pixelCombination.outputDeltaY = 3;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 2;
        pixelCombination.outputDeltaY = 4;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 3;
        pixelCombination.outputDeltaY = 0;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 3;
        pixelCombination.outputDeltaY = 1;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 3;
        pixelCombination.outputDeltaY = 2;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 3;
        pixelCombination.outputDeltaY = 3;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 3;
        pixelCombination.outputDeltaY = 4;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 4;
        pixelCombination.outputDeltaY = 0;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 4;
        pixelCombination.outputDeltaY = 1;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, -1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(2, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 4;
        pixelCombination.outputDeltaY = 2;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 4;
        pixelCombination.outputDeltaY = 3;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 4;
        pixelCombination.outputDeltaY = 4;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 2));
            subset.push_back(std::pair<int, int>(2, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    return pixelCombinations;
}

std::vector<SOFIVirtualPixel> sofiPixelCombinations6()
{
    std::vector<SOFIVirtualPixel> pixelCombinations;
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 0;
        pixelCombination.outputDeltaY = 0;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, -1));
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, -1));
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 0;
        pixelCombination.outputDeltaY = 1;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, -1));
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 0;
        pixelCombination.outputDeltaY = 2;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, -1));
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 0;
        pixelCombination.outputDeltaY = 3;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 0;
        pixelCombination.outputDeltaY = 4;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 0;
        pixelCombination.outputDeltaY = 5;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 1;
        pixelCombination.outputDeltaY = 0;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, -1));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 1;
        pixelCombination.outputDeltaY = 1;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 1;
        pixelCombination.outputDeltaY = 2;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 1;
        pixelCombination.outputDeltaY = 3;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 1;
        pixelCombination.outputDeltaY = 4;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 1;
        pixelCombination.outputDeltaY = 5;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 2;
        pixelCombination.outputDeltaY = 0;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 2;
        pixelCombination.outputDeltaY = 1;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 2;
        pixelCombination.outputDeltaY = 2;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 2;
        pixelCombination.outputDeltaY = 3;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 2;
        pixelCombination.outputDeltaY = 4;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 2;
        pixelCombination.outputDeltaY = 5;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 3;
        pixelCombination.outputDeltaY = 0;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 3;
        pixelCombination.outputDeltaY = 1;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 3;
        pixelCombination.outputDeltaY = 2;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 3;
        pixelCombination.outputDeltaY = 3;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 3;
        pixelCombination.outputDeltaY = 4;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 3;
        pixelCombination.outputDeltaY = 5;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 4;
        pixelCombination.outputDeltaY = 0;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, -1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(2, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 4;
        pixelCombination.outputDeltaY = 1;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 4;
        pixelCombination.outputDeltaY = 2;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 4;
        pixelCombination.outputDeltaY = 3;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 4;
        pixelCombination.outputDeltaY = 4;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 4;
        pixelCombination.outputDeltaY = 5;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 5;
        pixelCombination.outputDeltaY = 0;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, -1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 5;
        pixelCombination.outputDeltaY = 1;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 5;
        pixelCombination.outputDeltaY = 2;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 5;
        pixelCombination.outputDeltaY = 3;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            subset.push_back(std::pair<int, int>(2, -1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 5;
        pixelCombination.outputDeltaY = 4;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIVirtualPixel pixelCombination;
        pixelCombination.outputDeltaX = 5;
        pixelCombination.outputDeltaY = 5;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            subset.push_back(std::pair<int, int>(2, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    return pixelCombinations;
}
