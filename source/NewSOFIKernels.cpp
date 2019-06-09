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

std::vector<SOFIVirtualPixel> SOFIVirtualPixelsForOrder(int order, const std::vector<int>& timeLags, const SOFIOptions::AllowablePixelCombinations allowablePixelCombinations, double pixelCombinationCutoff);
std::vector<SOFIVirtualPixel> SOFIAutoCumulantPixelsForOrder(int order, const std::vector<int>& timeLags);
std::string SummarizeSOFIVirtualPixels(const std::vector<SOFIVirtualPixel>& virtualPixels);

std::vector<SOFIKernel> KernelsForOrder(const int order, const std::vector<int>& timeLags, const SOFIOptions::AllowablePixelCombinations allowablePixelCombinations, const double pixelCombinationCutoff) {
    std::vector<SOFIKernel> kernels = SOFIVirtualPixelsToKernels(SOFIVirtualPixelsForOrder(order, timeLags, allowablePixelCombinations, pixelCombinationCutoff));
    
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

std::vector<SOFIKernel> AutoKernelsForOrder(const int order, const std::vector<int>& timeLags) {
    return SOFIVirtualPixelsToKernels(SOFIAutoCumulantPixelsForOrder(order, timeLags));
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
    kernel.pixelCombinations = virtualPixel.combinations;
    kernel.scores = virtualPixel.scores;
    for (auto it = virtualPixel.combinations.cbegin(); it != virtualPixel.combinations.cend(); ++it) {
        kernel.combinations.push_back(AllPartitions(*it));
    }
    return kernel;
}

void SortPixelCombination(PixelCombination& combination) {
    std::sort(combination.begin(), combination.end());
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
        pixelCombination.push_back(PixelCoordinate(i, i, 0));
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
            ss << "(" << pixelIt->dx << "," << pixelIt->dy << ", " << pixelIt->dt << ")";
        }
        ss << "} ";
    }
    ss << "\r";
    return ss.str();
}

std::string PrintSOFIVirtualPixels(const std::vector<SOFIVirtualPixel>& virtualPixels) {
    std::ostringstream ss;
    
    for (size_t i = 0; i < virtualPixels.size(); i++) {
        const SOFIVirtualPixel& singleVirtualPixel = virtualPixels[i];
        ss << "(" << singleVirtualPixel.outputDeltaX << "," << singleVirtualPixel.outputDeltaY << ")\r";
        for (size_t pixelIndex = 0; pixelIndex < singleVirtualPixel.combinations.size(); pixelIndex++) {
            ss << "\t" << singleVirtualPixel.scores[pixelIndex] << "\t";
            ss << "{";
            for (auto pixelPairIt = singleVirtualPixel.combinations[pixelIndex].cbegin(); pixelPairIt != singleVirtualPixel.combinations[pixelIndex].cend(); pixelPairIt++) {
                ss << "(" << pixelPairIt->dx << "," << pixelPairIt->dy << ", " << pixelPairIt->dt << ")";
            }
            ss << "}\r";
        }
    }
    
    return ss.str();
}

class PixelCombinationIterator {
public:
    PixelCombinationIterator(const std::vector<int>& lagTimes) :
        _lagTimes(lagTimes) {}
    virtual ~PixelCombinationIterator() {;}
    
    void nextPixelCombination(std::vector<PixelCoordinate>& pixelCombination);
    virtual bool exhaustedAllCombinations() const = 0;
private:
    virtual void _derivedNextPixelCombination(std::vector<PixelCoordinate>& pixelCombination) = 0;
    
    std::vector<int> _lagTimes;
};

void PixelCombinationIterator::nextPixelCombination(std::vector<PixelCoordinate>& pixelCombination) {
    if (_lagTimes.size() < pixelCombination.size())
        throw std::logic_error("not enough lag times");
    
    _derivedNextPixelCombination(pixelCombination);
    
    for (size_t i = 0; i < pixelCombination.size(); ++i) {
        pixelCombination[i].dt = _lagTimes[i];
    }
}

class PixelCombinationIterator_NoRepeats : public PixelCombinationIterator {
public:
    PixelCombinationIterator_NoRepeats(const int order, const std::vector<int>&  lagTimes);
    ~PixelCombinationIterator_NoRepeats() override {;}
    
    bool exhaustedAllCombinations() const override;
private:
    void _derivedNextPixelCombination(std::vector<PixelCoordinate>& pixelCombination) override;
    
    int _order;
    std::vector<PixelCoordinate> _inputPixels;
    unsigned int _nextPermutation;
    unsigned int _lastPermutation;
};

PixelCombinationIterator_NoRepeats::PixelCombinationIterator_NoRepeats(const int order, const std::vector<int>&  lagTimes) :
    PixelCombinationIterator(lagTimes),
    _order(order),
    _nextPermutation(0),
    _lastPermutation(0)
{
    if (!Within(order, kMinSofiOrder, kMaxSofiOrder))
        throw std::logic_error("unsupported order");
    
    int neighborhoodSize = 5;
    for (int i = -neighborhoodSize / 2; i <= neighborhoodSize / 2; i++) {
        for (int j = -neighborhoodSize / 2; j <= neighborhoodSize / 2; j++) {
            _inputPixels.push_back(PixelCoordinate(i, j, 0));
        }
    }
    
    _nextPermutation = (1 << order) - 1; // starting permutation of bits
    _lastPermutation = 1 << _inputPixels.size();
}

void PixelCombinationIterator_NoRepeats::_derivedNextPixelCombination(std::vector<PixelCoordinate>& pixelCombination) {
    int offset = 0;
    for (size_t j = 0; j < _inputPixels.size(); j++) {
        if ((_nextPermutation >> j) & 1) {
            pixelCombination.at(offset) = _inputPixels.at(j);
            offset++;
        }
    }

    // http://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
    unsigned int t = (_nextPermutation | (_nextPermutation - 1)) + 1;
    unsigned int w = t | ((((t & -t) / (_nextPermutation & -_nextPermutation)) >> 1) - 1);  // next permutation of bits
    _nextPermutation = w;
}

bool PixelCombinationIterator_NoRepeats::exhaustedAllCombinations() const {
    return (_nextPermutation >= _lastPermutation);
}

class PixelCombinationIterator_Repeats : public PixelCombinationIterator {
public:
    PixelCombinationIterator_Repeats(const int order, const std::vector<int>&  lagTimes);
    ~PixelCombinationIterator_Repeats() override {;}
    
    bool exhaustedAllCombinations() const override;
private:
    void _derivedNextPixelCombination(std::vector<PixelCoordinate>& pixelCombination) override;
    
    int _order;
    std::vector<PixelCoordinate> _inputPixels;
    unsigned int _nCoordinatesInNeighborhood;
    unsigned int _currentIndices;
    unsigned int _lastIndices;
};

PixelCombinationIterator_Repeats::PixelCombinationIterator_Repeats(const int order, const std::vector<int>&  lagTimes) :
    PixelCombinationIterator(lagTimes),
    _order(order),
    _currentIndices(0),
    _lastIndices(0)
{
    if (!Within(order, kMinSofiOrder, kMaxSofiOrder))
        throw std::logic_error("unsupported order");
    
    int neighborhoodSize = 5;
    for (int i = -neighborhoodSize / 2; i <= neighborhoodSize / 2; i++) {
        for (int j = -neighborhoodSize / 2; j <= neighborhoodSize / 2; j++) {
            _inputPixels.push_back(PixelCoordinate(i, j, 0));
        }
    }
    
    _nCoordinatesInNeighborhood = _inputPixels.size();
    
    
    _lastIndices = 1;
    for (int i = 0; i < order ;++i) {
        _lastIndices *= _nCoordinatesInNeighborhood;
    }
    _lastIndices -= 1;
}

void PixelCombinationIterator_Repeats::_derivedNextPixelCombination(std::vector<PixelCoordinate>& pixelCombination) {
    int currentIndices = _currentIndices;
    
    for (int i = 0; i < _order; ++i) {
        int pixelIndex = (currentIndices % _nCoordinatesInNeighborhood);
        currentIndices /= _nCoordinatesInNeighborhood;
        pixelCombination.at(i) = _inputPixels.at(pixelIndex);
    }
    
    _currentIndices += 1;
}

bool PixelCombinationIterator_Repeats::exhaustedAllCombinations() const {
    return (_currentIndices > _lastIndices);
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
    
    void addCombination(const std::vector<PixelCoordinate>& combination);
    SOFIVirtualPixel getCombination(double pixelCombinationCutoff) const;
    int getOutputDeltaX() const {return _outputDeltaX;}
    int getOutputDeltaY() const {return _outputDeltaY;}
    
private:
    int _outputDeltaX, _outputDeltaY;
    int _maxNCombinations;
    
    double _worstScore;
    int _worstScoreIndex;
    
    std::vector<std::vector<PixelCoordinate>> _pixelCombinations;
    std::vector<double> _scores;
    
};

void PixelCombinationAccumulator::addCombination(const std::vector<PixelCoordinate>& combination) {
    // calculate product of all pairwise distances in the combination. This will be the metric for determining
    // which pixel combinations are good. The lower the distance, the better.
    double score = 1.0;
    for (size_t j = 0; j < combination.size(); j++) {
        for (size_t i = j + 1; i < combination.size(); i++) {
            double sqDistance = square((double)combination[i].dx - (double)combination[j].dx) + square((double)combination[i].dy - (double)combination[j].dy);
            sqDistance = std::max(sqDistance, 0.2);
            score *= sqDistance;
        }
    }
    score = std::sqrt(score);
    
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
            _scores[_worstScoreIndex] = score;
            auto it = std::max_element(_scores.cbegin(), _scores.cend());
            _worstScore = *it;
            _worstScoreIndex = it - _scores.cbegin();
        }
    }
}

SOFIVirtualPixel PixelCombinationAccumulator::getCombination(double pixelCombinationCutoff) const {
    std::vector<int> sortIndex(_pixelCombinations.size());
    for (int i = 0; i < sortIndex.size(); i++)
        sortIndex[i] = i;
    
    std::sort(sortIndex.begin(), sortIndex.end(), [=](int i, int j) -> bool {return _scores[i] < _scores[j];});
    
    SOFIVirtualPixel result;
    result.outputDeltaX = _outputDeltaX;
    result.outputDeltaY = _outputDeltaY;
    double bestCombinationScore = _scores.at(sortIndex.at(0));
    double scoreLimit = pixelCombinationCutoff * bestCombinationScore;
    for (int i = 0; i < _pixelCombinations.size(); i++) {
        if (_scores.at(sortIndex.at(i)) > scoreLimit)
            break;
        result.combinations.push_back(_pixelCombinations.at(sortIndex.at(i)));
        result.scores.push_back(_scores.at(sortIndex.at(i)));
    }
    
    return result;
}

std::vector<SOFIVirtualPixel> SOFIVirtualPixelsForOrder(int order, const std::vector<int>& timeLags, const SOFIOptions::AllowablePixelCombinations allowablePixelCombinations, double pixelCombinationCutoff) {
    if (Clip(order, kMinSofiOrder, kMaxSofiOrder) != order)
        throw std::logic_error("unsupported order");
    if (timeLags.size() < order)
        throw std::logic_error("too few or too many time lags");
    
    int maxNCombinations = 32;
    int neighborhoodSize = 5;
    
    int nSubPixels = order * order;
    std::vector<PixelCombinationAccumulator> accumulators;
    for (int i = 0; i < nSubPixels; i++) {
        accumulators.push_back(PixelCombinationAccumulator(maxNCombinations, i / order, i % order));
    }
    
    std::shared_ptr<PixelCombinationIterator> pixelCombinationIterator;
    switch (allowablePixelCombinations) {
        case SOFIOptions::NonOverlappingPixels:
            pixelCombinationIterator = std::shared_ptr<PixelCombinationIterator>(new PixelCombinationIterator_NoRepeats(order, timeLags));
            break;
        case SOFIOptions::AllowOverlappingPixels:
            pixelCombinationIterator = std::shared_ptr<PixelCombinationIterator>(new PixelCombinationIterator_Repeats(order, timeLags));
            break;
        default:
            throw std::logic_error("unknown allowable pixel combinations");
            break;
    }
    
    std::vector<PixelCoordinate> pixelCombination(order);
    while (!pixelCombinationIterator->exhaustedAllCombinations()) {
        pixelCombinationIterator->nextPixelCombination(pixelCombination);
        
        int outputDeltaX = 0, outputDeltaY = 0;
        for (size_t j = 0; j < pixelCombination.size(); j++) {
            outputDeltaX += pixelCombination[j].dx;
            outputDeltaY += pixelCombination[j].dy;
        }
        
        if (!Within(outputDeltaX, 0, order - 1) || !Within(outputDeltaY, 0, order - 1))
            continue;
        
        for (size_t j = 0; j < pixelCombination.size(); j++) {
            pixelCombination[j].dt = timeLags[j];
        }
        
        assert(accumulators.at(outputDeltaX * order + outputDeltaY).getOutputDeltaX() == outputDeltaX);
        assert(accumulators.at(outputDeltaX * order + outputDeltaY).getOutputDeltaY() == outputDeltaY);
        accumulators.at(outputDeltaX * order + outputDeltaY).addCombination(pixelCombination);
    }
    
    std::vector<SOFIVirtualPixel> sofiPixelCombinations(accumulators.size());
    for (size_t i = 0; i < accumulators.size(); i++) {
        sofiPixelCombinations[i] = accumulators[i].getCombination(pixelCombinationCutoff);
    }
    
    return sofiPixelCombinations;
}

std::vector<SOFIVirtualPixel> SOFIAutoCumulantPixelsForOrder(int order, const std::vector<int>& timeLags) {
    if (timeLags.size() < order)
        throw std::logic_error("too few or too many time lags");
    std::vector<SOFIVirtualPixel> result(1);
    SOFIVirtualPixel& pixel = result.at(0);
    pixel.outputDeltaX = 0;
    pixel.outputDeltaY = 0;
    pixel.scores.resize(1, 0.0);
    pixel.combinations.resize(1);
    for (int i = 0; i < order; i+=1) {
        pixel.combinations.at(0).push_back(PixelCoordinate(0,0,timeLags[i]));
    }
    return result;
}

std::vector<SOFIVirtualPixel> VirtualPixelsFromExternalDescription(const Eigen::MatrixXd& pixelDescription, const std::vector<int>& timeLags) {
    int nCombinations = pixelDescription.rows();
    int nCols = pixelDescription.cols();
    int nPixelsInCombinations = nCols / 2 - 1;
    if (nCols < (2 + 2 * kMaxSofiOrder)) {
        throw std::runtime_error("pixel combinations specificiation must have more columns");
    }
    if (nCols > (2 + 2 * kMaxSofiOrder)) {
        throw std::runtime_error("pixel combinations specificiation must have fewer columns");
    }
    if ((nCols % 2) != 0) {
        throw std::runtime_error("pixel combinations specificiation must have an even number of columns");
    }
    if (nPixelsInCombinations != timeLags.size()) {
        throw std::runtime_error("number of time lags must match number of combinations");
    }
    
    // first two columns are the offsets for the virtual pixels. They must all be positive.
    for (int i = 0; i < nCombinations; i++) {
        if ((static_cast<int>(pixelDescription(i, 0)) < 0.0) || (static_cast<int>(pixelDescription(i, 1)) < 0)) {
            throw std::runtime_error("virtual pixel offsets must be non-negative");
        }
    }
    
    std::vector<std::pair<int, int>> knownCombinations;
    std::vector<SOFIVirtualPixel> virtualPixels;
    std::vector<PixelCoordinate> thisCombination(nCombinations);
    for (int i = 0; i < nCombinations; i++) {
        for (int combIndex = 0; combIndex < nCombinations; combIndex++) {
            thisCombination.at(combIndex) = PixelCoordinate(pixelDescription(i, 2 + 2 * combIndex), pixelDescription(i, 2 + 2 * combIndex + 1), timeLags.at(i));
        }
        std::pair<int, int> thisOffset(pixelDescription(i, 0), pixelDescription(i, 1));
        auto findLoc = std::find(knownCombinations.cbegin(), knownCombinations.cend(), thisOffset);
        if (findLoc != knownCombinations.cend()) {
            // we have already seen this virtual pixel
            int combIndex = findLoc - knownCombinations.cbegin();
            SOFIVirtualPixel& thisPixel = virtualPixels.at(combIndex);
            thisPixel.combinations.push_back(thisCombination);
        } else {
            // new virtual pixel
            SOFIVirtualPixel thisPixel;
            thisPixel.outputDeltaX = static_cast<int>(pixelDescription(i, 0));
            thisPixel.outputDeltaY = static_cast<int>(pixelDescription(i, 1));
            thisPixel.combinations.push_back(thisCombination);
            knownCombinations.push_back(thisOffset);
            virtualPixels.push_back(thisPixel);
        }
    }
    
    for (SOFIVirtualPixel& pixel : virtualPixels) {
        pixel.scores.resize(pixel.combinations.size(), 1.0);
    }
    
    return virtualPixels;
}

std::vector<SOFIKernel> KernelsFromExternalDescription(const Eigen::MatrixXd& kernelDescription, const std::vector<int>& timeLags) {
    return (SOFIVirtualPixelsToKernels(VirtualPixelsFromExternalDescription(kernelDescription, timeLags)));
}
