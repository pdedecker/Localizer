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

std::vector<SOFIKernel> PixelCombinationsToKernels(const std::vector<SOFIPixelCombination>& sofiPixelCombination);

std::vector<SOFIKernel> KernelsForOrder(const int order) {
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
    
    return kernels;
}

SOFIKernel PixelCombinationToKernel(const SOFIPixelCombination& sofiPixelCombination);

std::vector<SOFIKernel> PixelCombinationsToKernels(const std::vector<SOFIPixelCombination>& sofiPixelCombination) {
    std::vector<SOFIKernel> kernels;
    kernels.reserve(sofiPixelCombination.size());
    
    for (auto it = sofiPixelCombination.cbegin(); it != sofiPixelCombination.cend(); ++it) {
        kernels.push_back(PixelCombinationToKernel(*it));
    }
    
    return kernels;
}

SOFIKernel PixelCombinationToKernel(const SOFIPixelCombination& sofiPixelCombination) {
    SOFIKernel kernel;
    kernel.outputDeltaX = sofiPixelCombination.outputDeltaX;
    kernel.outputDeltaY = sofiPixelCombination.outputDeltaY;
    for (auto it = sofiPixelCombination.combinations.cbegin(); it != sofiPixelCombination.combinations.cend(); ++it) {
        kernel.combinations.push_back(AllPartitions(*it));
    }
    return kernel;
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

std::vector<SOFIPixelCombination> sofiPixelCombinations2()
{
    std::vector<SOFIPixelCombination> pixelCombinations;
    {
        SOFIPixelCombination pixelCombination;
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
        SOFIPixelCombination pixelCombination;
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
        SOFIPixelCombination pixelCombination;
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
        SOFIPixelCombination pixelCombination;
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
std::vector<SOFIPixelCombination> sofiPixelCombinations3()
{
    std::vector<SOFIPixelCombination> pixelCombinations;
    {
        SOFIPixelCombination pixelCombination;
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
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 0;
        pixelCombination.outputDeltaY = 1;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -2));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 0;
        pixelCombination.outputDeltaY = 2;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
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
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 1;
        pixelCombination.outputDeltaY = 1;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 1;
        pixelCombination.outputDeltaY = 2;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 2;
        pixelCombination.outputDeltaY = 0;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 2;
        pixelCombination.outputDeltaY = 1;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(2, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 2;
        pixelCombination.outputDeltaY = 2;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    return pixelCombinations;
}
std::vector<SOFIPixelCombination> sofiPixelCombinations4()
{
    std::vector<SOFIPixelCombination> pixelCombinations;
    {
        SOFIPixelCombination pixelCombination;
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
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 0;
        pixelCombination.outputDeltaY = 1;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -2));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 0;
        pixelCombination.outputDeltaY = 2;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 0;
        pixelCombination.outputDeltaY = 3;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
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
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
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
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 1;
        pixelCombination.outputDeltaY = 2;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
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
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -2));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 2;
        pixelCombination.outputDeltaY = 0;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
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
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(2, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 2;
        pixelCombination.outputDeltaY = 2;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 2;
        pixelCombination.outputDeltaY = 3;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
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
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
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
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(2, -1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
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
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, -1));
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
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(2, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 3;
        pixelCombination.outputDeltaY = 2;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 3;
        pixelCombination.outputDeltaY = 3;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 2));
            subset.push_back(std::pair<int, int>(2, -2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    return pixelCombinations;
}

std::vector<SOFIPixelCombination> sofiPixelCombinations5()
{
    std::vector<SOFIPixelCombination> pixelCombinations;
    {
        SOFIPixelCombination pixelCombination;
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
        SOFIPixelCombination pixelCombination;
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
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, -2));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 1));
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
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 0;
        pixelCombination.outputDeltaY = 2;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, -2));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 0;
        pixelCombination.outputDeltaY = 3;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 0;
        pixelCombination.outputDeltaY = 4;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
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
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-2, 0));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
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
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 1));
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
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 1;
        pixelCombination.outputDeltaY = 2;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -2));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -2));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, -1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 1;
        pixelCombination.outputDeltaY = 3;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -2));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 1;
        pixelCombination.outputDeltaY = 4;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 2;
        pixelCombination.outputDeltaY = 0;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(2, 0));
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
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(2, -1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
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
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, -2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
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
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
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
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
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
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
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
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
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
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 3;
        pixelCombination.outputDeltaY = 0;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, -1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
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
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(2, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -2));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 3;
        pixelCombination.outputDeltaY = 2;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
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
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
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
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 3;
        pixelCombination.outputDeltaY = 4;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
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
        SOFIPixelCombination pixelCombination;
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
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(2, -1));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
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
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 4;
        pixelCombination.outputDeltaY = 2;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, -2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 2));
            subset.push_back(std::pair<int, int>(2, -2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 4;
        pixelCombination.outputDeltaY = 3;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 2));
            subset.push_back(std::pair<int, int>(2, -2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, -1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 4;
        pixelCombination.outputDeltaY = 4;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            subset.push_back(std::pair<int, int>(2, -2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 2));
            subset.push_back(std::pair<int, int>(2, -1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    return pixelCombinations;
}

std::vector<SOFIPixelCombination> sofiPixelCombinations6()
{
    std::vector<SOFIPixelCombination> pixelCombinations;
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 0;
        pixelCombination.outputDeltaY = 0;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, -2));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, -2));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, -2));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 0;
        pixelCombination.outputDeltaY = 1;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, -2));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, -2));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, -2));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 0;
        pixelCombination.outputDeltaY = 2;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, -2));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, -2));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, -2));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 0;
        pixelCombination.outputDeltaY = 3;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, -2));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, -2));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, -2));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 0;
        pixelCombination.outputDeltaY = 4;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, -2));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, -2));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
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
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 0;
        pixelCombination.outputDeltaY = 5;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, -2));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
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
            subset.push_back(std::pair<int, int>(0, 2));
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
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, -2));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 1;
        pixelCombination.outputDeltaY = 1;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, -2));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
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
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 1;
        pixelCombination.outputDeltaY = 2;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -2));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, -1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 1;
        pixelCombination.outputDeltaY = 3;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 1;
        pixelCombination.outputDeltaY = 4;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 1;
        pixelCombination.outputDeltaY = 5;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 2;
        pixelCombination.outputDeltaY = 0;
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
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, -1));
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-2, 1));
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(2, -1));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 2;
        pixelCombination.outputDeltaY = 1;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
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
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, -1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
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
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 2;
        pixelCombination.outputDeltaY = 2;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
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
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -2));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -2));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
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
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -2));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -2));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 2;
        pixelCombination.outputDeltaY = 4;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -2));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 2;
        pixelCombination.outputDeltaY = 5;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 2));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 3;
        pixelCombination.outputDeltaY = 0;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(2, -1));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, -1));
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(2, 0));
            subset.push_back(std::pair<int, int>(2, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 3;
        pixelCombination.outputDeltaY = 1;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, -2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, -1));
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
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 3;
        pixelCombination.outputDeltaY = 2;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
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
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
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
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 3;
        pixelCombination.outputDeltaY = 4;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
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
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
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
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 2));
            subset.push_back(std::pair<int, int>(2, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 4;
        pixelCombination.outputDeltaY = 0;
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
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(2, -1));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(2, 0));
            subset.push_back(std::pair<int, int>(2, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-2, 0));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, -1));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
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
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -2));
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
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, -1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 4;
        pixelCombination.outputDeltaY = 2;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 0));
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
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, -2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 4;
        pixelCombination.outputDeltaY = 3;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
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
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 2));
            subset.push_back(std::pair<int, int>(2, -2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 4;
        pixelCombination.outputDeltaY = 4;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            subset.push_back(std::pair<int, int>(2, -2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 2));
            subset.push_back(std::pair<int, int>(2, -1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 4;
        pixelCombination.outputDeltaY = 5;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            subset.push_back(std::pair<int, int>(2, -1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 2));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
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
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 5;
        pixelCombination.outputDeltaY = 0;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, -1));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-1, 0));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(2, 0));
            subset.push_back(std::pair<int, int>(2, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(-2, 0));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(2, -1));
            subset.push_back(std::pair<int, int>(2, 0));
            subset.push_back(std::pair<int, int>(2, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
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
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, -1));
            subset.push_back(std::pair<int, int>(0, 0));
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, -1));
            subset.push_back(std::pair<int, int>(2, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 5;
        pixelCombination.outputDeltaY = 2;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 2));
            subset.push_back(std::pair<int, int>(2, -2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, -1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            subset.push_back(std::pair<int, int>(2, -2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 5;
        pixelCombination.outputDeltaY = 3;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            subset.push_back(std::pair<int, int>(2, -2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 2));
            subset.push_back(std::pair<int, int>(2, -1));
            pixelCombination.combinations.push_back(subset);
        }
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
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -2));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            subset.push_back(std::pair<int, int>(2, -1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 5;
        pixelCombination.outputDeltaY = 4;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            subset.push_back(std::pair<int, int>(2, -2));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            subset.push_back(std::pair<int, int>(2, -1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 2));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 1));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    {
        SOFIPixelCombination pixelCombination;
        pixelCombination.outputDeltaX = 5;
        pixelCombination.outputDeltaY = 5;
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            subset.push_back(std::pair<int, int>(2, -1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(1, 2));
            subset.push_back(std::pair<int, int>(2, 0));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 2));
            subset.push_back(std::pair<int, int>(2, 1));
            pixelCombination.combinations.push_back(subset);
        }
        {
            std::vector<std::pair<int, int> > subset;
            subset.push_back(std::pair<int, int>(0, 1));
            subset.push_back(std::pair<int, int>(0, 2));
            subset.push_back(std::pair<int, int>(1, -1));
            subset.push_back(std::pair<int, int>(1, 0));
            subset.push_back(std::pair<int, int>(1, 1));
            subset.push_back(std::pair<int, int>(2, 2));
            pixelCombination.combinations.push_back(subset);
        }
        pixelCombinations.push_back(pixelCombination);
    }
    return pixelCombinations;
}
