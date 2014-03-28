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

#ifndef NEWSOFIKERNELS_H
#define NEWSOFIKERNELS_H

#include <vector>
#include <utility>

//typedef std::vector<std::pair<int, int> > PixelCombination;
//typedef std::vector<PixelCombination> Subset;
//typedef std::vector<Subset> Partition;
//typedef std::vector<Partition> GroupOfPartitions;

typedef std::vector<std::pair<int, int> > PixelCombination;
typedef std::vector<PixelCombination> Partition;
typedef std::vector<Partition> GroupOfPartitions;

class SOFIPixelCombination {
public:
    int outputDeltaX;
    int outputDeltaY;
    std::vector<PixelCombination> combinations;
};

class SOFIKernel {
public:
    int outputDeltaX;
    int outputDeltaY;
    std::vector<GroupOfPartitions> combinations;
};

std::vector<SOFIKernel> KernelsForOrder(const int order);

class ComparePixelCombinations {
public:
    bool operator() (const PixelCombination& comb1, const PixelCombination& comb2) const;
};

GroupOfPartitions AllPartitions(PixelCombination pixelCombination);
std::string PrintPartition(const Partition& partition);
std::string ParititionTest();

std::vector<SOFIPixelCombination> sofiPixelCombinations2();
std::vector<SOFIPixelCombination> sofiPixelCombinations3();
std::vector<SOFIPixelCombination> sofiPixelCombinations4();
std::vector<SOFIPixelCombination> sofiPixelCombinations5();
std::vector<SOFIPixelCombination> sofiPixelCombinations6();

#endif
