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
#include <tuple>

#include "NewSOFI.h"

class PixelCoordinate;
typedef std::vector<PixelCoordinate> PixelCombination;
typedef std::vector<PixelCombination> Partition;
typedef std::vector<Partition> GroupOfPartitions;

class PixelCoordinate {
public:
    PixelCoordinate() : dx(-10), dy(-10), dt(-10) {}
    PixelCoordinate(int xx, int yy, int tt) : dx(xx), dy(yy), dt(tt) {}
    
    int dx;
    int dy;
    int dt;
};

class SOFIVirtualPixel {
public:
    int outputDeltaX;
    int outputDeltaY;
    std::vector<PixelCombination> combinations;
    std::vector<double> scores;
};

class SOFIKernel {
public:
    int outputDeltaX;
    int outputDeltaY;
    std::vector<GroupOfPartitions> combinations;
    std::vector<PixelCombination> pixelCombinations;
    std::vector<double> scores;
    std::vector<double> requestedCombinationWeights;
};

std::vector<SOFIKernel> KernelsForOrder(const int order, const std::vector<int>& timeLags, const SOFIOptions::AllowablePixelCombinations allowablePixelCombinations, const double pixelCombinationCutoff);
std::vector<SOFIKernel> AutoKernelsForOrder(const int order, const std::vector<int>& timeLags);
std::vector<SOFIKernel> KernelsFromExternalDescription(const Eigen::MatrixXd& kernelDescription, const std::vector<int>& timeLags);

GroupOfPartitions AllPartitions(PixelCombination pixelCombination);
std::string PrintPartition(const Partition& partition);
std::string ParititionTest();

inline bool operator<(const PixelCoordinate& lhs, const PixelCoordinate& rhs) {
    if (lhs.dx == rhs.dx) {
        if (lhs.dy == rhs.dy) {
            return (lhs.dt < rhs.dt);
        } else {
            return (lhs.dy < rhs.dy);
        }
    } else {
        return (lhs.dx < rhs.dx);
    }
}

inline bool operator==(const PixelCoordinate& lhs, const PixelCoordinate& rhs) {
    return ((lhs.dx == rhs.dx) && (lhs.dy == rhs.dy) && (lhs.dt == rhs.dt));
}

class ComparePixelCombinations {
public:
    bool operator() (const PixelCombination& comb1, const PixelCombination& comb2) const;
};

#endif
