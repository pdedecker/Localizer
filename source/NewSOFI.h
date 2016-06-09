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

#ifndef NEWSOFI_H
#define NEWSOFI_H

#include <vector>
#include <memory>
#include <functional>

#include "Defines.h"
#include "Storage.h"
#include "FileIO.h"
#include "ProgressReporting.h"
#include "SOFIFrameVerifiers.h"

class SOFIOptions {
public:
    enum AllowablePixelCombinations {
        NonOverlappingPixels,
        AllowOverlappingPixels
    };
    
    SOFIOptions() :
        nFramesToSkip(0),
        nFramesToInclude(-1),
        batchSize(0),
        wantCrossCumulant(true),
        allowablePixelCombinations(NonOverlappingPixels),
        doPixelationCorrection(true),
        alsoCorrectVariance(true),
        wantAverageImage(false),
        wantJackKnife(false),
        pixelCombinationCutoff(1.0),
        haveExternalPixelCombination(false),
        wantDebugMessages(false)
    {
        orders.push_back(2);
    }
    
    std::vector<int> orders;
    int nFramesToSkip;
    int nFramesToInclude;
    int batchSize;
    bool wantCrossCumulant;
    AllowablePixelCombinations allowablePixelCombinations;
    std::vector<int> lagTimes;
    bool doPixelationCorrection;
    bool alsoCorrectVariance;
    bool wantAverageImage;
    ImagePtr averageImage;
    bool wantJackKnife;
    std::function<ExternalImageBuffer(int,int,int,double,double,double,int,int,bool)> jackKnifeAllocator;
    std::vector<ExternalImageBuffer> jackKnifeImages;
    double pixelCombinationCutoff;
    bool haveExternalPixelCombination;
    Eigen::MatrixXd externalPixelCombinations;
    std::vector<double> pixelCombinationWeights;
    bool wantDebugMessages;
    
    int outputUsedBoundaryMargin;
};

void DoNewSOFI(std::shared_ptr<ImageLoader> imageLoader, SOFIOptions& options, std::shared_ptr<ProgressReporter> progressReporter, std::vector<ImagePtr>& sofiOutputImages);

Eigen::MatrixXd PixelCombinationsForOrderAsMatrix(const int order, const SOFIOptions::AllowablePixelCombinations allowablePixelCombinations, const double pixelCombinationCutoff);

#endif
