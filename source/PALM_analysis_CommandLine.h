/*
 Copyright 2008-2011 Peter Dedecker.
 
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
 with analysis programs such as Igor Pro or Matlab, 
 the licensors of this Program grant you additional permission 
 to convey the resulting work.
 */

#include <string>
#include <iostream>
#include <algorithm>
#include <sstream>
#include "boost/program_options.hpp"
#include "PALM_analysis_defines.h"
#include "PALM_analysis_segmentation.h"
#include "PALM_analysis_ParticleFinding.h"
#include "PALM_analysis_Localization.h"
#include "PALM_analysis_storage.h"
#include "PALM_analysis.h"
#include "PALM_analysis_FileIO.h"
#include "PALM_analysis_ProgressReporting.h"


// the main function
int main(int argc, char *argv[]);

// function that will parse the string command line arguments
// and convert them to integer constants
boost::shared_ptr<ThresholdImage_Preprocessor> GetPreProcessorType(std::string name);
boost::shared_ptr<ThresholdImage_Postprocessor> GetPostProcessorType(std::string name);
boost::shared_ptr<ThresholdImage> GetSegmentationType(std::string name, double pfa, double threshold, double psfWidth, double smoothSigmaFactor);
boost::shared_ptr<ParticleFinder> GetParticleFinderType(std::string name);
boost::shared_ptr<ParticleVerifier> GetParticleVerifierType(std::string name, double psfWidth, double sigma);
boost::shared_ptr<FitPositions> GetPositionsFitter(std::string name, double psfWidth);

// function that will guess the CCD file type and return an image loader
boost::shared_ptr<ImageLoader> GetImageLoader(std::string filePath);

// function that takes the CCD file path and return an output string for the localized positions
std::string GetOutputPositionsFilePath(std::string dataFilePath);
