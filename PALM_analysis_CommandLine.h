/*
 *  PALM_analysis_CommandLine.h
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 30/12/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
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
#include "PALM_analysis_Processing.h"


// the main function
int main(int argc, char *argv[]);

// function that will parse the string command line arguments
// and convert them to the appropriate objects
boost::shared_ptr<ThresholdImage_Preprocessor> GetPreProcessorType(std::string name);
boost::shared_ptr<ThresholdImage_Postprocessor> GetPostProcessorType(std::string name);
boost::shared_ptr<ThresholdImage> GetSegmentationType(std::string name, double pfa, double threshold, double psfWidth);
boost::shared_ptr<ParticleFinder> GetParticleFinderType(std::string name);
boost::shared_ptr<ParticleVerifier> GetParticleVerifierType(std::string name, double psfWidth, double sigma);
boost::shared_ptr<FitPositions> GetPositionsFitter(std::string name, double psfWidth);
boost::shared_ptr<CCDImagesProcessor> GetCCDImagesProcessor(std::string name, size_t nFramesAveraging, double cameraMultiplier, double cameraOffset);

// function that will guess the CCD file type and return an image loader
boost::shared_ptr<ImageLoader> GetImageLoader(std::string filePath);

// function that will provide an appropriate output writer for CCD processing
boost::shared_ptr<ImageOutputWriter> GetImageOutputWriter(std::string processMethodName, std::string outputFilePath, size_t compression);

// function that takes the CCD file path and return an output string for the localized positions
std::string GetOutputPositionsFilePath(std::string dataFilePath);
std::string GetOutputProcessedImagesFilePath(std::string dataFilePath);
