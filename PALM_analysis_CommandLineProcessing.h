/*
 *  PALM_analysis_CommandLineProcessing.h
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 01/06/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include <string>
#include "boost/program_options.hpp"
#include "PALM_analysis.h"
#include "PALM_analysis_defines.h"
#include "PALM_analysis_FileIO.h"
#include "PALM_analysis_Processing.h"


// the main function
int main(int argc, char *argv[]);

// function that will parse the string command line arguments
// and convert them to the appropriate objects
boost::shared_ptr<CCDImagesProcessor> GetCCDImagesProcessor(std::string name, boost::shared_ptr<PALMAnalysisProgressReporter> progressReporter, size_t nFramesAveraging, double cameraMultiplier, double cameraOffset);

// function that will guess the CCD file type and return an image loader
boost::shared_ptr<ImageLoader> GetImageLoader(std::string filePath);

// function that will provide an appropriate output writer for CCD processing
boost::shared_ptr<ImageOutputWriter> GetImageOutputWriter(std::string processMethodName, int originalStorageFormat, std::string requestedFormat, std::string outputFilePath, size_t compression);

// function that takes the CCD file path and return an output string for the processed images
std::string GetOutputProcessedImagesFilePath(std::string dataFilePath, std::string outputFormat);
