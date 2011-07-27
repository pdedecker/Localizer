/*
 *  PALM_analysis_IgorUtilities.h
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 30/12/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef PALM_ANALYSIS_IGORUTILITIES
#define PALM_ANALYSIS_IGORUTILITIES

#include <string>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include "boost/smart_ptr.hpp"

#include "XOPStandardHeaders.h"
#include <eigen3/Eigen/Eigen>
#include "PALM_analysis.h"
#include "PALM_analysis_ProgressReporting.h"
#include "PALM_analysis_defines.h"
#include "PALM_analysis_storage.h"
#include "PALM_analysis_FileIO.h"

class ImageLoader;

/**
 * Igor passes strings as a handle, which needs special conversion
 * also try to determine the camera type.
 * The difficulty here is in trying to distinguish between a file path
 * and a path to an Igor wave. So do the following: first try to treat the
 * the string as a path to a wave. If that wave doesn't exist then assume it's a
 * file path.
 */
void GetFilePathAndCameraType(std::string& inputFilePath, std::string &filePath, size_t &cameraType);

int GetFileStorageType(std::string &filePath);

boost::shared_ptr<ImageLoader> GetImageLoader(size_t camera_type, std::string& data_file_path);

// Routines that return information on CCD files and image frames to Igor
int LoadPartialCCDImage(ImageLoader *image_loader, size_t firstImage, size_t nImagesRequested, int overwrite, DataFolderAndName destination, 
						   boost::shared_ptr<ProgressReporter> progressReporter);

int ParseCCDHeaders(ImageLoader *image_loader);

// Routines that process CCD files and return data to Igor
waveHndl construct_summed_intensity_trace(ImageLoader *image_loader, DataFolderAndName outputWaveParams, 
										  long startX, long startY, long endX, long endY, 
										  boost::shared_ptr<ProgressReporter> progressReporter);

waveHndl construct_average_intensity_trace(ImageLoader *image_loader, DataFolderAndName outputWaveParams, 
										   long startX, long startY, long endX, long endY, 
										   boost::shared_ptr<ProgressReporter> progressReporter);

waveHndl construct_average_image(ImageLoader *image_loader, DataFolderAndName outputWaveParams, 
								 long startX, long startY, long endX, long endY,
								 boost::shared_ptr<ProgressReporter> progressReporter);

waveHndl calculateVarianceImage(ImageLoader *image_loader, DataFolderAndName outputWaveParams, 
								long startX, long startY, long endX, long endY,
								boost::shared_ptr<ProgressReporter> progressReporter);

// Routines that can fetch and make waves from datafolders
waveHndl FetchWaveUsingFullPath(std::string wavePath);
waveHndl MakeWaveUsingFullPath(std::string wavePath, long *dimensionSizes, int type, int overwrite);

// Routines to convert between handles and C strings
std::string ConvertHandleToString(Handle handle);

std::string ConvertPathToNativePath(std::string filePath);

// routines to convert data from and to Igor format
waveHndl CopyVectorToIgorDPWave(boost::shared_ptr<std::vector<double> > vec, std::string waveName);

ImagePtr CopyIgorDPWaveToMatrix(waveHndl wave);

waveHndl CopyMatrixToIgorDPWave(ImagePtr matrix, std::string waveName);

#endif
