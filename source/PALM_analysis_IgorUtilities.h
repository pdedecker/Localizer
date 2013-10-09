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

std::shared_ptr<ImageLoader> GetImageLoader(size_t camera_type, std::string& data_file_path);

// Routines that return information on CCD files and image frames to Igor
int LoadPartialCCDImage(ImageLoader *image_loader, size_t firstImage, size_t nImagesRequested, int overwrite, DataFolderAndName destination, 
                        std::shared_ptr<ProgressReporter> progressReporter);

int ParseCCDHeaders(ImageLoader *image_loader);

// Routines that process CCD files and return data to Igor
std::vector<double> ConstructIntensityTrace(std::shared_ptr<ImageLoader> imageLoader, std::shared_ptr<ProgressReporter> progressReporter, int startX, int startY, int endX, int endY, bool doAverage);

std::vector<double> ConstructSummedIntensityTrace(std::shared_ptr<ImageLoader> imageLoader, std::shared_ptr<ProgressReporter> progressReporter, int startX, int startY, int endX, int endY);

std::vector<double> ConstructAverageIntensityTrace(std::shared_ptr<ImageLoader> imageLoader, std::shared_ptr<ProgressReporter> progressReporter, int startX, int startY, int endX, int endY);

waveHndl construct_summed_intensity_trace(ImageLoader *image_loader, DataFolderAndName outputWaveParams, 
                                          long startX, long startY, long endX, long endY,
                                          std::shared_ptr<ProgressReporter> progressReporter);

waveHndl construct_average_intensity_trace(ImageLoader *image_loader, DataFolderAndName outputWaveParams, 
                                           long startX, long startY, long endX, long endY,
                                           std::shared_ptr<ProgressReporter> progressReporter);

waveHndl construct_average_image(ImageLoader *image_loader, DataFolderAndName outputWaveParams, 
                                 long startX, long startY, long endX, long endY,
                                 std::shared_ptr<ProgressReporter> progressReporter);

waveHndl calculateVarianceImage(ImageLoader *image_loader, DataFolderAndName outputWaveParams, 
                                long startX, long startY, long endX, long endY,
                                std::shared_ptr<ProgressReporter> progressReporter);

// Routines that can fetch and make waves from datafolders
waveHndl FetchWaveUsingFullPath(std::string wavePath);
waveHndl MakeWaveUsingFullPath(std::string wavePath, CountInt *dimensionSizes, int type, int overwrite);

// Routines to convert between handles and C strings
std::string ConvertHandleToString(Handle handle);

std::string ConvertPathToNativePath(std::string filePath);

// routines to convert data from and to Igor format
waveHndl CopyVectorToIgorDPWave(std::shared_ptr<std::vector<double> > vec, std::string waveName);
waveHndl CopyVectorToIgorDPWave(const std::vector<double>& vec, DataFolderAndName outputWaveParams);

ImagePtr CopyIgorDPWaveToMatrix(waveHndl wave);

waveHndl CopyMatrixToIgorDPWave(ImagePtr matrix, std::string waveName);
waveHndl CopyMatrixToIgorDPWave(ImagePtr matrix, DataFolderAndName dataFolderAndName);

#endif
