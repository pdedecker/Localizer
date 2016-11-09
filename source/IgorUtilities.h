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
#include "boost/smart_ptr.hpp"

#include "XOPStandardHeaders.h"
#include <eigen3/Eigen/Eigen>
#include "PALMAnalysis.h"
#include "Defines.h"
#include "Storage.h"

class ImageLoader;
class ImageOutputWriter;

// Routines that return information on CCD files and image frames to Igor
int LoadPartialCCDImage(ImageLoader *image_loader, int firstImage, int nImagesRequested, int overwrite, DataFolderAndName destination, 
                        std::shared_ptr<ProgressReporter> progressReporter);
int ParseCCDHeaders(ImageLoader *image_loader);
void WriteImagesToDisk(std::shared_ptr<ImageLoader> imageLoader, std::shared_ptr<ImageOutputWriter> outputWriter, std::shared_ptr<ProgressReporter> progressReporter, int firstImage, int nImagesToWrite);

// Routines that process CCD files and return data to Igor
std::vector<double> ConstructIntensityTrace(std::shared_ptr<ImageLoader> imageLoader, std::shared_ptr<ProgressReporter> progressReporter, bool doAverage);

std::vector<double> ConstructSummedIntensityTrace(std::shared_ptr<ImageLoader> imageLoader, std::shared_ptr<ProgressReporter> progressReporter);

std::vector<double> ConstructAverageIntensityTrace(std::shared_ptr<ImageLoader> imageLoader, std::shared_ptr<ProgressReporter> progressReporter);

ImagePtr ConstructAverageImage(std::shared_ptr<ImageLoader> imageLoader, std::shared_ptr<ProgressReporter> progressReporter);

ImagePtr ConstructVarianceImage(std::shared_ptr<ImageLoader> imageLoader, std::shared_ptr<ProgressReporter> progressReporter);

// Routines that can fetch and make waves from datafolders
waveHndl FetchWaveUsingFullPath(std::string wavePath);
waveHndl MakeWaveUsingFullPath(std::string wavePath, CountInt *dimensionSizes, int type, int overwrite);

// Routines to convert between handles and C strings
std::string ConvertHandleToString(Handle handle);

std::string ConvertPathToNativePath(std::string filePath);

// routines to convert data from and to Igor format
waveHndl CopyVectorToIgorDPWave(const std::vector<double>& vec, std::string waveName);
waveHndl CopyVectorToIgorDPWave(const std::vector<double>& vec, DataFolderAndName outputWaveParams);
std::vector<double> IgorWaveToVector(waveHndl wav);

ImagePtr CopyIgorDPWaveToMatrix(waveHndl wave);

waveHndl CopyMatrixToIgorDPWave(ImagePtr matrix, std::string waveName);
waveHndl CopyMatrixToIgorDPWave(const Eigen::MatrixXd& matrix, std::string waveName);
waveHndl CopyMatrixToIgorDPWave(ImagePtr matrix, DataFolderAndName dataFolderAndName);
waveHndl CopyStackToIgorDPWave(std::vector<ImagePtr> stack, DataFolderAndName dataFolderAndName);

void SetWaveNote(waveHndl wav, const std::string& note);

void PrintToHistory(const std::string& str);

#endif
