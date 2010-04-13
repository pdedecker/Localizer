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

#include "XOPStandardHeaders.h"
#include <boost/numeric/ublas/matrix.hpp>
#include "PALM_analysis_defines.h"
#include "PALM_analysis_storage.h"
#include "PALM_analysis_FileIO.h"

namespace ublas = boost::numeric::ublas;

class ImageLoader;

// Routines that return information on CCD files and image frames to Igor
int load_partial_ccd_image(ImageLoader *image_loader, size_t n_start, size_t n_end);

int parse_ccd_headers(ImageLoader *image_loader);

// Routines that process CCD files and return data to Igor
int construct_summed_intensity_trace(ImageLoader *image_loader, DataFolderAndName outputWaveParams, long startX, long startY, long endX, long endY);

int construct_average_image(ImageLoader *image_loader, DataFolderAndName outputWaveParams, long startX, long startY, long endX, long endY);

void calculateStandardDeviationImage(ImageLoader *image_loader, DataFolderAndName outputWaveParams, long startX, long startY, long endX, long endY);

// Routines that can fetch and make waves from datafolders
waveHndl FetchWaveUsingFullPath(std::string wavePath);
waveHndl MakeWaveUsingFullPath(std::string wavePath, long *dimensionSizes, int type, int overwrite);

// Routines to convert between handles and C strings
int ConvertHandleToString(Handle handle, std::string& convertedString);

int ConvertHandleToFilepathString(Handle handle, std::string &output_string);

// routines to convert data from and to Igor format
waveHndl CopyVectorToIgorDPWave(boost::shared_ptr<std::vector<double> > vec, std::string waveName);

boost::shared_ptr<ublas::matrix<double> > CopyIgorDPWaveToMatrix(waveHndl wave);

waveHndl CopyMatrixToIgorDPWave(boost::shared_ptr<ublas::matrix<double> > matrix, std::string waveName);

waveHndl copy_PALMMatrix_float_to_IgorFPWave(boost::shared_ptr<ublas::matrix<float> > matrix, std::string waveName);

#endif
