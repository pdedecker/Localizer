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

#include "XOPStandardHeaders.h"

#include "PALM_analysis_defines.h"
#include "PALM_analysis_storage.h"
#include "PALM_analysis_FileIO.h"

// Routines that return information on CCD files and image frames to Igor
int load_partial_ccd_image(ImageLoader *image_loader, size_t n_start, size_t n_end);

int parse_ccd_headers(ImageLoader *image_loader);

// Routines that process CCD files and return data to Igor
int construct_summed_intensity_trace(ImageLoader *image_loader, std::string output_wave_name, long startX, long startY, long endX, long endY);

int construct_average_image(ImageLoader *image_loader, std::string output_wave_name, long startX, long startY, long endX, long endY);

void calculateStandardDeviationImage(ImageLoader *image_loader, std::string output_wave_name, long startX, long startY, long endX, long endY);


// Routines to convert between handles and C strings
int ConvertHandleToString(Handle handle, std::string& convertedString);

int ConvertHandleToFilepathString(Handle handle, std::string &output_string);

// routines to convert data from and to Igor format
waveHndl copy_vector_to_IgorDPWave(boost::shared_ptr<std::vector<double> > vec, std::string waveName);

boost::shared_ptr<PALMMatrix<double> > copy_IgorDPWave_to_gsl_matrix(waveHndl wave);

waveHndl copy_PALMMatrix_to_IgorDPWave(boost::shared_ptr<PALMMatrix<double> > matrix, std::string waveName);

waveHndl copy_PALMMatrix_float_to_IgorFPWave(boost::shared_ptr<PALMMatrix<float> > matrix, std::string waveName);

boost::shared_ptr<PALMVolume <double> > copy_IgorDPWave_to_gsl_volume(waveHndl wave);

waveHndl copy_PALMVolume_to_IgorDPWave(boost::shared_ptr<PALMVolume <double> > volume, std::string waveName);

waveHndl copy_PALMVolume_ushort_to_IgorUINT16wave(boost::shared_ptr<PALMVolume<unsigned short> > volume, std::string waveName);

#endif
