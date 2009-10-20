/*
 *  PALM_analysis_IgorXOP.h
 *  PALM analysis
 *
 *  Created by Peter Dedecker on 05/01/2008.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef PALM_ANALYSIS_IGOR_XOP
#define PALM_ANALYSIS_IGOR_XOP

#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include "tiffio.h"
#include "boost/smart_ptr.hpp"
#include "XOPStandardHeaders.h"
#include "PALM_analysis.h"
#include "PALM_analysis_defines.h"
#include "PALM_analysis_errors.h"
#include "PALM_analysis_Processing.h"
#include "PALM_analysis_FileIO.h"
#include "PALM_analysis_segmentation.h"
#include "PALM_analysis_Localization.h"
#include "PALM_analysis_ParticleFinding.h"
#include "PALM_analysis_PALMImages.h"



using namespace std;

HOST_IMPORT int main(IORecHandle ioRecHandle);

int ConvertHandleToFilepathString(Handle handle, size_t cameraType, string &output_string);

boost::shared_ptr<PALMMatrix<double> > copy_IgorDPWave_to_gsl_matrix(waveHndl wave);

waveHndl copy_PALMMatrix_to_IgorDPWave(boost::shared_ptr<PALMMatrix<double> > matrix, string waveName);
waveHndl copy_PALMMatrix_float_to_IgorFPWave(boost::shared_ptr<PALMMatrix<float> > matrix, string waveName);

boost::shared_ptr<ImageLoader> get_image_loader_for_camera_type(size_t camera_type, string data_file_path, size_t cache_size = N_SIMULTANEOUS_IMAGE_LOADS);

boost::shared_ptr<PALMVolume <double> > copy_IgorDPWave_to_gsl_volume(waveHndl wave);

waveHndl copy_PALMVolume_to_IgorDPWave(boost::shared_ptr<PALMVolume <double> > volume, string waveName);

waveHndl copy_PALMVolume_ushort_to_IgorUINT16wave(boost::shared_ptr<PALMVolume<unsigned short> > volume, string waveName);

#endif

