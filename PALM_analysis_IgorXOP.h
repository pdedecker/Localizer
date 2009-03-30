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
#include "boost/smart_ptr.hpp"
#include "XOPStandardHeaders.h"
#include "PALM_analysis.h"
#include "PALM_analysis_defines.h"
#include "PALM_analysis_errors.h"
#include "PALM_analysis_FileIO.h"



using namespace std;

HOST_IMPORT int main(IORecHandle ioRecHandle);

int convert_handle_to_string(Handle handle, string &output_string);

boost::shared_ptr<encap_gsl_matrix> copy_IgorDPWave_to_gsl_matrix(waveHndl wave);

waveHndl copy_gsl_matrix_to_IgorDPWave(boost::shared_ptr<encap_gsl_matrix> matrix, string waveName);

boost::shared_ptr<ImageLoader> get_image_loader_for_camera_type(unsigned long camera_type, string data_file_path, unsigned long cache_size = N_SIMULTANEOUS_IMAGE_LOADS);

boost::shared_ptr<encap_gsl_volume> copy_IgorDPWave_to_gsl_volume(waveHndl wave);

waveHndl copy_gsl_volume_to_IgorDPWave(boost::shared_ptr<encap_gsl_volume> volume, string waveName);

waveHndl copy_gsl_volume_ushort_to_IgorUINT16wave(boost::shared_ptr<encap_gsl_volume_ushort> volume, string waveName);

#endif

