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


// Routines to convert between handles and C strings
int ConvertHandleToString(Handle handle, std::string& convertedString);

int ConvertHandleToFilepathString(Handle handle, string &output_string);

// routines to convert data from and to Igor format
boost::shared_ptr<PALMMatrix<double> > copy_IgorDPWave_to_gsl_matrix(waveHndl wave);

waveHndl copy_PALMMatrix_to_IgorDPWave(boost::shared_ptr<PALMMatrix<double> > matrix, string waveName);

waveHndl copy_PALMMatrix_float_to_IgorFPWave(boost::shared_ptr<PALMMatrix<float> > matrix, string waveName);

boost::shared_ptr<PALMVolume <double> > copy_IgorDPWave_to_gsl_volume(waveHndl wave);

waveHndl copy_PALMVolume_to_IgorDPWave(boost::shared_ptr<PALMVolume <double> > volume, string waveName);

waveHndl copy_PALMVolume_ushort_to_IgorUINT16wave(boost::shared_ptr<PALMVolume<unsigned short> > volume, string waveName);

#endif
