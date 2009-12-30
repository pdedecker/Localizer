/*
 *  PALM_analysis_IgorUtilities.cpp
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 30/12/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "PALM_analysis_IgorUtilities.h"

int ConvertHandleToString(Handle handle, std::string& convertedString) {
	int err;
	
	// determine the type of positions being passed
	size_t stringLength = GetHandleSize(handle);
	boost::scoped_array<char> CStringWaveNote(new char[stringLength + 1]);
	
	err = GetCStringFromHandle(handle, CStringWaveNote.get(), stringLength);
	if (err != 0)
		return err;
	
	// save the wavenote as a std::string
	convertedString.assign(CStringWaveNote.get());
	
	return 0;
}

int ConvertHandleToFilepathString(Handle handle, string &output_path) {
	int err;
	char handle_char[1024];
	char handle_char_POSIX[1024];
	
	err = GetCStringFromHandle(handle, handle_char, 1023);
	if (err != 0) {
		return err;
	}
	
#ifdef _MACINTOSH_
	err = WinToMacPath(handle_char);
	if (err != 0) {
		return err;
	}
	
	
	err = HFSToPosixPath(handle_char, handle_char_POSIX, 0);
	if (err != 0) {
		return err;
	}
	output_path.assign(handle_char_POSIX);
#endif
#ifdef _WINDOWS_
	err = MacToWinPath(handle_char);
	if (err != 0) {
		return err;
	}
	output_path.assign(handle_char);
#endif
	
	return 0;
	
}



boost::shared_ptr<PALMMatrix<double> > copy_IgorDPWave_to_gsl_matrix(waveHndl wave) {
	// copy a Igor wave into a new gsl_matrix
	
	int err;
	long numDimensions; 
	long dimensionSizes[MAX_DIMENSIONS+1];
	size_t x_size, y_size;
	long indices[MAX_DIMENSIONS];
	double value[2];
	
	
	err = MDGetWaveDimensions(wave, &numDimensions, dimensionSizes);
	if (err != 0) {
		throw err;
	}
	if (numDimensions != 2) {
		throw INCOMPATIBLE_DIMENSIONING;
	}
	
	x_size = dimensionSizes[0];
	y_size = dimensionSizes[1];
	
	boost::shared_ptr<PALMMatrix<double> > matrix(new PALMMatrix<double>(x_size, y_size));
	
	for (size_t j = 0; j < y_size; ++j) {
		for (size_t i = 0; i < x_size; ++i) {
			indices[0] = i;
			indices[1] = j;
			
			err = MDGetNumericWavePointValue(wave, indices, value);
			if (err != 0) {
				throw err;
			}
			
			matrix->set(i, j, value[0]);
		}
	}
	
	return matrix;
}

waveHndl copy_PALMMatrix_to_IgorDPWave(boost::shared_ptr<PALMMatrix<double> > matrix, string waveName) {
	
	waveHndl DPWave;
	
	int err;
	long indices[MAX_DIMENSIONS];
	long dimensionSizes[MAX_DIMENSIONS+1];
	double value[2];
	
	// special case:
	// if the matrix is NULL (such as when there are no positions found)
	// then we return an empty wave
	if (matrix.get() == NULL) {
		dimensionSizes[0] = 0;
		dimensionSizes[1] = 0;
		dimensionSizes[2] = 0;
		
		err = MDMakeWave(&DPWave, waveName.c_str(), NULL, dimensionSizes, NT_FP64, 1);
		if (err != 0) {
			throw err;
		}
		
		return DPWave;
		
	}
	
	
	size_t x_size = (size_t)matrix->getXSize();
	size_t y_size = (size_t)matrix->getYSize();
	
	dimensionSizes[0] = x_size;
	dimensionSizes[1] = y_size;
	dimensionSizes[2] = 0;
	
	err = MDMakeWave(&DPWave, waveName.c_str(), NULL, dimensionSizes, NT_FP64, 1);
	if (err != 0) {
		throw err;
	}
	
	for (size_t j = 0; j < y_size; ++j) {
		for (size_t i = 0; i < x_size; ++i) {
			indices[0] = i;
			indices[1] = j;
			
			value[0] = (*matrix)(i, j);
			
			err = MDSetNumericWavePointValue(DPWave, indices, value);
			if (err != 0) {
				throw err;
			}
		}
	}
	
	return DPWave;
}

waveHndl copy_PALMMatrix_float_to_IgorFPWave(boost::shared_ptr<PALMMatrix<float> > matrix, string waveName) {
	
	waveHndl DPWave;
	
	int err;
	long indices[MAX_DIMENSIONS];
	long dimensionSizes[MAX_DIMENSIONS+1];
	double value[2];
	
	// special case:
	// if the matrix is NULL (such as when there are no positions found)
	// then we return an empty wave
	if (matrix.get() == NULL) {
		dimensionSizes[0] = 0;
		dimensionSizes[1] = 0;
		dimensionSizes[2] = 0;
		
		err = MDMakeWave(&DPWave, waveName.c_str(), NULL, dimensionSizes, NT_FP32, 1);
		if (err != 0) {
			throw err;
		}
		
		return DPWave;
		
	}
	
	
	size_t x_size = (size_t)matrix->getXSize();
	size_t y_size = (size_t)matrix->getYSize();
	
	dimensionSizes[0] = x_size;
	dimensionSizes[1] = y_size;
	dimensionSizes[2] = 0;
	
	err = MDMakeWave(&DPWave, waveName.c_str(), NULL, dimensionSizes, NT_FP32, 1);
	if (err != 0) {
		throw err;
	}
	
	for (size_t j = 0; j < y_size; ++j) {
		for (size_t i = 0; i < x_size; ++i) {
			indices[0] = i;
			indices[1] = j;
			
			value[0] = (*matrix)(i, j);
			
			err = MDSetNumericWavePointValue(DPWave, indices, value);
			if (err != 0) {
				throw err;
			}
		}
	}
	
	return DPWave;
}


boost::shared_ptr<PALMVolume <double> > copy_IgorDPWave_to_gsl_volume(waveHndl wave) {
	int err;
	long indices[MAX_DIMENSIONS];
	long dimensionSizes[MAX_DIMENSIONS+1];
	long numDimensions;
	double value[2];
	size_t x_size, y_size, z_size;
	
	err = MDGetWaveDimensions(wave, &numDimensions, dimensionSizes);
	if (err != 0) {
		throw err;
	}
	if (numDimensions != 3) {
		throw INCOMPATIBLE_DIMENSIONING;
	}
	
	x_size = dimensionSizes[0];
	y_size = dimensionSizes[1];
	z_size = dimensionSizes[2];
	
	boost::shared_ptr<PALMVolume <double> > volume(new PALMVolume <double>(x_size, y_size, z_size));
	
	for (size_t k = 0; k < z_size; ++k)
		for (size_t j = 0; j < y_size; ++j) {
			for (size_t i = 0; i < x_size; ++i) {
				indices[0] = i;
				indices[1] = j;
				indices[2] = k;
				
				err = MDGetNumericWavePointValue(wave, indices, value);
				if (err != 0) {
					throw err;
				}
				
				(*volume)(i, j, k) = value[0];
			}
		}
	
	return volume;
}

waveHndl copy_PALMVolume_to_IgorDPWave(boost::shared_ptr<PALMVolume<double> > volume, string waveName) {
	
	waveHndl DPWave;
	
	int err;
	long indices[MAX_DIMENSIONS];
	long dimensionSizes[MAX_DIMENSIONS+1];
	double value[2];
	
	
	size_t x_size = (size_t)volume->getXSize();
	size_t y_size = (size_t)volume->getYSize();
	size_t z_size = (size_t)volume->getZSize();
	
	dimensionSizes[0] = x_size;
	dimensionSizes[1] = y_size;
	dimensionSizes[2] = z_size;
	dimensionSizes[3] = 0;
	
	err = MDMakeWave(&DPWave, waveName.c_str(), NULL, dimensionSizes, NT_FP64, 1);
	if (err != 0) {
		throw err;
	}
	
	for (size_t k = 0; k < z_size; ++k) {
		for (size_t j = 0; j < y_size; ++j) {
			for (size_t i = 0; i < x_size; ++i) {
				indices[0] = i;
				indices[1] = j;
				indices[2] = k;
				
				value[0] = (*volume)(i, j, k);
				
				err = MDSetNumericWavePointValue(DPWave, indices, value);
				if (err != 0) {
					throw err;
				}
			}
		}
	}
	
	return DPWave;
}


waveHndl copy_PALMVolume_ushort_to_IgorUINT16wave(boost::shared_ptr<PALMVolume<unsigned short> > volume, string waveName) {
	
	waveHndl DPWave;
	
	int err;
	long indices[MAX_DIMENSIONS];
	long dimensionSizes[MAX_DIMENSIONS+1];
	double value[2];
	
	
	size_t x_size = (size_t)volume->getXSize();
	size_t y_size = (size_t)volume->getYSize();
	size_t z_size = (size_t)volume->getZSize();
	
	dimensionSizes[0] = x_size;
	dimensionSizes[1] = y_size;
	dimensionSizes[2] = z_size;
	dimensionSizes[3] = 0;
	
	err = MDMakeWave(&DPWave, waveName.c_str(), NULL, dimensionSizes, NT_I16 | NT_UNSIGNED, 1);
	if (err != 0) {
		throw err;
	}
	
	for (size_t k = 0; k < z_size; ++k) {
		for (size_t j = 0; j < y_size; ++j) {
			for (size_t i = 0; i < x_size; ++i) {
				indices[0] = i;
				indices[1] = j;
				indices[2] = k;
				
				value[0] = (double)(*volume)(i, j, k);
				
				err = MDSetNumericWavePointValue(DPWave, indices, value);
				if (err != 0) {
					throw err;
				}
			}
		}
	}
	
	return DPWave;
}