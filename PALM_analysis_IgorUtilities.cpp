/*
 *  PALM_analysis_IgorUtilities.cpp
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 30/12/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "PALM_analysis_IgorUtilities.h"


int load_partial_ccd_image(ImageLoader *image_loader, size_t n_start, size_t n_end) {
	size_t total_n_images = image_loader->get_total_number_of_images();
	size_t n_images_to_load;
	size_t x_size, y_size;
	int storage_type, waveType;
	waveHndl output_wave;
	long dimension_sizes[MAX_DIMENSIONS + 1];
	long indices[MAX_DIMENSIONS];
	
	int result;
	boost::shared_ptr<PALMMatrix<double> > current_image;
	double current_value;
	double current_value_array[2];
	
	if (n_start > n_end)
		throw END_SHOULD_BE_LARGER_THAN_START(std::string("When loading part of a CCD file a the starting index was larger than the ending index"));
	
	// how many images should we load?
	if (n_end <= total_n_images) {
		n_images_to_load = n_end - n_start + 1;
	} else {
		n_images_to_load = total_n_images - n_start;
		n_end = total_n_images - 1;
	}
	
	x_size = image_loader->getXSize();
	y_size = image_loader->getYSize();
	storage_type = image_loader->getStorageType();
	
	result = SetOperationNumVar("V_numberOfImages", total_n_images);
	if (result != 0)
		throw result;
	result = SetOperationNumVar("V_xSize", x_size);
	if (result != 0)
		throw result;
	result = SetOperationNumVar("V_ySize", y_size);
	if (result != 0)
		throw result;
	result = SetOperationNumVar("V_firstImageLoaded", n_start);
	if (result != 0)
		throw result;
	result = SetOperationNumVar("V_lastImageLoaded", n_end);
	if (result != 0)
		throw result;
	
	
	// make the output wave that will store the data
	dimension_sizes[0] = (long)x_size;
	dimension_sizes[1] = (long)y_size;
	dimension_sizes[2] = (long)n_images_to_load;
	dimension_sizes[3] = 0;
	
	// decide on the storage type to use
	// we use the storage type that matches that of the original frames
	switch (storage_type) {
		case STORAGE_TYPE_INT4:
		case STORAGE_TYPE_UINT4:
		case STORAGE_TYPE_INT8:
			waveType = NT_I8;
			break;
		case STORAGE_TYPE_UINT8:
			waveType = NT_I8 | NT_UNSIGNED;
			break;
		case STORAGE_TYPE_INT16:
			waveType = NT_I16;
			break;
		case STORAGE_TYPE_UINT16:
			waveType = NT_I16 | NT_UNSIGNED;
			break;
		case STORAGE_TYPE_INT32:
			waveType = NT_I32;
			break;
		case STORAGE_TYPE_UINT32:
			waveType = NT_I32 | NT_UNSIGNED;
			break;
		case STORAGE_TYPE_FP32:
			waveType = NT_FP32;
			break;
		case STORAGE_TYPE_FP64:
			waveType = NT_FP64;
			break;
		default:
			waveType = NT_FP64;
	}
	
	result = MDMakeWave(&output_wave, "M_CCDFrames", NULL, dimension_sizes, waveType, 1);
	
	if (result != 0) {
		throw result;
	}
	
	// load the data
	for (size_t i = n_start; i <= n_end; i++) {
		current_image = image_loader->get_nth_image(i);
		
		indices[2] = i - n_start;
		
		// store the data in the output wave
		for (size_t k = 0; k < y_size; k++) {
			for (size_t j = 0; j < x_size; j++) {
				current_value = (*current_image)(j, k);
				
				indices[0] = j;
				indices[1] = k;
				
				current_value_array[0] = current_value;
				
				result = MDSetNumericWavePointValue(output_wave, indices, current_value_array);
				if (result != 0) {
					throw result;
				}
			}
		}
		
	}
	
	
	return 0;
}


int parse_ccd_headers(ImageLoader *image_loader) {
	size_t total_n_images = image_loader->get_total_number_of_images();
	size_t x_size = image_loader->getXSize();
	size_t y_size = image_loader->getYSize();
	
	int result;
	
	result = SetOperationNumVar("V_numberOfImages", total_n_images);
	if (result != 0)
		throw result;
	result = SetOperationNumVar("V_xSize", x_size);
	if (result != 0)
		throw result;
	result = SetOperationNumVar("V_ySize", y_size);
	if (result != 0)
		throw result;
	result = SetOperationNumVar("V_firstImageLoaded", -1);
	if (result != 0)
		throw result;
	result = SetOperationNumVar("V_lastImageLoaded", -1);
	if (result != 0)
		throw result;
	
	return 0;
}

int construct_summed_intensity_trace(ImageLoader *image_loader, string output_wave_name, long startX, long startY, long endX, long endY) {
	size_t n_images = image_loader->get_total_number_of_images();
	size_t x_size = image_loader->getXSize();
	size_t y_size = image_loader->getYSize();
	
	boost::shared_ptr<PALMMatrix<double> > current_image;
	double summed_intensity;
	
	waveHndl output_wave;
	long dimension_sizes[MAX_DIMENSIONS + 1];
	dimension_sizes[0] = (long)n_images;
	dimension_sizes[1] = 0;
	int result;
	
	// check if we want the full frame
	if (startX == -1) {
		startX = 0;
		startY = 0;
		endX = x_size - 1;
		endY = y_size - 1;
	} else {
		if ((startX >= x_size) || (endX >= x_size) || (startY >= y_size) || (endY >= y_size)) {
			throw kBadROIDimensions;
		}
	}
	
	// try to allocate a buffer that will hold the intensity trace
	boost::scoped_array<double> intensity_trace_buffer(new double[n_images]);
	
	for (size_t i = 0; i < n_images; i++) {
		summed_intensity = 0;
		current_image = image_loader->get_nth_image(i);
		
		// calculate the total sum of the image
		for (size_t k = startY; k <= endY; k++) {
			for (size_t j = startX; j <= endX; j++) {
				summed_intensity += (*current_image)(j, k);
			}
		}
		// store the contents in the buffer
		intensity_trace_buffer[i] = summed_intensity;
	}
	
	// try to create the output wave
	result = MDMakeWave(&output_wave, output_wave_name.c_str(), NULL, dimension_sizes, NT_FP64, 1);
	if (result != 0) {
		throw result;
	}
	
	// write the output data to the wave
	result = MDStoreDPDataInNumericWave(output_wave, intensity_trace_buffer.get());
	if (result != 0) {
		throw result;
	}
	
	return 0;
}

int construct_average_image(ImageLoader *image_loader, string output_wave_name, long startX, long startY, long endX, long endY) {
	size_t n_images = image_loader->get_total_number_of_images();
	size_t x_size = image_loader->getXSize();
	size_t y_size = image_loader->getYSize();
	
	long xRange, yRange;
	
	if (endX > x_size - 1)
		endX = x_size - 1;
	if (endY > y_size - 1)
		endY = y_size - 1;
	
	if ((startX < 0) || (startY < 0) || (startY < 0) || (startY < 0)) {	// we want the full frame
		startX = 0;
		startY = 0;
		endX = x_size - 1;
		endY = y_size - 1;
	}
	
	xRange = endX - startX + 1;
	yRange = endY - startY + 1;
	
	boost::shared_ptr<PALMMatrix<double> > current_image;
	boost::shared_ptr<PALMMatrix<double> > average_image(new PALMMatrix<double>(xRange, yRange));
	
	waveHndl output_wave;
	long dimension_sizes[MAX_DIMENSIONS + 1];
	dimension_sizes[0] = xRange;
	dimension_sizes[1] = yRange;
	dimension_sizes[2] = 0;
	long indices[MAX_DIMENSIONS];
	double current_value[2];
	int result;
	
	average_image->set_all(0);
	
	
	for (size_t i = 0; i < n_images; i++) {
		current_image = image_loader->get_nth_image(i);
		
		// add the values of the newly loaded image to the average image
		(*average_image) += (*current_image);
	}
	
	// divide by the number of images
	(*average_image) = (*average_image).DivideByScalar(n_images);
	
	// try to create the output wave
	result = MDMakeWave(&output_wave, output_wave_name.c_str(), NULL, dimension_sizes, NT_FP64, 1);
	if (result != 0) {
		throw result;
	}
	
	// write the output data to the wave
	for (size_t i = startX; i <= endX; ++i) {
		for (size_t j = startY; j <= endY; ++j) {
			current_value[0] = (*average_image)(i, j);
			indices[0] = i;
			indices[1] = j;
			result = MDSetNumericWavePointValue(output_wave, indices, current_value);
			if (result != 0) {
				throw result;
			}
		}
	}
	
	return 0;
}


void calculateStandardDeviationImage(ImageLoader *image_loader, string output_wave_name, long startX, long startY, long endX, long endY) {
	size_t n_images = image_loader->get_total_number_of_images();
	size_t x_size = image_loader->getXSize();
	size_t y_size = image_loader->getYSize();
	int result;
	waveHndl output_wave;
	long dimension_sizes[MAX_DIMENSIONS + 1];
	long indices[MAX_DIMENSIONS];
	
	double current_value[2];
	
	long xRange, yRange;
	
	if (endX > x_size - 1)
		endX = x_size - 1;
	if (endY > y_size - 1)
		endY = y_size - 1;
	
	if ((startX < 0) || (startY < 0) || (startY < 0) || (startY < 0)) {	// we want the full frame
		startX = 0;
		startY = 0;
		endX = x_size - 1;
		endY = y_size - 1;
	}
	
	xRange = endX - startX + 1;
	yRange = endY - startY + 1;
	
	boost::scoped_ptr<PALMMatrix<double> > stdDevImage(new PALMMatrix<double>(xRange, yRange));
	boost::scoped_ptr<PALMMatrix<double> > average_image(new PALMMatrix<double>(xRange, yRange));
	boost::shared_ptr<PALMMatrix<double> > current_image;
	
	average_image->set_all(0);
	stdDevImage->set_all(0);
	
	// construct an average image
	for (size_t i = 0; i < n_images; i++) {
		current_image = image_loader->get_nth_image(i);
		
		// add the values of the newly loaded image to the average image
		(*average_image) += (*current_image);
	}
	
	// divide by the number of images
	(*average_image) = (*average_image).DivideByScalar(n_images);
	
	// now loop over the images again, calculating the standard deviation of each pixel
	for (size_t i = 0; i < n_images; i++) {
		current_image = image_loader->get_nth_image(i);
		
		// add the deviation of the newly loaded image from the mean to the stddev image
		(*stdDevImage) += (((*current_image) - (*average_image)) * ((*current_image) - (*average_image)));
	}
	
	// divide by the number of images to get the average deviation, and take the square root
	(*stdDevImage) = (*stdDevImage).DivideByScalar(n_images);
	(*stdDevImage) = (*stdDevImage).RaiseToPower(0.5);
	
	// try to create the output wave
	dimension_sizes[0] = xRange;
	dimension_sizes[1] = yRange;
	dimension_sizes[2] = 0;
	result = MDMakeWave(&output_wave, output_wave_name.c_str(), NULL, dimension_sizes, NT_FP64, 1);
	if (result != 0) {
		throw result;
	}
	
	// write the output data to the wave
	for (size_t i = startX; i <= endX; ++i) {
		for (size_t j = startY; j <= endX; ++j) {
			current_value[0] = (*stdDevImage)(i, j);
			indices[0] = i;
			indices[1] = j;
			result = MDSetNumericWavePointValue(output_wave, indices, current_value);
			if (result != 0) {
				throw result;
			}
		}
	}
}


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