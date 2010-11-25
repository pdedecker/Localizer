/*
 *  PALM_analysis_IgorUtilities.cpp
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 30/12/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "PALM_analysis_IgorUtilities.h"


int load_partial_ccd_image(ImageLoader *image_loader, size_t n_start, size_t n_end, DataFolderAndName destination) {
	size_t total_n_images = image_loader->GetNImages();
	size_t n_images_to_load;
	size_t x_size, y_size;
	int storage_type, waveType;
	waveHndl output_wave;
	long dimension_sizes[MAX_DIMENSIONS + 1];
	long indices[MAX_DIMENSIONS];
	int destWaveCreated;
	
	int result;
	boost::shared_ptr<ublas::matrix<double> > current_image;
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
	
	result = MDMakeWave(&output_wave, destination.name, destination.dfH, dimension_sizes, waveType, 1);
	if (result != 0)
		throw result;
	
	// load the data
	for (size_t i = n_start; i <= n_end; i++) {
		current_image = image_loader->readImage(i);
		
		indices[2] = i - n_start;
		
		// store the data in the output wave
		for (size_t j = 0; j < x_size; j++) {
			for (size_t k = 0; k < y_size; k++) {
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
	size_t total_n_images = image_loader->GetNImages();
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

waveHndl construct_summed_intensity_trace(ImageLoader *image_loader, DataFolderAndName outputWaveParams, long startX, long startY, long endX, long endY) {
	size_t n_images = image_loader->GetNImages();
	size_t x_size = image_loader->getXSize();
	size_t y_size = image_loader->getYSize();
	
	boost::shared_ptr<ublas::matrix<double> > current_image;
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
		if ((startX > endX) || (startY > endY))
			throw kBadROIDimensions;
	}
	
	// try to allocate a buffer that will hold the intensity trace
	boost::scoped_array<double> intensity_trace_buffer(new double[n_images]);
	
	for (size_t i = 0; i < n_images; i++) {
		summed_intensity = 0;
		current_image = image_loader->readImage(i);
		
		// calculate the total sum of the image
		for (size_t j = startX; j <= endX; j++) {
			for (size_t k = startY; k <= endY; k++) {
				summed_intensity += (*current_image)(j, k);
			}
		}
		// store the contents in the buffer
		intensity_trace_buffer[i] = summed_intensity;
	}
	
	// try to create the output wave
	result = MDMakeWave(&output_wave, outputWaveParams.name, outputWaveParams.dfH, dimension_sizes, NT_FP64, 1);
	if (result != 0)
		throw result;
	
	// write the output data to the wave
	result = MDStoreDPDataInNumericWave(output_wave, intensity_trace_buffer.get());
	if (result != 0) {
		throw result;
	}
	
	return output_wave;
}

waveHndl construct_average_intensity_trace(ImageLoader *image_loader, DataFolderAndName outputWaveParams, long startX, long startY, long endX, long endY) {
	
	int numDimensions;
	CountInt dimensionSizes[MAX_DIMENSIONS+1];
	IndexInt indices[MAX_DIMENSIONS];
	double value[2];
	int result;
	
	// start by calculating the summed intensity
	waveHndl outputWave = construct_summed_intensity_trace(image_loader, outputWaveParams, startX, startY, endX, endY);
	
	MDGetWaveDimensions(outputWave, &numDimensions, dimensionSizes);
	
	size_t nFrames = dimensionSizes[0];
	int xSize = endX - startX + 1;
	int ySize = endY - startY + 1;
	double nPixels = xSize * ySize;
	
	for (size_t i = 0; i < nFrames; ++i) {
		indices[0] = i;
		result = MDGetNumericWavePointValue(outputWave, indices, value);
		if (result != 0) {
			throw result;
		}
		value[0] /= nPixels;
		MDSetNumericWavePointValue(outputWave, indices, value);
	}
	
	return outputWave;
	
}

waveHndl construct_average_image(ImageLoader *image_loader, DataFolderAndName outputWaveParams, long startX, long startY, long endX, long endY) {
	size_t n_images = image_loader->GetNImages();
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
	
	boost::shared_ptr<ublas::matrix<double> > current_image;
	boost::shared_ptr<ublas::matrix<double> > average_image(new ublas::matrix<double>(xRange, yRange));
	
	waveHndl output_wave;
	long dimension_sizes[MAX_DIMENSIONS + 1];
	dimension_sizes[0] = xRange;
	dimension_sizes[1] = yRange;
	dimension_sizes[2] = 0;
	long indices[MAX_DIMENSIONS];
	double current_value[2];
	int result;
	
	std::fill(average_image->data().begin(), average_image->data().end(), double(0.0));
	
	
	for (size_t i = 0; i < n_images; i++) {
		current_image = image_loader->readImage(i);
		
		// add the values of the newly loaded image to the average image
		(*average_image) += (*current_image);
	}
	
	// divide by the number of images
	*average_image /= (double)n_images;
	
	// try to create the output wave
	result = MDMakeWave(&output_wave, outputWaveParams.name, outputWaveParams.dfH, dimension_sizes, NT_FP64, 1);
	if (result != 0)
		throw result;
	
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
	
	return output_wave;
}


waveHndl calculateStandardDeviationImage(ImageLoader *image_loader, DataFolderAndName outputWaveParams, long startX, long startY, long endX, long endY) {
	size_t n_images = image_loader->GetNImages();
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
	
	boost::scoped_ptr<ublas::matrix<double> > stdDevImage(new ublas::matrix<double>(xRange, yRange));
	boost::scoped_ptr<ublas::matrix<double> > average_image(new ublas::matrix<double>(xRange, yRange));
	boost::shared_ptr<ublas::matrix<double> > current_image;
	
	std::fill(average_image->data().begin(), average_image->data().end(), double(0.0));
	std::fill(stdDevImage->data().begin(), stdDevImage->data().end(), double(0.0));
	
	// construct an average image
	for (size_t i = 0; i < n_images; i++) {
		current_image = image_loader->readImage(i);
		
		// add the values of the newly loaded image to the average image
		(*average_image) += (*current_image);
	}
	
	// divide by the number of images
	*average_image /= (double)n_images;
	
	// now loop over the images again, calculating the standard deviation of each pixel
	for (size_t i = 0; i < n_images; i++) {
		current_image = image_loader->readImage(i);
		
		// add the deviation of the newly loaded image from the mean to the stddev image
		(*stdDevImage) = (*stdDevImage) + element_prod((*current_image) - (*average_image), (*current_image) - (*average_image));
	}
	
	// divide by the number of images to get the average deviation, and take the square root
	*stdDevImage /= (double)n_images;
	std::transform(stdDevImage->data().begin(), stdDevImage->data().end(), stdDevImage->data().begin(), static_cast<double(*)(double)>(sqrt));
	//(*stdDevImage) = (*stdDevImage).RaiseToPower(0.5);
	
	// try to create the output wave
	dimension_sizes[0] = xRange;
	dimension_sizes[1] = yRange;
	dimension_sizes[2] = 0;
	result = MDMakeWave(&output_wave, outputWaveParams.name, outputWaveParams.dfH, dimension_sizes, NT_FP64, 1);
	if (result != 0)
		throw result;
	
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
	
	return output_wave;
}

waveHndl FetchWaveUsingFullPath(std::string wavePath) {
	waveHndl fetchedWave;
	DataFolderHandle dataFolder;
	size_t wavePathOffset, position;
	std::string dataFolderPath;
	int err;
	
	// if the wavePath really is just a name, do not do any advanced processing
	if (wavePath.find(':') == std::string::npos) {
		fetchedWave = FetchWaveFromDataFolder(NULL, wavePath.c_str());
		if (fetchedWave == NULL) {
			throw NOWAV;
		}
	} else {
		// we found a ':' in the string
		// therefore we're being passed a wave location that also includes datafolder information
		// we need to parse this string to get out the data folder handle and the wave handle
		wavePathOffset = wavePath.rfind(':');
		if (wavePathOffset != 0) {	// check for the case where the user specifies something like ":wavePath"
			dataFolderPath = wavePath.substr(0, wavePathOffset);
		} else {
			dataFolderPath = ":";
		}
		err = GetNamedDataFolder(NULL, dataFolderPath.c_str(), &dataFolder);
		if (err != 0)
			throw err;
		
		// retain only the wavePath part, discard the information about the data folder
		wavePath = wavePath.substr(wavePathOffset + 1, wavePath.length() - wavePathOffset - 1);
		
		// if the wavePath part is quoted ('') then Igor chokes on this
		// so we remove the quotes
		while ((position = wavePath.find('\'')) != std::string::npos) {
			wavePath.erase(position, 1);
		}
		
		fetchedWave = FetchWaveFromDataFolder(dataFolder, wavePath.c_str());
		if (fetchedWave == NULL) {
			throw NOWAV;
		}
	}
	
	return fetchedWave;
}

waveHndl MakeWaveUsingFullPath(std::string wavePath, long *dimensionSizes, int type, int overwrite) {
	waveHndl createdWave;
	DataFolderHandle dataFolder;
	size_t wavePathOffset, position;
	std::string dataFolderPath;
	int err;
	
	// if the wavePath really is just a name, do not do any advanced processing
	if (wavePath.find(':') == std::string::npos) {
		err = MDMakeWave(&createdWave, wavePath.c_str(), NULL, dimensionSizes, type, overwrite);
		if (err != 0) {
			throw err;
		}
	} else {
		// we found a ':' in the string
		// therefore we're being passed a wave location that also includes datafolder information
		// we need to parse this string to get out the data folder handle and the wave handle
		wavePathOffset = wavePath.rfind(':');
		if (wavePathOffset != 0) {	// check for the case where the user specifies something like ":wavePath"
			dataFolderPath = wavePath.substr(0, wavePathOffset);
		} else {
			dataFolderPath = ":";
		}
		err = GetNamedDataFolder(NULL, dataFolderPath.c_str(), &dataFolder);
		if (err != 0)
			throw err;
		
		// retain only the wavePath part, discard the information about the data folder
		wavePath = wavePath.substr(wavePathOffset + 1, wavePath.length() - wavePathOffset - 1);
		
		// if the wavePath part is quoted ('') then Igor chokes on this
		// so we remove the quotes
		while ((position = wavePath.find('\'')) != std::string::npos) {
			wavePath.erase(position, 1);
		}
		
		err = MDMakeWave(&createdWave, wavePath.c_str(), dataFolder, dimensionSizes, type, overwrite);
		if (err != 0) {
			throw err;
		}
	}
	
	return createdWave;
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

int ConvertHandleToFilepathString(Handle handle, std::string &output_path) {
	int err;
	char handle_char[1024];
	
	err = GetCStringFromHandle(handle, handle_char, 1023);
	if (err != 0) {
		return err;
	}
	
#ifdef MACIGOR
	err = WinToMacPath(handle_char);
	if (err != 0) {
		return err;
	}
	
	char posixPATH[MAX_PATH_LEN+1];
	
	err = HFSToPosixPath(handle_char, posixPATH, 0);
	if (err != 0) {
		return err;
	}
	output_path.assign(posixPATH);
#endif
#ifdef WINIGOR
	err = MacToWinPath(handle_char);
	if (err != 0) {
		return err;
	}
	output_path.assign(handle_char);
#endif
	
	return 0;
	
}

waveHndl CopyVectorToIgorDPWave(boost::shared_ptr<std::vector<double> > vec, std::string waveName) {
	waveHndl DPWave;
	
	int err;
	long indices[MAX_DIMENSIONS];
	long dimensionSizes[MAX_DIMENSIONS+1];
	double value[2];
	
	
	size_t nElements = vec->size();
	
	dimensionSizes[0] = nElements;
	dimensionSizes[1] = 0;
	dimensionSizes[2] = 0;
	
	DPWave = MakeWaveUsingFullPath(waveName, dimensionSizes, NT_FP64, 1);
	
	indices[1] = 0;
	for (size_t i = 0; i < nElements; ++i) {
		indices[0] = i;
		
		value[0] = (*vec)[i];
		
		err = MDSetNumericWavePointValue(DPWave, indices, value);
		if (err != 0) {
			throw err;
		}
	}
	
	return DPWave;
}


boost::shared_ptr<ublas::matrix<double> > CopyIgorDPWaveToMatrix(waveHndl wave) {
	// copy a Igor wave into a new gsl_matrix
	
	int err;
	int numDimensions; 
	CountInt dimensionSizes[MAX_DIMENSIONS+1];
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
	
	boost::shared_ptr<ublas::matrix<double> > matrix(new ublas::matrix<double>(x_size, y_size));
	
	for (size_t i = 0; i < x_size; ++i) {
		for (size_t j = 0; j < y_size; ++j) {
			indices[0] = i;
			indices[1] = j;
			
			err = MDGetNumericWavePointValue(wave, indices, value);
			if (err != 0) {
				throw err;
			}
			
			(*matrix)(i, j) = value[0];
		}
	}
	
	return matrix;
}

waveHndl CopyMatrixToIgorDPWave(boost::shared_ptr<ublas::matrix<double> > matrix, std::string waveName) {
	
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
		
		DPWave = MakeWaveUsingFullPath(waveName, dimensionSizes, NT_FP64, 1);
		
		return DPWave;
		
	}
	
	
	size_t x_size = (size_t)matrix->size1();
	size_t y_size = (size_t)matrix->size2();
	
	dimensionSizes[0] = x_size;
	dimensionSizes[1] = y_size;
	dimensionSizes[2] = 0;
	
	DPWave = MakeWaveUsingFullPath(waveName, dimensionSizes, NT_FP64, 1);
	
	for (size_t i = 0; i < x_size; ++i) {
		for (size_t j = 0; j < y_size; ++j) {
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


