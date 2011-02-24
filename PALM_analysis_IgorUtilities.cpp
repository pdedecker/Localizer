/*
 *  PALM_analysis_IgorUtilities.cpp
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 30/12/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "PALM_analysis_IgorUtilities.h"

void GetFilePathAndCameraType(std::string &inputFilePath, std::string &filePath, size_t &cameraType) {
	int isWavePath = 1;
	
	try {
		FetchWaveUsingFullPath(inputFilePath);
	}
	// if the wave pointed to by the path does not the exist
	// (as would be the case if it was a file path)
	// then FetchWaveUsingFullPath will throw exceptions
	catch (...) {
		isWavePath = 0;
	}
	
	if (isWavePath == 1) {
		filePath = inputFilePath;
		cameraType = CAMERA_TYPE_IGOR_WAVE;
		return;
	}
	
	// if we're still here then it's a path to a file
	// first we need to try to convert the path to the
	// appropriate format, if required
	std::string convertedPath = ConvertPathToNativePath(inputFilePath);
	
	filePath = convertedPath;
	cameraType = GetFileStorageType(filePath);
	return;
}

int GetFileStorageType(std::string &filePath) {
	size_t startOfExtension = filePath.rfind('.');
	if (startOfExtension == size_t(-1)) {
		// the filepath does not appear to contain an extension
		throw std::runtime_error("Unable to deduce the file type");
	}
	
	std::string extension = filePath.substr(startOfExtension);
	if ((extension.length() < 3) || (extension.length() > 4)) {
		throw std::runtime_error("Unable to deduce the file type");
	}
	
	if (boost::algorithm::iequals(extension, "spe"))
		return CAMERA_TYPE_WINSPEC;
	if (boost::algorithm::iequals(extension, "sif"))
		return CAMERA_TYPE_ANDOR;
	if (boost::algorithm::iequals(extension, "his"))
		return CAMERA_TYPE_HAMAMATSU;
	if (boost::algorithm::iequals(extension, "pde"))
		return CAMERA_TYPE_PDE;
	if (boost::algorithm::iequals(extension, "tif"))
		return CAMERA_TYPE_TIFF;
	if (boost::algorithm::iequals(extension, "tiff"))
		return CAMERA_TYPE_TIFF;
	if (boost::algorithm::iequals(extension, "lsm"))
		return CAMERA_TYPE_TIFF;
	
	// if we're still here then the extension was not recognized
	throw std::runtime_error("Unable to deduce the file type");
}

boost::shared_ptr<ImageLoader> GetImageLoader(size_t camera_type, std::string& data_file_path) {
	boost::shared_ptr<ImageLoader> image_loader;
	size_t estimatedCameraType;
	std::string convertedFilePath;
	
	// the camera type might be unknown
	if (camera_type == (size_t)-1) {
		GetFilePathAndCameraType(data_file_path, convertedFilePath, estimatedCameraType);
	} else {
		estimatedCameraType = camera_type;
		if (estimatedCameraType != CAMERA_TYPE_IGOR_WAVE) {
			convertedFilePath = ConvertPathToNativePath(data_file_path);
		} else {
			convertedFilePath = data_file_path;
		}
	}
	
	switch (estimatedCameraType) {
		case CAMERA_TYPE_WINSPEC:	// spe files
			image_loader = boost::shared_ptr<ImageLoader>(new ImageLoaderSPE(convertedFilePath));
			break;
		case CAMERA_TYPE_ANDOR:
			image_loader = boost::shared_ptr<ImageLoader>(new ImageLoaderAndor(convertedFilePath));
			break;
		case CAMERA_TYPE_HAMAMATSU:
			image_loader = boost::shared_ptr<ImageLoader>(new ImageLoaderHamamatsu(convertedFilePath));
			break;
		case CAMERA_TYPE_TIFF:	// 3 is reserved for TIFF files
			image_loader = boost::shared_ptr<ImageLoader>(new ImageLoaderTIFF(convertedFilePath));
			break;
		case CAMERA_TYPE_PDE:
			image_loader = boost::shared_ptr<ImageLoader>(new ImageLoaderPDE(convertedFilePath));
			break;
		case CAMERA_TYPE_ZEISS:	// Zeiss lsm files
			image_loader = boost::shared_ptr<ImageLoader>(new ImageLoaderTIFF(convertedFilePath));
			break;
		case CAMERA_TYPE_IGOR_WAVE: // Matrix wave in Igor
			image_loader = boost::shared_ptr<ImageLoader>(new ImageLoaderIgor(convertedFilePath));
			break;
		default:
			throw std::runtime_error("Unsupported CCD file type (/Y flag)");
			break;
	}
	
	return image_loader;
	
}

int load_partial_ccd_image(ImageLoader *image_loader, size_t n_start, size_t n_end, DataFolderAndName destination) {
	size_t total_n_images = image_loader->GetNImages();
	size_t n_images_to_load;
	size_t x_size, y_size;
	int storage_type;
	
	int result;
	boost::shared_ptr<Eigen::MatrixXd> current_image;
	
	if (n_start > n_end)
		throw END_SHOULD_BE_LARGER_THAN_START(std::string("When loading part of a CCD file the start index was larger than the end index"));
	
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
	
	// allocate the object that will write the data to Igor
	IgorImageOutputWriter waveWriter(destination, n_images_to_load, 1, storage_type);
	
	// load the data and write it to Igor
	for (size_t i = n_start; i <= n_end; i++) {
		current_image = image_loader->readImage(i);
		waveWriter.write_image(current_image);
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
	
	boost::shared_ptr<Eigen::MatrixXd> current_image;
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
		for (size_t k = startY; k <= endY; k++) {
			for (size_t j = startX; j <= endX; j++) {
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
	
	boost::shared_ptr<Eigen::MatrixXd> current_image;
	boost::shared_ptr<Eigen::MatrixXd> average_image(new Eigen::MatrixXd((int)xRange, (int)yRange));
	
	waveHndl output_wave;
	long dimension_sizes[MAX_DIMENSIONS + 1];
	dimension_sizes[0] = xRange;
	dimension_sizes[1] = yRange;
	dimension_sizes[2] = 0;
	int result;
	
	average_image->setConstant(0.0);
	
	
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
	result = MDStoreDPDataInNumericWave(output_wave, average_image->data());
	if (result != 0)
		throw result;
	
	return output_wave;
}


waveHndl calculateStandardDeviationImage(ImageLoader *image_loader, DataFolderAndName outputWaveParams, long startX, long startY, long endX, long endY) {
	size_t n_images = image_loader->GetNImages();
	size_t x_size = image_loader->getXSize();
	size_t y_size = image_loader->getYSize();
	int result;
	waveHndl output_wave;
	long dimension_sizes[MAX_DIMENSIONS + 1];
	
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
	
	boost::scoped_ptr<Eigen::MatrixXd> stdDevImage(new Eigen::MatrixXd((int)xRange, (int)yRange));
	boost::scoped_ptr<Eigen::MatrixXd> average_image(new Eigen::MatrixXd((int)xRange, (int)yRange));
	boost::shared_ptr<Eigen::MatrixXd> current_image;
	
	average_image->setConstant(0.0);
	stdDevImage->setConstant(0.0);
	
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
		(*stdDevImage) += ((*current_image) - (*average_image)).cwise().square();
	}
	
	// divide by the number of images to get the average deviation, and take the square root
	*stdDevImage /= (double)n_images;
	*stdDevImage = (*stdDevImage).cwise().sqrt();
	
	// try to create the output wave
	dimension_sizes[0] = xRange;
	dimension_sizes[1] = yRange;
	dimension_sizes[2] = 0;
	result = MDMakeWave(&output_wave, outputWaveParams.name, outputWaveParams.dfH, dimension_sizes, NT_FP64, 1);
	if (result != 0)
		throw result;
	
	// write the output data to the wave
	result = MDStoreDPDataInNumericWave(output_wave, stdDevImage->data());
	if (result != 0)
		throw result;
	
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


std::string ConvertHandleToString(Handle handle) {
	int err;
	std::string convertedString;
	
	size_t stringLength = GetHandleSize(handle);
	boost::scoped_array<char> cString(new char[stringLength + 1]);
	
	err = GetCStringFromHandle(handle, cString.get(), stringLength);
	if (err != 0)
		throw err;
	
	// save the wavenote as a std::string
	convertedString.assign(cString.get());
	
	return convertedString;
}

std::string ConvertPathToNativePath(std::string filePath) {
	int err;
	std::string convertedPath;
	
	boost::scoped_array<char> cString(new char[filePath.size() + 1]);
	memcpy(cString.get(), filePath.c_str(), filePath.size() * sizeof(char));
	cString[filePath.size()] = '\0';
	
	
	
	#ifdef MACIGOR
	err = WinToMacPath(cString.get());
	if (err != 0) {
		throw err;
	}
	
	char posixPATH[MAX_PATH_LEN+1];
	
	err = HFSToPosixPath(cString.get(), posixPATH, 0);
	if (err != 0) {
		throw err;
	}
	convertedPath.assign(posixPATH);
	#endif
	#ifdef WINIGOR
	err = MacToWinPath(cString.get());
	if (err != 0) {
		throw err;
	}
	convertedPath.assign(cString.get());
	#endif
	
	return convertedPath;
	
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


boost::shared_ptr<Eigen::MatrixXd> CopyIgorDPWaveToMatrix(waveHndl wave) {
	// copy a Igor wave into a new gsl_matrix
	
	int err;
	int numDimensions; 
	CountInt dimensionSizes[MAX_DIMENSIONS+1];
	size_t x_size, y_size;
	
	
	err = MDGetWaveDimensions(wave, &numDimensions, dimensionSizes);
	if (err != 0) {
		throw err;
	}
	if (numDimensions != 2) {
		throw INCOMPATIBLE_DIMENSIONING;
	}
	
	x_size = dimensionSizes[0];
	y_size = dimensionSizes[1];
	
	boost::shared_ptr<Eigen::MatrixXd> matrix(new Eigen::MatrixXd((int)x_size, (int)y_size));
	
	err = MDGetDPDataFromNumericWave(wave, matrix->data());
	if (err != 0)
		throw err;
	
	return matrix;
}

waveHndl CopyMatrixToIgorDPWave(boost::shared_ptr<Eigen::MatrixXd> matrix, std::string waveName) {
	
	waveHndl DPWave;
	
	int err;
	long dimensionSizes[MAX_DIMENSIONS+1];
	
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
	
	
	size_t x_size = (size_t)matrix->rows();
	size_t y_size = (size_t)matrix->cols();
	
	dimensionSizes[0] = x_size;
	dimensionSizes[1] = y_size;
	dimensionSizes[2] = 0;
	
	DPWave = MakeWaveUsingFullPath(waveName, dimensionSizes, NT_FP64, 1);
	
	err = MDStoreDPDataInNumericWave(DPWave, matrix->data());
	if (err != 0)
		throw err;
	
	return DPWave;
}


