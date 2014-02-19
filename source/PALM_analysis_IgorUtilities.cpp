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

    std::string extension = filePath.substr(startOfExtension + 1);
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
	if (boost::algorithm::iequals(extension, "btf"))
        return CAMERA_TYPE_TIFF;
	if (boost::algorithm::iequals(extension, "tf8"))
        return CAMERA_TYPE_TIFF;
    if (boost::algorithm::iequals(extension, "lsm"))
        return CAMERA_TYPE_TIFF;

    // if we're still here then the extension was not recognized
    throw std::runtime_error("Unable to deduce the file type");
}

std::shared_ptr<ImageLoader> GetImageLoader(size_t camera_type, std::string& data_file_path) {
    std::shared_ptr<ImageLoader> image_loader;
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
			image_loader = std::shared_ptr<ImageLoader>(new ImageLoaderSPE(convertedFilePath));
			break;
		case CAMERA_TYPE_ANDOR:
			image_loader = std::shared_ptr<ImageLoader>(new ImageLoaderAndor(convertedFilePath));
			break;
		case CAMERA_TYPE_HAMAMATSU:
			image_loader = std::shared_ptr<ImageLoader>(new ImageLoaderHamamatsu(convertedFilePath));
			break;
		case CAMERA_TYPE_TIFF:
			image_loader = std::shared_ptr<ImageLoader>(new ImageLoaderTIFF(convertedFilePath));
			break;
		case CAMERA_TYPE_PDE:
			image_loader = std::shared_ptr<ImageLoader>(new ImageLoaderPDE(convertedFilePath));
			break;
		case CAMERA_TYPE_ZEISS:	// Zeiss lsm files
			image_loader = std::shared_ptr<ImageLoader>(new ImageLoaderTIFF(convertedFilePath));
			break;
		case CAMERA_TYPE_IGOR_WAVE: // Matrix wave in Igor
			image_loader = std::shared_ptr<ImageLoader>(new ImageLoaderIgor(convertedFilePath));
			break;
		case CAMERA_TYPE_MULTIFILE_TIFF:
			image_loader = std::shared_ptr<ImageLoader>(new ImageLoaderMultiFileTIFF(convertedFilePath));
			break;
		default:
			throw std::runtime_error("Unsupported CCD file type (/Y flag)");
			break;
    }

    return image_loader;

}

int LoadPartialCCDImage(ImageLoader *image_loader, int firstImage, int nImagesRequested, int overwrite, DataFolderAndName destination,
                        std::shared_ptr<ProgressReporter> progressReporter) {
    int nImages = image_loader->getNImages();
    int maxNImagesToLoad, nImagesToLoad;
    int x_size, y_size;
    int storage_type;

    int result;
    ImagePtr current_image;

    if (firstImage >= nImages) {
        throw std::runtime_error("the requested starting image is larger than the number of images available in the file");
    }

    // how many images should we load?
    maxNImagesToLoad = nImages - firstImage;
    if (nImagesRequested == -1) {
        nImagesToLoad = maxNImagesToLoad;
    } else {
        if (nImagesRequested > maxNImagesToLoad) {
            nImagesToLoad = maxNImagesToLoad;
        } else {
            nImagesToLoad = nImagesRequested;
        }
    }

    // only report progress if more than 50 frames are requested
    int doProgress = 0, progressStatus = 0;
    if (nImagesToLoad >= 50)
        doProgress = 1;

    if (doProgress)
        progressReporter->CalculationStarted();

    x_size = image_loader->getXSize();
    y_size = image_loader->getYSize();
    storage_type = image_loader->getStorageType();

    result = SetOperationNumVar("V_numberOfImages", nImages);
    if (result != 0)
        throw int(result);
    result = SetOperationNumVar("V_xSize", x_size);
    if (result != 0)
        throw int(result);
    result = SetOperationNumVar("V_ySize", y_size);
    if (result != 0)
        throw int(result);
    result = SetOperationNumVar("V_storageType", storage_type);
    if (result != 0)
        throw int(result);
    result = SetOperationNumVar("V_firstImageLoaded", firstImage);
    if (result != 0)
        throw int(result);
    result = SetOperationNumVar("V_lastImageLoaded", firstImage + nImagesToLoad - 1);
    if (result != 0)
        throw int(result);
	result = SetOperationNumVar("V_fileType", image_loader->getFileType());
	if (result != 0)
		throw int(result);

    // allocate the object that will write the data to Igor
    IgorImageOutputWriter waveWriter(destination, nImagesToLoad, overwrite, storage_type);

    // load the data and write it to Igor
    int nextFrameRead;
    image_loader->spoolTo(firstImage);
    for (int i = firstImage; i < firstImage + nImagesToLoad; i++) {
        if (doProgress && (i % 10 == 0)) {
            progressStatus =  progressReporter->UpdateCalculationProgress(i - firstImage, nImagesToLoad);
            if (progressStatus != 0) {
                progressReporter->CalculationAborted();
                throw USER_ABORTED("");
            }
        }

        current_image = image_loader->readNextImage(nextFrameRead);
        assert(nextFrameRead == i);
        waveWriter.write_image(current_image);
        if (CheckAbort(0))
            return 0;
    }

    if (doProgress)
        progressReporter->CalculationDone();

    return 0;
}


int ParseCCDHeaders(ImageLoader *image_loader) {
    int total_n_images = image_loader->getNImages();
    int x_size = image_loader->getXSize();
    int y_size = image_loader->getYSize();
    int storageType = image_loader->getStorageType();

    int result;

    result = SetOperationNumVar("V_numberOfImages", total_n_images);
    if (result != 0)
        throw int(result);
    result = SetOperationNumVar("V_xSize", x_size);
    if (result != 0)
        throw int(result);
    result = SetOperationNumVar("V_ySize", y_size);
    if (result != 0)
        throw int(result);
    result = SetOperationNumVar("V_storageType", storageType);
    if (result != 0)
        throw int(result);
    result = SetOperationNumVar("V_firstImageLoaded", -1);
    if (result != 0)
        throw int(result);
    result = SetOperationNumVar("V_lastImageLoaded", -1);
    if (result != 0)
        throw int(result);
	result = SetOperationNumVar("V_fileType", image_loader->getFileType());
	if (result != 0)
		throw int(result);

    return 0;
}

std::vector<double> ConstructIntensityTrace(std::shared_ptr<ImageLoader> imageLoader, std::shared_ptr<ProgressReporter> progressReporter, bool doAverage) {
    int nImages = imageLoader->getNImages();
    int xSize = imageLoader->getXSize();
    int ySize = imageLoader->getYSize();
    int nPixelsInImage = xSize * ySize;
    
    double summedIntensity;
    std::vector<double> trace(nImages);
    
    progressReporter->CalculationStarted();
    int progressStatus;
    imageLoader->rewind();
    for (int i = 0; i < nImages; ++i) {
        if (i % 20 == 0) {
            progressStatus = progressReporter->UpdateCalculationProgress(i, nImages);
            if (progressStatus != 0) {
                progressReporter->CalculationAborted();
                throw USER_ABORTED("");
            }
        }
        
        ImagePtr currentImage = imageLoader->readNextImage();
        summedIntensity = 0.0;
        for (int k = 0; k < ySize; k++) {
            for (int j = 0; j < xSize; j++) {
                summedIntensity += (*currentImage)(j, k);
            }
        }
        trace[i] = summedIntensity;
    }
    
    if (doAverage) {
        for (int i = 0; i < nImages; ++i) {
            trace[i] /= static_cast<double>(nPixelsInImage);
        }
    }
    
    return trace;
}

std::vector<double> ConstructSummedIntensityTrace(std::shared_ptr<ImageLoader> imageLoader, std::shared_ptr<ProgressReporter> progressReporter) {
    return ConstructIntensityTrace(imageLoader, progressReporter, false);
}

std::vector<double> ConstructAverageIntensityTrace(std::shared_ptr<ImageLoader> imageLoader, std::shared_ptr<ProgressReporter> progressReporter) {
    return ConstructIntensityTrace(imageLoader, progressReporter, true);
}

ImagePtr ConstructAverageImage(std::shared_ptr<ImageLoader> imageLoader, std::shared_ptr<ProgressReporter> progressReporter) {
    int nImages = imageLoader->getNImages();
    int xSize = imageLoader->getXSize();
    int ySize = imageLoader->getYSize();

    ImagePtr current_image;
    ImagePtr average_image(new Image(xSize, ySize));

    average_image->setConstant(0.0);

    int progressStatus;
    progressReporter->CalculationStarted();
    int nextFrameRead;
    for (int i = 0; i < nImages; i++) {
        if (i % 20 == 0) {
            progressStatus = progressReporter->UpdateCalculationProgress(i, nImages);
            if (progressStatus != 0) {
                progressReporter->CalculationAborted();
                throw USER_ABORTED("");
            }
        }

        current_image = imageLoader->readNextImage(nextFrameRead);
        assert(i == nextFrameRead);
        
        (*average_image) += (*current_image);
    }

    *average_image /= static_cast<double>(nImages);

    progressReporter->CalculationDone();

    return average_image;
}


ImagePtr ConstructVarianceImage(std::shared_ptr<ImageLoader> imageLoader, std::shared_ptr<ProgressReporter> progressReporter) {
    int nImages = imageLoader->getNImages();
    int xSize = imageLoader->getXSize();
    int ySize = imageLoader->getYSize();

    std::shared_ptr<Image> varianceImage(new Image(xSize, ySize));
    boost::scoped_ptr<Image> average_image(new Image(xSize, ySize));
    ImagePtr current_image;

    average_image->setConstant(0.0);
    varianceImage->setConstant(0.0);

    progressReporter->CalculationStarted();
    int progressStatus;
    int nextFrameRead;
    // construct an average image
    for (int i = 0; i < nImages; i++) {
        if (i % 10 == 0) {
            progressStatus = progressReporter->UpdateCalculationProgress(i / 2, nImages);
            if (progressStatus != 0) {
                progressReporter->CalculationAborted();
                throw USER_ABORTED("");
            }
        }

        current_image = imageLoader->readNextImage(nextFrameRead);
        assert(nextFrameRead == i);

        (*average_image) += (*current_image);
    }

    *average_image /= static_cast<double>(nImages);

    // now loop over the images again, calculating the standard deviation of each pixel
    for (int i = 0; i < nImages; i++) {
        if (i % 10 == 0) {
            progressStatus = progressReporter->UpdateCalculationProgress((i + nImages) / 2, nImages);
            if (progressStatus != 0) {
                progressReporter->CalculationAborted();
                throw USER_ABORTED("");
            }
        }

        current_image = imageLoader->readImage(i);

        // add the deviation of the newly loaded image from the mean to the stddev image
        (*varianceImage) += ((*current_image) - (*average_image)).array().square().matrix();
    }

    *varianceImage /= static_cast<double>(nImages);

    progressReporter->CalculationDone();

    return varianceImage;
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

waveHndl MakeWaveUsingFullPath(std::string wavePath, CountInt *dimensionSizes, int type, int overwrite) {
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

    if (filePath.size() > MAX_PATH_LEN)
        throw PATH_TOO_LONG;

    char cString[MAX_PATH_LEN+1];
    memcpy(cString, filePath.c_str(), filePath.size() * sizeof(char));
    cString[filePath.size()] = '\0';



#ifdef MACIGOR
    err = WinToMacPath(cString);
    if (err != 0) {
        throw err;
    }

    char posixPATH[MAX_PATH_LEN+1];

    err = HFSToPosixPath(cString, posixPATH, 0);
    if (err != 0) {
        throw err;
    }
    convertedPath.assign(posixPATH);
#endif
#ifdef WINIGOR
    err = MacToWinPath(cString);
    if (err != 0) {
        throw err;
    }
    convertedPath.assign(cString);
#endif

    return convertedPath;

}

waveHndl CopyVectorToIgorDPWave(std::shared_ptr<std::vector<double> > vec, std::string waveName) {
    waveHndl DPWave;

    int err;
    CountInt indices[MAX_DIMENSIONS];
    CountInt dimensionSizes[MAX_DIMENSIONS+1];
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

waveHndl CopyVectorToIgorDPWave(const std::vector<double>& vec, DataFolderAndName outputWaveParams) {
    waveHndl wave;
    int nElements = vec.size();
    
    int err;
    CountInt dimensionSizes[MAX_DIMENSIONS+1];
    dimensionSizes[0] = nElements;
    dimensionSizes[1] = 0;
    
    err = MDMakeWave(&wave, outputWaveParams.name, outputWaveParams.dfH, dimensionSizes, NT_FP64, 1);
    if (err)
        throw err;
    
    void* waveDataPtr = WaveData(wave);
    memcpy(waveDataPtr, &(vec[0]), nElements * sizeof(double));
    return wave;
}


ImagePtr CopyIgorDPWaveToMatrix(waveHndl wave) {
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

    ImagePtr matrix(new Image((int)x_size, (int)y_size));

    err = MDGetDPDataFromNumericWave(wave, matrix->data());
    if (err != 0)
        throw err;

    return matrix;
}

waveHndl CopyMatrixToIgorDPWave(ImagePtr matrix, std::string waveName) {

    waveHndl DPWave;

    int err;
    CountInt dimensionSizes[MAX_DIMENSIONS+1];

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

waveHndl CopyMatrixToIgorDPWave(ImagePtr matrix, DataFolderAndName dataFolderAndName) {
	waveHndl outputWave;
	int err = 0;
	
	CountInt dimensionSizes[MAX_DIMENSIONS + 1];
	
	// special case:
    // if the matrix is NULL (such as when there are no positions found)
    // then we return an empty wave
    if (matrix.get() == NULL) {
        dimensionSizes[0] = 0;
		err = MDMakeWave(&outputWave, dataFolderAndName.name, dataFolderAndName.dfH, dimensionSizes, NT_FP64, 1);
		if (err != 0)
			throw err;
		
		
        return outputWave;
		
    }
	
	dimensionSizes[0] = matrix->rows();
	dimensionSizes[1] = matrix->cols();
	dimensionSizes[2] = 0;
	err = MDMakeWave(&outputWave, dataFolderAndName.name, dataFolderAndName.dfH, dimensionSizes, NT_FP64, 1);
	if (err != 0)
		throw err;
	
	err = MDStoreDPDataInNumericWave(outputWave, matrix->data());
	return outputWave;
}


