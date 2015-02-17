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

#include "IgorUtilities.h"

#include <algorithm>
#include <boost/algorithm/string.hpp>

#include "FileIO.h"
#include "ProgressReporting.h"

int LoadPartialCCDImage(ImageLoader *image_loader, int firstImage, int nImagesRequested, int overwrite, DataFolderAndName destination,
                        std::shared_ptr<ProgressReporter> progressReporter) {
    int nImages = image_loader->getNImages();
    int maxNImagesToLoad, nImagesToLoad;
    int x_size, y_size;

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
    LocalizerStorageType storageType = image_loader->getStorageType();

    result = SetOperationNumVar("V_numberOfImages", nImages);
    if (result != 0)
        throw int(result);
    result = SetOperationNumVar("V_xSize", x_size);
    if (result != 0)
        throw int(result);
    result = SetOperationNumVar("V_ySize", y_size);
    if (result != 0)
        throw int(result);
    result = SetOperationNumVar("V_storageType", static_cast<int>(storageType));
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
    IgorImageOutputWriter waveWriter(destination, nImagesToLoad, overwrite, storageType);

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

void WriteImagesToDisk(std::shared_ptr<ImageLoader> imageLoader, std::shared_ptr<ImageOutputWriter> outputWriter, std::shared_ptr<ProgressReporter> progressReporter, int firstImage, int nImagesToWrite) {
    
    int nImages = imageLoader->getNImages();
    firstImage = Clip(firstImage, 0, nImages - 1);
    nImagesToWrite = (nImagesToWrite > 0) ? (std::min(nImages - firstImage, nImagesToWrite)) : nImages;
    
    progressReporter->CalculationStarted();
    imageLoader->spoolTo(firstImage);
    for (int n = 0; n < nImagesToWrite; ++n) {
        if ((n % 20) == 0) {
            int doAbort = progressReporter->UpdateCalculationProgress(n + 1, nImagesToWrite);
            if (doAbort)
                throw USER_ABORTED("user abort");
        }
        ImagePtr image = imageLoader->readNextImage();
        outputWriter->write_image(image);
    }
    progressReporter->CalculationDone();
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
    std::unique_ptr<Image> average_image(new Image(xSize, ySize));
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
    std::unique_ptr<char[]> cString(new char[stringLength + 1]);

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

waveHndl CopyVectorToIgorDPWave(const std::vector<double>& vec, std::string waveName) {
    waveHndl DPWave;

    int err;
    CountInt indices[MAX_DIMENSIONS];
    CountInt dimensionSizes[MAX_DIMENSIONS+1];
    double value[2];


    size_t nElements = vec.size();

    dimensionSizes[0] = nElements;
    dimensionSizes[1] = 0;
    dimensionSizes[2] = 0;

    DPWave = MakeWaveUsingFullPath(waveName, dimensionSizes, NT_FP64, 1);

    indices[1] = 0;
    for (size_t i = 0; i < nElements; ++i) {
        indices[0] = i;

        value[0] = vec[i];

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

waveHndl CopyStackToIgorDPWave(std::vector<ImagePtr> stack, DataFolderAndName dataFolderAndName) {
    IgorImageOutputWriter outputWriter(dataFolderAndName, stack.size(), 1, kFP64);
    for (auto it = stack.cbegin(); it != stack.cend(); ++it) {
        outputWriter.write_image(*it);
    }
    
    return outputWriter.getWave();
}

void PrintToHistory(const std::string& str) {
    size_t offset = 0;
    while (offset < str.size()) {
        std::string subStr = str.substr(offset, MAXCMDLEN - 2);
        XOPNotice(subStr.c_str());
        offset += subStr.size();
    }
    XOPNotice("\r");
}
