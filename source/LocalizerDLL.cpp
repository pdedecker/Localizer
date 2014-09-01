/*
 Copyright 2008-2014 Peter Dedecker.
 
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

#include "LocalizerDLL.h"

#include <memory>
#include <gsl/gsl_errno.h>

#include "FileIO.h"

int LocalizerFileInfo(char* filePath,            // in: path to the file on disk
                      int* nImagesInFile,        // out: number of images
                      int* nRows,                // out: number of rows
                      int* nCols,                // out: number of cols
                      int* storageType) {        // out: pixel storage type
    gsl_set_error_handler_off();	// we will handle GSL errors ourselves
    TIFFSetErrorHandler(NULL);      // we will handle libtiff errors ourselves
    try {
        std::shared_ptr<ImageLoader> imageLoader = GetImageLoader(std::string(filePath));
        *nImagesInFile = imageLoader->getNImages();
        *nRows = imageLoader->getXSize();
        *nCols = imageLoader->getYSize();
        *storageType = imageLoader->getStorageType();
    }
    catch (...) {
        return -1;
    }
    
    return 0;
}

EXPORT int LocalizerLoadImages(char* filePath,          // in: path to the file on disk
                               int nImagesToSkip,       // in: number of images to skip
                               int nImagesToLoad,       // in: number of images to load -- -1 to load up to end
                               double* imageData) {     // in: pointer to allocate buffer -- must contain at least (nRows * nCols * nImages) doubles
    gsl_set_error_handler_off();	// we will handle GSL errors ourselves
    TIFFSetErrorHandler(NULL);      // we will handle libtiff errors ourselves
    try {
        std::shared_ptr<ImageLoader> imageLoader = GetImageLoader(std::string(filePath));
        ImageLoaderWrapper imageLoaderWrapper(imageLoader);
        imageLoaderWrapper.setImageRange(nImagesToSkip, nImagesToLoad);
        int nRows = imageLoaderWrapper.getXSize();
        int nCols = imageLoaderWrapper.getYSize();
        int nImagesTotal = imageLoaderWrapper.getNImages();
        int nPixelsPerImage = nRows * nCols;
        int nBytesPerImage = nPixelsPerImage * sizeof(double);
        if (nImagesToLoad < 0) {
            nImagesToLoad = nImagesTotal - nImagesToSkip;
        }
        
        for (int n = 0; n < nImagesToLoad; ++n) {
            ImagePtr image = imageLoaderWrapper.readImage(n);
            memcpy(imageData, image->data(), nBytesPerImage);
            imageData += nPixelsPerImage;
        }
    }
    catch (...) {
        return -1;
    }
    
    return nImagesToLoad;
}
