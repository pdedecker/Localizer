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

#ifndef Localizer_LocalizerDLL_h
#define Localizer_LocalizerDLL_h

#ifdef _WIN32
    #define EXPORT __declspec( dllexport )
#else
    #define EXPORT __attribute__((visibility("default")))
#endif

extern "C" {
    
    // returns non-zero for error
    EXPORT int LocalizerFileInfo(char* filePath,            // in: path to the file on disk
                                 int* nImagesInFile,        // out: number of images
                                 int* nRows,                // out: number of rows
                                 int* nCols,                // out: number of cols
                                 int* storageType);         // out: pixel storage type (data will always be loaded as doubles)
    
    // return number of images loaded, or negative value for error.
    EXPORT int LocalizerLoadImages(char* filePath,          // in: path to the file on disk
                                   int nImagesToSkip,       // in: number of images to skip
                                   int nImagesToLoad,       // in: number of images to load -- -1 to load up to end
                                   double* imageData);      // in: pointer to allocate buffer -- must contain at least (nRows * nCols * nImages) doubles
    
}
#endif
