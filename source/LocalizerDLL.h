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

#ifndef Localizer_LocalizerDLL_h
#define Localizer_LocalizerDLL_h

#include <string>

#include "boost/smart_ptr.hpp"

#include "PALM_analysis.h"

#define EXPORT __attribute__((visibility("default")))

EXPORT int DoLocalizationAnalysis(char *filePath,           // path to the file on disk
                                  double pfa,               // pfa for GLRT
                                  double psfWidth,          // width of the psf
                                  double** positionsArray,  // pointer to pointer to double that will be set to 2D array with results
                                  int* nPositions,          // number of positions in result
                                  int* nColumns);           // number of columns in the array

std::shared_ptr<ImageLoader> GetImageLoader(std::string filePath);

EXPORT void LocalizerFreeArray(double *arrayPtr);           // function to clean up arrays allocated by this program

#endif
