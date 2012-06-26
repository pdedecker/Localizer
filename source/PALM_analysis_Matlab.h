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

#include <string>
#include "boost/algorithm/string.hpp"

#include "mex.h"

#include "PALM_analysis_defines.h"

class ImageLoader;

void mexFunction(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs);
void MatlabLocalization(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs);
void MatlabTestSegmentation(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs);

std::string GetMatlabString(const mxArray* array);
boost::shared_ptr<ImageLoader> GetImageLoader(std::string& data_file_path);
int GetFileStorageType(std::string &filePath);
mxArray* ConvertImageToArray(ImagePtr image);
