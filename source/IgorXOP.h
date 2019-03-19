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

#ifndef PALM_ANALYSIS_IGOR_XOP
#define PALM_ANALYSIS_IGOR_XOP

#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <sstream>

#include <eigen3/Eigen/Eigen>
#include "XOPStandardHeaders.h"
#include "PALMAnalysis.h"
#include "Defines.h"
#include "Errors.h"
#include "Processing.h"
#include "FileIO.h"
#include "segmentation.h"
#include "Localization.h"
#include "ParticleFinding.h"
#include "PALMImages.h"
#include "IgorUtilities.h"
#include "PositionsProcessing.h"
#include "ProgressReporting.h"
#include "SOFI.h"


HOST_IMPORT int XOPMain(IORecHandle ioRecHandle);

std::shared_ptr<ImageLoader> get_image_loader_for_camera_type(size_t camera_type, std::string data_file_path);

#endif

