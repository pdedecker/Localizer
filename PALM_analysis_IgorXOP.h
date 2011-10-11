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
#include "boost/smart_ptr.hpp"
#include <eigen3/Eigen/Eigen>
#include "XOPStandardHeaders.h"
#include "PALM_analysis.h"
#include "PALM_analysis_defines.h"
#include "PALM_analysis_errors.h"
#include "PALM_analysis_Processing.h"
#include "PALM_analysis_FileIO.h"
#include "PALM_analysis_segmentation.h"
#include "PALM_analysis_Localization.h"
#include "PALM_analysis_ParticleFinding.h"
#include "PALM_analysis_PALMImages.h"
#include "PALM_analysis_IgorUtilities.h"
#include "PALM_analysis_PositionsProcessing.h"
#include "PALM_analysis_ProgressReporting.h"
#include "PALM_analysis_SOFI.h"


HOST_IMPORT int XOPMain(IORecHandle ioRecHandle);

boost::shared_ptr<ImageLoader> get_image_loader_for_camera_type(size_t camera_type, std::string data_file_path);

#endif

