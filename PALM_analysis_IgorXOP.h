/*
 *  PALM_analysis_IgorXOP.h
 *  PALM analysis
 *
 *  Created by Peter Dedecker on 05/01/2008.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef PALM_ANALYSIS_IGOR_XOP
#define PALM_ANALYSIS_IGOR_XOP

#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <sstream>
#include "boost/smart_ptr.hpp"
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




HOST_IMPORT int main(IORecHandle ioRecHandle);

boost::shared_ptr<ImageLoader> get_image_loader_for_camera_type(size_t camera_type, std::string data_file_path);

#endif

