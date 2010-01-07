/*
 *  PALM_analysis_PositionsProcessing.h
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 06/01/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


/*** routines that deal with processing localization results to produce some output ***/
#include <vector>
#include <cmath>
#include "boost/smart_ptr.hpp"
#include "PALM_analysis_storage.h"
#include "PALM_analysis.h"
#include "PALM_analysis_defines.h"


/**
 * Given a set of input positions, calculate the K-function to analyze clustering
*/

boost::shared_ptr<std::vector<double> > CalculateKFunctionClustering(boost::shared_ptr<LocalizedPositionsContainer> positions,
																	 double binWidth, size_t nBins);
