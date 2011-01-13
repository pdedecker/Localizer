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
#include <Eigen/Eigen>
#include "PALM_analysis_storage.h"
#include "PALM_analysis.h"
#include "PALM_analysis_defines.h"

// the visual studio compiler defines min and max preprocessor macros, but gcc doesn't
// so we have to undefine those for compilation on windows since otherwise the compiler
// throws an error on std::min and std::max
#undef min
#undef max

/**
 * Given a set of input positions, calculate the L-function to analyze clustering
 *
 * This function is a wrapper around the implementation provided in the 'spatial' package for R,
 * based on Venables, W. N. and Ripley, B. D. (2002) Modern Applied Statistics with S. Fourth edition. Springer.
 * Their code is provided below with only very minor changes (make xl0 etc. nonglobal).
 */

boost::shared_ptr<std::vector<double> > CalculateLFunctionClustering(boost::shared_ptr<LocalizedPositionsContainer> positions,
																	 double plotFullScale, size_t nBins, double lowerX, double upperX,
																	 double lowerY, double upperY);

/**
 * The code to calculate the L function, verbatim from the R spatial package by Ripley et al
 *
 *	x	double array containing the x coordinates
 *	y	double array containing the y coordinates
 *	npt	the number of positions
 *	k	the number of points in the output (bins)
 *	h	double array that will be filled with the l function (allocated to k values)
 *	fs	the scale of the plot, that is, the range of the x coordinates of the l function
 *	xu0	the coordinates of the rectangle bounding the points
 */
void VR_sp_pp2(double *x, double *y, size_t *npt, size_t *k,
			   double *h, double *fs, double xu0, double xl0,
			   double yu0, double yl0);

/**
 * Function that corrects for edge effects in the calculated L function, verbatim from Ripley et al
 */
double VR_edge(double x, double y, double a, double xu0, double xl0,
			double yu0, double yl0);
