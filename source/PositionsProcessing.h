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
 with analysis programs such as Igor Pro or Matlab, or to acquire
 data from or control hardware related to an experimental measurement,
 the licensors of this Program grant you additional permission
 to convey the resulting work.
 */

/*** routines that deal with processing localization results to produce some output ***/
#include <vector>
#include <cmath>
#include "Storage.h"
#include "PALMAnalysis.h"
#include "Defines.h"

// the visual studio compiler defines min and max preprocessor macros, but gcc doesn't
// so we have to undefine those for compilation on windows since otherwise the compiler
// throws an error on std::min and std::max
#undef min
#undef max

/**
 * Given a set of input positions, calculate the K-function or pairwise correlation to analyze clustering
 *
 * The K-function code is adapted from the implementation provided in the 'spatial' package for R,
 * based on Venables, W. N. and Ripley, B. D. (2002) Modern Applied Statistics with S. Fourth edition. Springer.
 * Their code is provided below with only very minor changes (make xl0 etc. nonglobal).
 */

std::vector<double> CalculateClustering(bool useKFunction, std::shared_ptr<LocalizedPositionsContainer> positions,
                                        double calculationRange, size_t nBins, double lowerX, double upperX,
                                        double lowerY, double upperY,
                                        std::shared_ptr<LocalizedPositionsContainer> positions2 = std::shared_ptr<LocalizedPositionsContainer>());

/**
 * The code to calculate the L function, copied with small modifications from the R spatial package by Ripley et al
 *
 *	xCoordinates1       array containing the x coordinates of probe 1
 *	yCoordinates1       array containing the y coordinates of probe 1
 *  xCoordinates2       array containing the x coordinates of probe 2   (pass xCoordinates1 for non-bivariate)
 *	yCoordinates2       array containing the y coordinates of probe 2   (pass yCoordinates1 for non-bivariate)
 *	nPoints1            number of coordinates in 1
 *	nPoints2            number of coordinates in 2
 *  nBins               requested number of output bins
 *  calculationRange    maximum distance between points considered
 *  lowerX
 *  upperX              coordinates of bounding rectangle
 *  lowerY
 *  upperY
 *  isKFunction         calculates K function is true, pairwise correlation otherwise
 */

std::vector<double> VR_sp_pp2(const double *xCoordinates1, const double *yCoordinates1, const double* xCoordinates2, const double* yCoordinates2,
                              size_t nPoints1, size_t nPoints2, size_t nBins,
                              double calculationRange, double upperX, double lowerX,
                              double upperY, double lowerY, bool isKFunction);


