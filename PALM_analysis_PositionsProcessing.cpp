/*
 *  PALM_analysis_PositionsProcessing.cpp
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 06/01/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "PALM_analysis_PositionsProcessing.h"

boost::shared_ptr<std::vector<double> > CalculateLFunctionClustering(boost::shared_ptr<LocalizedPositionsContainer> positions,
																	 double calculationRange, size_t nBins, double lowerX, double upperX,
																	 double lowerY, double upperY) {
	size_t nPositions = positions->getNPositions();
	double currentX, currentY;
	
	// get the x and y boundaries for the positions if they have not been provided
	if ((lowerX == 0) && (upperX == 0) && (lowerY == 0) && (lowerX == 0)) {
		lowerX = 1e200;
		upperX = -1;
		lowerY = 1e200; 
		upperY = -1;
		
		for (size_t i = 0; i < nPositions; ++i) {
			currentX = positions->getXPosition(i);
			currentY = positions->getYPosition(i);
			
			if (currentX < lowerX)
				lowerX = currentX;
			if (currentX > upperX)
				upperX = currentX;
			if (currentY < lowerY)
				lowerY = currentY;
			if (currentY > upperY)
				upperY = currentY;
		}
	}
	
	// provide that positions as arrays suitable for the actual calculation
	boost::scoped_array<double> xPositions(new double[nPositions]);
	boost::scoped_array<double> yPositions(new double[nPositions]);
	boost::shared_ptr<std::vector<double> > lFunction(new std::vector<double>(nBins));
	
	for (size_t i = 0; i < nPositions; ++i) {
		xPositions[i] = positions->getXPosition(i);
		yPositions[i] = positions->getYPosition(i);
	}
	
	// run the actual calculation
	VR_sp_pp2(xPositions.get(), yPositions.get(), nPositions, &nBins,
			  &((*lFunction)[0]), calculationRange, upperX, lowerX,
			  upperY, lowerY);
	
    // return the calculated function
	return lFunction;
	
}

void VR_sp_pp2(double *xCoordinates, double *yCoordinates, size_t nPoints, size_t *nBins,
			   double *outputArray, double calculationRange, double upperX, double lowerX,
			   double upperY, double lowerY) {
    int   nIncludedBins, i, j, ib;
    double xSize, ySize, xi, yi, sqrtArea, g, dm, alm;
    double sqDistance, x1, y1, sqMaxDistance, effectiveCalculationRange, binsPerDistance;
	
    // testinit();
    xSize = upperX - lowerX;
    ySize = upperY - lowerY;
    sqrtArea = sqrt(xSize * ySize);
    dm = calculationRange;
    g = 2.0 / (nPoints * nPoints);
    effectiveCalculationRange = std::min(calculationRange, 0.5 * sqrt(xSize * xSize + ySize * ySize));
    binsPerDistance = *nBins / calculationRange;
    nIncludedBins = floor(binsPerDistance * effectiveCalculationRange + 1e-3);
    *nBins = nIncludedBins;
    sqMaxDistance = effectiveCalculationRange * effectiveCalculationRange;
    for (i = 0; i < *nBins; i++)
		outputArray[i] = 0.0;
	
    for (i = 1; i < nPoints; i++) {
		xi = xCoordinates[i];
		yi = yCoordinates[i];
		for (j = 0; j < i; j++) {
			x1 = xCoordinates[j] - xi;
			y1 = yCoordinates[j] - yi;
			sqDistance = x1 * x1 + y1 * y1;
			if (sqDistance < sqMaxDistance) {
				sqDistance = sqrt(sqDistance);
				dm = std::min(sqDistance, dm);
				ib = floor(binsPerDistance * sqDistance);
				if (ib < nIncludedBins)
					outputArray[ib] += g * (VR_edge(xi, yi, sqDistance, upperX, lowerX, upperY, lowerY) + VR_edge(xCoordinates[j], yCoordinates[j], sqDistance, upperX, lowerX, upperY, lowerY));
			}
		}
    }
    sqDistance = 0.0;
    alm = 0.0;
    for (i = 0; i < nIncludedBins; i++) {
		sqDistance += outputArray[i];
		outputArray[i] = sqrt(sqDistance / M_PI) * sqrtArea;
    }
}

double VR_edge(double x, double y, double a, double xu0, double xl0,
			double yu0, double yl0) {
    double b, c, c1, c2, r[6], w;
    int   i;
	
    w = x - xl0;
    if (w > y - yl0) w = y - yl0;
    if (w > xu0 - x) w = xu0 - x;
    if (w > yu0 - y) w = yu0 - y;
    if (a <= w) return (0.5);
    r[4] = r[0] = x - xl0;
    r[5] = r[1] = yu0 - y;
    r[2] = xu0 - x;
    r[3] = y - yl0;
    b = 0.0;
    for (i = 1; i <= 4; i++)
		if (r[i] < a) {
			if (r[i] == 0.0)
				b += M_PI;
			else {
				c = acos(r[i] / a);
				c1 = atan(r[i - 1] / r[i]);
				c2 = atan(r[i + 1] / r[i]);
				b += std::min(c, c1);
				b += std::min(c, c2);
			}
		}
    if (b < 6.28)
		return (1.0 / (2.0 - b / M_PI));
    return (0.0);
}
