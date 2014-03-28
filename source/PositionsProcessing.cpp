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

#include "PALM_analysis_PositionsProcessing.h"

std::shared_ptr<std::vector<double> > CalculateLFunctionClustering(std::shared_ptr<LocalizedPositionsContainer> positions,
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
	std::unique_ptr<double[]> xPositions(new double[nPositions]);
	std::unique_ptr<double[]> yPositions(new double[nPositions]);
	std::shared_ptr<std::vector<double> > lFunction(new std::vector<double>(nBins));
	
	for (size_t i = 0; i < nPositions; ++i) {
		xPositions[i] = positions->getXPosition(i);
		yPositions[i] = positions->getYPosition(i);
	}
	
	// run the actual calculation
	VR_sp_pp2(xPositions.get(), yPositions.get(), nPositions, &nBins,
			  &((*lFunction)[0]), calculationRange, upperX, lowerX,
			  upperY, lowerY);
    
    // if the requested calculation range is too high then VR_sp_pp2 will
    // have automatically reduced it. It reports the correct number of points
    // in nBins
    if (nBins != lFunction->size()) {
        std::shared_ptr<std::vector<double> > resizedLFunction(new std::vector<double>(nBins));
        for (size_t i = 0; i < nBins; i+=1) {
            (*resizedLFunction)[i] = (*lFunction)[i];
        }
        lFunction = resizedLFunction;
    }
	
    // return the calculated function
	return lFunction;
	
}

void VR_sp_pp2(double *xCoordinates, double *yCoordinates, size_t nPoints, size_t *nBins,
			   double *outputArray, double calculationRange, double upperX, double lowerX,
			   double upperY, double lowerY) {
    size_t nIncludedBins;
    int   ib;
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
    for (size_t i = 0; i < *nBins; i++)
		outputArray[i] = 0.0;
	
    for (size_t i = 1; i < nPoints; i++) {
		xi = xCoordinates[i];
		yi = yCoordinates[i];
		for (size_t j = 0; j < i; j++) {
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
    for (size_t i = 0; i < nIncludedBins; i++) {
		sqDistance += outputArray[i];
		outputArray[i] = sqrt(sqDistance / M_PI) * sqrtArea;
    }
}

double VR_edge(double x, double y, double a, double xu0, double xl0,
			double yu0, double yl0) {
    double b, c, c1, c2, r[6], w;
    int   i;
	
	// set w to the distance from the point to
	// the closest edge
    w = x - xl0;
    if (w > y - yl0) w = y - yl0;
    if (w > xu0 - x) w = xu0 - x;
    if (w > yu0 - y) w = yu0 - y;
	
	// if the distance between the points
	// is less than the distance to the closest edge
	// then the entire circle is within the sample region
    if (a <= w) return (0.5);
	
	
    r[4] = r[0] = x - xl0;
    r[5] = r[1] = yu0 - y;
    r[2] = xu0 - x;
    r[3] = y - yl0;
    b = 0.0;
    for (i = 1; i <= 4; i++)
		if (r[i] < a) {	// the distance from this point to the edge
						// is closer than the radius of the circle
						// so some part of it is outside the region
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
