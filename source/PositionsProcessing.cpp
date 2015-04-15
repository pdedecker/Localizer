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

#include "PositionsProcessing.h"

#include "tbb/tbb.h"

/**
 * Function that corrects for edge effects in the calculated L function, copied with small modifications from Ripley et al
 */
double VR_edge(double x, double y, double pointDistance, double xu0, double xl0,
               double yu0, double yl0);

std::vector<double> CalculateClustering(bool useKFunction, std::shared_ptr<LocalizedPositionsContainer> positions,
                                        double calculationRange, size_t nBins, double lowerX, double upperX,
                                        double lowerY, double upperY,
                                        std::shared_ptr<LocalizedPositionsContainer> positions2) {
    bool haveExplicitRegion = !((lowerX == 0) && (upperX == 0) && (lowerY == 0) && (lowerX == 0));
    size_t nPositions = positions->getNPositions();
    bool isBivariate = positions2.get() != NULL;
    
	// provide that positions as arrays suitable for the actual calculation
	std::unique_ptr<double[]> xPositions(new double[nPositions]);
	std::unique_ptr<double[]> yPositions(new double[nPositions]);
	std::vector<double> calculatedFunction;
	
	for (size_t i = 0; i < nPositions; ++i) {
		xPositions[i] = positions->getXPosition(i);
		yPositions[i] = positions->getYPosition(i);
	}
    
    // get the x and y boundaries for the positions if they have not been provided
    if (!haveExplicitRegion) {
        MinAndMax(xPositions.get(), nPositions, lowerX, upperX);
        MinAndMax(yPositions.get(), nPositions, lowerY, upperY);
    }
	
	// run the actual calculation
    if (!isBivariate) {
        calculatedFunction = VR_sp_pp2(xPositions.get(), yPositions.get(), xPositions.get(), yPositions.get(),
                  nPositions, nPositions, nBins, calculationRange, upperX, lowerX,
                  upperY, lowerY, useKFunction);
    } else {
        size_t nPositions2 = positions2->getNPositions();
        std::unique_ptr<double[]> xPositions2(new double[nPositions2]);
        std::unique_ptr<double[]> yPositions2(new double[nPositions2]);
        for (size_t i = 0; i < nPositions2; ++i) {
            xPositions2[i] = positions2->getXPosition(i);
            yPositions2[i] = positions2->getYPosition(i);
        }
        if (!haveExplicitRegion) {
            double lowerX2, upperX2, lowerY2, upperY2;
            MinAndMax(xPositions2.get(), nPositions2, lowerX2, upperX2);
            MinAndMax(yPositions2.get(), nPositions2, lowerY2, upperY2);
            lowerX = std::min(lowerX, lowerX2);
            lowerY = std::min(lowerY, lowerY2);
            upperX = std::max(upperX, upperX2);
            upperY = std::max(upperY, upperY2);
        }
        calculatedFunction = VR_sp_pp2(xPositions.get(), yPositions.get(), xPositions2.get(), yPositions2.get(),
                  nPositions, nPositions2, nBins, calculationRange, upperX, lowerX,
                  upperY, lowerY, useKFunction);
    }
	
    // return the calculated function
	return calculatedFunction;
	
}
std::vector<double> VR_sp_pp2(const double *xCoordinates1, const double *yCoordinates1, const double* xCoordinates2, const double* yCoordinates2,
                             size_t nPoints1, size_t nPoints2, size_t nBins,
                             double calculationRange, double upperX, double lowerX,
                             double upperY, double lowerY, bool isKFunction) {
    
    const bool isBivariate = !((xCoordinates1 == xCoordinates2) && (yCoordinates1 == yCoordinates2));
    
    double xSize = upperX - lowerX;
    double ySize = upperY - lowerY;
    double g = 2.0;
    double binWidth = calculationRange / static_cast<double>(nBins);
    double binsPerDistance = 1.0 / binWidth;
    double sqMaxDistance = calculationRange * calculationRange;
    
    std::vector<double> result;
    result = tbb::parallel_reduce(tbb::blocked_range<size_t>(0, nPoints1), std::vector<double>(nBins, 0.0),
                                  [=](const tbb::blocked_range<size_t>& r, const std::vector<double>& initialValue) -> std::vector<double> {
                                      std::vector<double> accum(initialValue);
                                      for (size_t i = r.begin(); i != r.end(); i++) {
                                          double xi = xCoordinates1[i];
                                          double yi = yCoordinates1[i];
                                          size_t upperIndex = (isBivariate) ? nPoints2 : i;
                                          for (size_t j = 0; j < upperIndex; j++) {
                                              double dx = xCoordinates2[j] - xi;
                                              double dy = yCoordinates2[j] - yi;
                                              double sqDistance = dx * dx + dy * dy;
                                              if (sqDistance < sqMaxDistance) {
                                                  double distance = sqrt(sqDistance);
                                                  size_t ib = floor(binsPerDistance * distance);
                                                  if (ib < nBins) {
                                                      if (isBivariate) {
                                                          accum[ib] += g * VR_edge(xi, yi, distance, upperX, lowerX, upperY, lowerY);
                                                      } else {
                                                          accum[ib] += g * (VR_edge(xi, yi, distance, upperX, lowerX, upperY, lowerY) + VR_edge(xCoordinates1[j], yCoordinates1[j], distance, upperX, lowerX, upperY, lowerY));
                                                      }
                                                  }
                                              }
                                          }
                                      }
                                      return accum;
                                  },
                                  [=](const std::vector<double>& result1, const std::vector<double>& result2) -> std::vector<double> {
                                      std::vector<double> accumulated(result1);
                                      for (size_t i = 0; i < accumulated.size(); i++) {
                                          accumulated[i] += result2[i];
                                      }
                                      return accumulated;
                                  });
    
    if (isKFunction) {
        double sqrtArea = sqrt(xSize * ySize);
        double accum = 0.0;
        for (size_t i = 0; i < result.size(); i++) {
            accum += result[i];
            result[i] = sqrt(accum / (M_PI * nPoints1 * nPoints2)) * sqrtArea;
        }
    } else {    // pairwise correlation
        double probeDensity2 = (double)nPoints2 / (xSize * ySize);
        for (size_t i = 0; i < result.size(); i++) {
            double r = i * binWidth;
            result[i] = result[i] / (M_PI * ((r + binWidth) * (r + binWidth) - (r * r)) * probeDensity2 * nPoints1);
        }
    }
    
    return result;
}

double VR_edge(double x, double y, double pointDistance, double xu0, double xl0,
			double yu0, double yl0) {
    double b, c, c1, c2, r[6];
	
	// set w to the distance from the point to
	// the closest edge
    double w = std::min(std::min(x - xl0, y - yl0), std::min(xu0 - x, yu0 - y));
	
	// if the distance between the points
	// is less than the distance to the closest edge
	// then the entire circle is within the sample region
    if (pointDistance <= w) return (0.5);
	
	
    r[4] = r[0] = x - xl0;
    r[5] = r[1] = yu0 - y;
    r[2] = xu0 - x;
    r[3] = y - yl0;
    b = 0.0;
    for (int i = 1; i <= 4; i++)
		if (r[i] < pointDistance) {	// the distance from this point to the edge
						// is closer than the radius of the circle
						// so some part of it is outside the region
			if (r[i] == 0.0)
				b += M_PI;
			else {
				c = acos(r[i] / pointDistance);
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
