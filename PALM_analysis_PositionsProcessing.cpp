/*
 *  PALM_analysis_PositionsProcessing.cpp
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 06/01/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "PALM_analysis_PositionsProcessing.h"


boost::shared_ptr<std::vector<double> > CalculateKFunctionClustering(boost::shared_ptr<LocalizedPositionsContainer> positions, double binWidth, size_t nBins) {
	assert(binWidth > 0);
	assert(nBins > 0);
	size_t nPositions = positions->getNPositions();
	double currentX, currentY;
	
	double minX = 1e200, maxX = -1, minY = 1e200, maxY = -1;
	// get the x and y boundaries for the positions
	for (size_t i = 0; i < nPositions; ++i) {
		currentX = positions->getXPosition(i);
		currentY = positions->getYPosition(i);
		
		if (currentX < minX)
			minX = currentX;
		if (currentX > maxX)
			maxX = currentX;
		if (currentY < minY)
			minY = currentY;
		if (currentY > maxY)
			maxY = currentY;
	}
	
	// not all points are included since the points on the edge
	// do not have a full set of neighbours
	double requiredMargin = (nBins + 1) * binWidth;
	double lowerXLimit = minX + requiredMargin;
	double upperXLimit = maxX - requiredMargin;
	double lowerYLimit = minY + requiredMargin;
	double upperYLimit = maxY - requiredMargin;
	
	double effectiveWidth = maxX - minX;
	double effectiveHeight = maxY - minY;
	
	// the main loop
	double distance;
	double maxDistance = nBins * binWidth;
	size_t startBin;
	boost::shared_ptr<std::vector<double> > kFunction(new std::vector<double>());
	kFunction->resize(nBins);
	for (size_t i = 0; i < nBins; ++i) {
		(*kFunction)[i] = 0;
	}
	
	for (size_t i = 0; i < nPositions; ++i) {
		for (size_t j = i + 1; j < nPositions; ++j) {
			
			distance = sqrt((positions->getXPosition(i) - positions->getXPosition(j)) * (positions->getXPosition(i) - positions->getXPosition(j))
							+ (positions->getYPosition(i) - positions->getYPosition(j)) * (positions->getYPosition(i) - positions->getYPosition(j)));
			
			if (distance > maxDistance)
				continue;
			
			startBin = distance / binWidth;
			for (size_t k = startBin; k < nBins; ++k) {
				(*kFunction)[k] = (*kFunction)[k] + 2;
			}
		}
	}
	
	// introduce a correction for edge effects as found in http://e.marcon.free.fr/download/GeneralizingRipleysKFunctionToInhomogeneousPopulations.pdf
	// equation 15
	double lowerBinLimit, upperBinLimit, binCenter;
	double correctionFactor;
	for (size_t i = 0; i < nBins; ++i) {
		lowerBinLimit = i * binWidth;
		upperBinLimit = (i + 1) * binWidth;
		binCenter = (upperBinLimit - lowerBinLimit) / 2.0;
		
		correctionFactor = 1 - 4.0 / (3.0 * PI) * (binCenter / (maxX - minX) + binCenter / (maxY - minY));
		correctionFactor += (11.0 / 3.0 / PI - 1) * (binCenter * binCenter / (maxX - minX) / (maxY - minY));
		
		(*kFunction)[i] = (*kFunction)[i] / correctionFactor;
	}
		
	// apply the required normalizations
	double effectiveArea = effectiveWidth * effectiveHeight;
	for (size_t i = 0; i < nBins; ++i) {
		(*kFunction)[i] = (*kFunction)[i] * effectiveArea / (nPositions * nPositions);
	}
	
	return kFunction;
	
}
