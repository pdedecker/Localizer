/*
 *  PALM_analysis_PALMImages.cpp
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 01/10/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "PALM_analysis_PALMImages.h"

boost::shared_ptr<PALMMatrix<float> > PALMBitmapImageCalculator::CalculateImage(boost::shared_ptr<PALMMatrix<double> > positions, size_t xSize, 
																				size_t ySize, size_t imageWidth, size_t imageHeight) {
	double fittedXPos, fittedYPos, fittedAmplitude;
	double centerX, centerY, calculatedAmplitude;
	long startX, endX, startY, endY;
	double sqrt2pi = sqrt(2 * pi);
	double distanceXSquared, distanceYSquared, currentIntensity;
	
	size_t nPositions = positions->getXSize();
	boost::shared_ptr<PALMMatrix<float> > outputImage(new PALMMatrix<float> (imageWidth, imageHeight));
	outputImage->set_all(0);
	
	double imageWidthScaleFactor = (double)(imageWidth - 1) / (double)xSize;
	double imageHeightScaleFactor = (double)(imageHeight - 1) / (double)ySize;
	double integratedFittedIntensity;
	
	for (size_t n = 0; n < nPositions; ++n) {
		fittedAmplitude = positions->get(n, 1);	// warning: magic numbers
		fittedDeviation = positions->get(n, 2);
		fittedXPos = positions->get(n, 3);
		fittedYPos = positions->get(n, 4);
		integratedFittedIntensity = sqrt2pi * fittedDeviation * fittedAmplitude;
		
		calculatedDeviation = (calculatedDeviationCalculator->getcalculatedDeviation(positions, n) * imageScaleFactor);
		calculatedAmplitude = integratedFittedIntensity / (sqrt2pi * calculatedDeviation);	// keep the integrated intensity the same
		
		if ((fittedAmplitude < 0) || (fittedX < 0) || (fittedX >= xSize) || (fittedY < 0) || (fittedY >= ySize)) {
			continue;
		}
		
		centerX = (size_t)(fittedX * imageWidthScaleFactor + 0.5);
		centerY = (size_t)(fittedY * imageHeightScaleFactor + 0.5);
		
		startX = floor((double)centerX - 4.0 * calculatedDeviation);	// only run the calculation over a subset of the image surrounding the position
		startY = floor((double)centerY - 4.0 * calculatedDeviation);
		endX = ceil((double)centerX + 4.0 * calculatedDeviation);
		endY = ceil((double)centerY + 4.0 * calculatedDeviation);
		
		if (startX < 0)
			startX = 0;
		if (endX >= imageWidth)
			endX = imageWidth - 1;
		if (startY < 0)
			startY = 0;
		if (endY >= imageHeight)
			endY = imageHeight - 1;
		
		for (size_t i = startX; i <= endX; ++i) {
			for (size_t j = startY; j < endY; ++j) {
				distanceXSquared = ((double)i - (double)centerX) * ((double)i - (double)centerX);
				distanceYSquared = ((double)j - (double)centerY) * ((double)j - (double)centerY);
				currentIntensity = calculatedAmplitude * exp(- (distanceXSquared + distanceYSquared) / (2 * calculatedDeviation * calculatedDeviation));
				
				outputImage(i, j) += currentIntensity;
			}
		}
	}
	return outputImage;
}
