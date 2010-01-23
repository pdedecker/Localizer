/*
 *  PALM_analysis_PALMImages.cpp
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 01/10/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "PALM_analysis_PALMImages.h"

NormalCDFLookupTable::NormalCDFLookupTable() {
	// allocate an array of 1001 doubles containing the CDF of a normal distribution with stddev 5 between -5 and 5 (inclusive)
	this->cdfTable = boost::shared_array<double> (new double[1001]);
	for (size_t i = 0; i < 1001; ++i) {
		this->cdfTable[i] = gsl_cdf_gaussian_P(-5.0 + 0.01 * (double)i, 1.0);
	}
}

double NormalCDFLookupTable::getNormalCDF(double x, double sigma) {
	assert (this->cdfTable.get() != NULL);
	
	double rescaledX, lowerBracketX;
	size_t lowerBracket, upperBracket;
	
	// rescale the requested x to a distribution with stddev 1
	rescaledX = x / sigma;
	
	if (rescaledX < -5.0)
		return 0;
	if (rescaledX > 5.0)
		return 1;
	
	lowerBracket = (size_t)((rescaledX + 5.0) / 0.01);
	upperBracket = lowerBracket + 1;
	lowerBracketX = -5.0 + lowerBracket * 0.01;
	
	
	// obtain the output value by linear interpolation
	return cdfTable[lowerBracket] + (rescaledX - lowerBracketX) / 0.01 * (cdfTable[upperBracket] - cdfTable[lowerBracket]);
}

boost::shared_ptr<PALMMatrix<float> > PALMBitmapImageCalculator::CalculateImage(boost::shared_ptr<LocalizedPositionsContainer> positions, size_t xSize, 
																				size_t ySize, size_t imageWidth, size_t imageHeight) {
	int status;
	double fittedXPos, fittedYPos, fittedIntegral;
	double centerX, centerY, calculatedIntegral, calculatedDeviation;
	size_t startX, endX, startY, endY;
	double integralX, integralY;
	double lowerXEdgeDistance, higherXEdgeDistance, lowerYEdgeDistance, higherYEdgeDistance;
	
	size_t nPositions = positions->getNPositions();
	boost::shared_ptr<PALMMatrix<float> > outputImage(new PALMMatrix<float> (imageWidth, imageHeight));
	outputImage->set_all(0);
	
	double imageWidthScaleFactor = (double)(imageWidth - 1) / (double)xSize;
	double imageHeightScaleFactor = (double)(imageHeight - 1) / (double)ySize;
	
	// update the progress reporter
	this->progressReporter->CalculationStarted();
	
	for (size_t n = 0; n < nPositions; ++n) {
		fittedIntegral = positions->getIntegral(n);
		fittedXPos = positions->getXPosition(n);
		fittedYPos = positions->getYPosition(n);
		
		if ((fittedXPos < 0) || (fittedXPos >= xSize) || (fittedYPos < 0) || (fittedYPos >= ySize)) {
			continue;
		}
		
		if (n%100 == 0) {
			// every 100 iterations provide a progress update
			this->progressReporter->UpdateCalculationProgress((double)n / (double)nPositions * 100.0);
			
			#ifdef WITH_IGOR
			status = CheckAbort(0);
			if (status == -1) {
				// abort the calculation
				// just return the image that has been calculated now
				this->progressReporter->CalculationAborted();
				return outputImage;
			}
			#endif
		}
		
		calculatedDeviation = (this->devationCalculator->getDeviation(positions, n) * imageWidthScaleFactor);
		
		// the amplitude to use when constructing the bitmap depends on the chosen weighing method
		switch (this->emitterWeighingMethod) {
			case PALMBITMAP_EMITTERWEIGHING_SAME:
				calculatedIntegral = 1;
				break;
			case PALMBITMAP_EMITTERWEIGHING_INTEGRAL:
				calculatedIntegral = fittedIntegral;
				break;
			default:
				throw (std::runtime_error("Unrecognized emitter weighing method while calculating a PALM bitmap"));
		}
		
		centerX = (size_t)(fittedXPos * imageWidthScaleFactor + 0.5);
		centerY = (size_t)(fittedYPos * imageHeightScaleFactor + 0.5);
		
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
				// take into account that each pixel really should contain the integral of the Gaussian
				// for that calculate the distance of the edge of every pixel to the emitter (in 1D)
				lowerXEdgeDistance = centerX - (double)i - 0.5;
				higherXEdgeDistance = centerX - (double)i + 0.5;
				
				lowerYEdgeDistance = centerY - (double)j - 0.5;
				higherYEdgeDistance = centerY - (double)j + 0.5;
				
				integralX = this->cdfTable.getNormalCDF(higherXEdgeDistance, calculatedDeviation) - this->cdfTable.getNormalCDF(lowerXEdgeDistance, calculatedDeviation);
				integralY = this->cdfTable.getNormalCDF(higherYEdgeDistance, calculatedDeviation) - this->cdfTable.getNormalCDF(lowerYEdgeDistance, calculatedDeviation);
				
				(*outputImage)(i, j) = (*outputImage)(i, j) + integralX * integralY * calculatedIntegral;
			}
		}
	}
	
	this->progressReporter->CalculationDone();
	
	return outputImage;
}
