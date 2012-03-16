/*
 Copyright 2008-2011 Peter Dedecker.
 
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

#include "PALM_analysis_PALMImages.h"
#include <gsl/gsl_cdf.h>

NormalCDFLookupTable::NormalCDFLookupTable() {
	// allocate an array of 1001 doubles containing the CDF of a normal distribution
    this->lowerLimit = -10.0;
    this->upperLimit = 10.0;
    this->stride = 0.01;
    
    size_t nValues = (upperLimit - lowerLimit) / stride + 1;
    
	this->cdfTable = boost::shared_array<double> (new double[nValues]);
	for (size_t i = 0; i < nValues; ++i) {
		this->cdfTable[i] = gsl_cdf_gaussian_P(this->lowerLimit + static_cast<double>(i) * stride, 1.0);
	}
}

double NormalCDFLookupTable::getNormalCDF(double x, double sigma) {
    double rescaledX;
    
    // rescale the requested x to a distribution with stddev 1
    rescaledX = x / sigma;
    
    if (rescaledX < this->lowerLimit)
        return 0.0;
    if (rescaledX > this->upperLimit)
        return 1.0;
    
    // perform linear interpolation using a weighted average
    double fractionalIndex = (rescaledX - this->lowerLimit) / stride;
    size_t lowerIndex = std::floor(fractionalIndex);
    return cdfTable[lowerIndex] * (fractionalIndex - static_cast<double>(lowerIndex)) + cdfTable[lowerIndex + 1] * (static_cast<double>(lowerIndex + 1) - fractionalIndex);
}

ImagePtr PALMBitmapImageCalculator::CalculateImage(boost::shared_ptr<LocalizedPositionsContainer> positions, size_t xSize, 
																				size_t ySize, size_t imageWidth, size_t imageHeight) {
	int progressStatus;
	double fittedXPos, fittedYPos, fittedIntegral;
	double centerX, centerY, calculatedIntegral, calculatedDeviation;
	long startX, endX, startY, endY;
	double integralX, integralY;
    double halfPixelSizeX = 0.5, halfPixelSizeY = 0.5, shiftOfThisPixelX, shiftOfThisPixelY;
	
	size_t nPositions = positions->getNPositions();
	ImagePtr outputImage(new Image((int)imageWidth, (int)imageHeight));
	outputImage->setConstant(0.0);
	
	double imageWidthScaleFactor = static_cast<double>(imageWidth) / static_cast<double>(xSize);
	double imageHeightScaleFactor = static_cast<double>(imageHeight) / static_cast<double>(ySize);
	
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
			progressStatus = this->progressReporter->UpdateCalculationProgress((double)n / (double)nPositions * 100.0, 100.0);
			if (progressStatus != 0) {
				// abort the calculation
				// just return the image that has been calculated now
				this->progressReporter->CalculationAborted();
				return outputImage;
			}
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
		
		centerX = fittedXPos * imageWidthScaleFactor;
		centerY = fittedYPos * imageHeightScaleFactor;
		
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
		
		for (size_t j = startY; j < endY; ++j) {
			for (size_t i = startX; i <= endX; ++i) {
				// take into account that each pixel really should contain the integral of the Gaussian
				// how much is this pixel shifted with respect to the center of the emitter?
                shiftOfThisPixelX = static_cast<double>(i) - centerX;
                shiftOfThisPixelY = static_cast<double>(j) - centerY;
                
                integralX = this->cdfTable.getNormalCDF(shiftOfThisPixelX + halfPixelSizeX, calculatedDeviation) - this->cdfTable.getNormalCDF(shiftOfThisPixelX - halfPixelSizeX, calculatedDeviation);
                integralY = this->cdfTable.getNormalCDF(shiftOfThisPixelY + halfPixelSizeY, calculatedDeviation) - this->cdfTable.getNormalCDF(shiftOfThisPixelY - halfPixelSizeY, calculatedDeviation);
                //integralX = gsl_cdf_gaussian_P(shiftOfThisPixelX + halfPixelSizeX, calculatedDeviation) - gsl_cdf_gaussian_P(shiftOfThisPixelX - halfPixelSizeX, calculatedDeviation);
                //integralY = gsl_cdf_gaussian_P(shiftOfThisPixelY + halfPixelSizeY, calculatedDeviation) - gsl_cdf_gaussian_P(shiftOfThisPixelY - halfPixelSizeY, calculatedDeviation);
				
				(*outputImage)(i, j) += integralX * integralY * calculatedIntegral;
			}
		}
	}
	
	this->progressReporter->CalculationDone();
	
	return outputImage;
}
