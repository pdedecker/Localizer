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

#include "PALMImages.h"
#include <gsl/gsl_sf.h>
#include <tbb/tbb.h>

ERFLookupTable::ERFLookupTable() {
    this->lowerLimit = -4.0;
    this->upperLimit = 4.0;
    this->stride = 0.0001;
    
    int nValues = (upperLimit - lowerLimit) / stride + 1;
    
	this->erfTable.resize(nValues);
    
    tbb::parallel_for(tbb::blocked_range<int>(0, nValues), [=](const tbb::blocked_range<int>& r) {
        for (int i = r.begin(); i < r.end(); ++i) {
            this->erfTable[i] = gsl_sf_erf(this->lowerLimit + static_cast<double>(i) * stride);
        }
	});
}

ERFLookupTable::~ERFLookupTable() {
    
}

double ERFLookupTable::erf(double x) const {
    if (x < this->lowerLimit)
        return -1.0;
    if (x > this->upperLimit)
        return 1.0;
    
    // perform linear interpolation using a weighted average
    double fractionalIndex = (x - this->lowerLimit) / stride;
    int index = std::floor(fractionalIndex);
    double fraction = fractionalIndex - static_cast<double>(index);
    double erf = erfTable[index] * (1.0 - fraction) + erfTable[index + 1] * fraction;
    return erf;
}

ImagePtr PALMBitmapImageCalculator::CalculateImage(std::shared_ptr<LocalizedPositionsContainer> positions, size_t xSize, 
																				size_t ySize, double outputScaleFactor) {
    if (outputScaleFactor <= 0.0) {
        throw std::runtime_error("found negative or zero outputScaleFactor");
    }
    
    int imageWidth = static_cast<double>(xSize) * outputScaleFactor;
    int imageHeight = static_cast<double>(ySize) * outputScaleFactor;
	size_t nPositions = positions->getNPositions();
	ImagePtr outputImage(new Image((int)imageWidth, (int)imageHeight));
	outputImage->setConstant(0.0);
	
	double scaleFactor = static_cast<double>(imageWidth - 1) / static_cast<double>(xSize - 1);
    double scaledPixelSize = 1.0 / scaleFactor;
    double halfPixelSize = 0.5 * scaledPixelSize;
	
	// update the progress reporter
	this->progressReporter->CalculationStarted();
	
	for (size_t n = 0; n < nPositions; ++n) {
		double fittedIntegral = positions->getIntegral(n);
		double fittedXPos = positions->getXPosition(n);
		double fittedYPos = positions->getYPosition(n);

		if ((fittedIntegral == 0.0) && (this->emitterWeighingMethod == PALMBITMAP_EMITTERWEIGHING_INTEGRAL)) {
			continue;
		}
		
		if (n%100 == 0) {
			// progress update
			int progressStatus = this->progressReporter->UpdateCalculationProgress((double)n / (double)nPositions * 100.0, 100.0);
			if (progressStatus != 0) {
				// abort the calculation
				// just return the image that has been calculated now
				this->progressReporter->CalculationAborted();
				return outputImage;
			}
		}
		
		double calculatedDeviation = (this->devationCalculator->getDeviation(positions, n));
        double s = std::sqrt(2 * calculatedDeviation * calculatedDeviation);
		
		// the amplitude to use when constructing the bitmap depends on the chosen weighing method
        double calculatedIntegral = 1.0;
		switch (this->emitterWeighingMethod) {
			case PALMBITMAP_EMITTERWEIGHING_SAME:
				calculatedIntegral = 1.0;
				break;
			case PALMBITMAP_EMITTERWEIGHING_INTEGRAL:
				calculatedIntegral = fittedIntegral;
				break;
			default:
				throw (std::runtime_error("Unrecognized emitter weighing method while calculating a PALM bitmap"));
		}
		
		double centerRow = fittedXPos * scaleFactor;
		double centerCol = fittedYPos * scaleFactor;
		
		int startRow = floor((double)centerRow - 4.0 * calculatedDeviation * scaleFactor);	// only run the calculation over a subset of the image surrounding the position
		int startCol = floor((double)centerCol - 4.0 * calculatedDeviation * scaleFactor);
		int endRow = ceil((double)centerRow + 4.0 * calculatedDeviation * scaleFactor);
		int endCol = ceil((double)centerCol + 4.0 * calculatedDeviation * scaleFactor);
		
		if (startRow < 0)
			startRow = 0;
		if (endRow >= imageWidth)
			endRow = imageWidth - 1;
		if (startCol < 0)
			startCol = 0;
		if (endCol >= imageHeight)
			endCol = imageHeight - 1;
		
        for (int j = startCol; j <= endCol; ++j) {
            for (int i = startRow; i <= endRow; ++i) {
                double minX = static_cast<double>(i) * scaledPixelSize - halfPixelSize;
                double maxX = minX + scaledPixelSize;
                double minY = static_cast<double>(j) * scaledPixelSize - halfPixelSize;
                double maxY = minY + scaledPixelSize;
                double x0 = fittedXPos;
                double y0 = fittedYPos;
                
                double erfMinX = erfTable.erf((minX - x0) / s);
                double erfMaxX = erfTable.erf((maxX - x0) / s);
                double erfMinY = erfTable.erf((minY - y0) / s);
                double erfMaxY = erfTable.erf((maxY - y0) / s);
                double integral = 0.25 * (erfMinX * erfMinY - erfMaxX * erfMinY - erfMinX * erfMaxY + erfMaxX * erfMaxY);
                
                (*outputImage)(i, j) += integral * calculatedIntegral;
            }
        }
	}
	
	this->progressReporter->CalculationDone();
	
	return outputImage;
}

double PALMBitmapImageDeviationCalculator_GaussianMask::getDeviation(std::shared_ptr<LocalizedPositionsContainer> positions, size_t index) {
	double integral = positions->getIntegral(index);
	double background = positions->getBackground(index);
	
	if ((integral == -1.0) || (background == -1.0))
		throw std::runtime_error("the used positions do not provide a PSF width and/or background estimate (use a different localization algorithm)");
		
	background = floor((background - cameraOffset) / cameraMultiplicationFactor + 0.5);
	background = (background < 0.0) ? 0.0 : background;
	integral = floor(integral / cameraMultiplicationFactor + 0.5);
	integral = (integral < 0.0) ? 0.0 : integral;
	
	double deviation;
	double psfWidthSquared = PSFWidth * PSFWidth;
	deviation = (psfWidthSquared + 1.0 / 12.0) / integral;
	//deviation += 8 * M_PI * psfWidthSquared * psfWidthSquared * background * background / (integral * integral);
	return sqrt(deviation);
}
