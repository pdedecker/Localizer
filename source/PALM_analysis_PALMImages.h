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

#ifndef PALM_ANALYSIS_PALMIMAGES_H
#define PALM_ANALYSIS_PALMIMAGES_H

#include "PALM_analysis.h"
#include "PALM_analysis_defines.h"
#include "boost/smart_ptr.hpp"
#include "boost/thread.hpp"
#include "PALM_analysis_storage.h"

class PALMBitmapImageDeviationCalculator;
class NormalCDFLookupTable;

/**
 * A class that serves as a fast way of obtaining the cumulative distribution function of a normal distribution 
 * with standard deviation sigma and mean equal to zero
 *
 */
class NormalCDFLookupTable {
public:
	NormalCDFLookupTable();
	~NormalCDFLookupTable() {;}
	
	double getNormalCDF(double x, double sigma);
	
protected:
    double lowerLimit, upperLimit;
    double stride;
	boost::shared_array<double> cdfTable;
};

/**
 * @brief A class responsible for calculating a bitmap PALM output image by adding an individual Gaussian for every position to the image.
 * The output image is returned as a Image.
 * 
 * PALMBitmapImageCalculator takes a set of positions in a LocalizedPositionsContainer and creates an output image (2D) with sizes (xSize, ySize).
 * For every position it adds a Gaussian to the output image with an integrated contribution that is either equal to to integrated intensity
 * of the fitted emission, or equal to one, depending on the value of emitterWeighingMethod_rhs. The standard deviation of the Gaussian
 * is determined by which of the PALMBitmapImageDeviationCalculator is passed into the constructor.
 */
class PALMBitmapImageCalculator {
public:
	PALMBitmapImageCalculator(std::shared_ptr<PALMBitmapImageDeviationCalculator> devationCalculator_rhs, int emitterWeighingMethod_rhs,
							  std::shared_ptr<ProgressReporter> progressReporter_rhs) {
		devationCalculator = devationCalculator_rhs; emitterWeighingMethod = emitterWeighingMethod_rhs; progressReporter = progressReporter_rhs;
	}
	~PALMBitmapImageCalculator() {;}
	
	ImagePtr CalculateImage(std::shared_ptr<LocalizedPositionsContainer> positions, size_t xSize, 
														 size_t ySize, size_t imageWidth, size_t imageHeight);
	
	
protected:
	std::shared_ptr<PALMBitmapImageDeviationCalculator> devationCalculator;
	std::shared_ptr<ProgressReporter> progressReporter;
	NormalCDFLookupTable cdfTable;
	int emitterWeighingMethod;
};

/**
 * @brief An abstract base class from which the other PALMBitmapImageDeviationCalculators derive
 */
class PALMBitmapImageDeviationCalculator {
public:
	PALMBitmapImageDeviationCalculator() {;}
	virtual ~PALMBitmapImageDeviationCalculator() {;}
	
	virtual double getDeviation(std::shared_ptr<LocalizedPositionsContainer> positions, size_t index) = 0;
};

/**
 * @brief Assumes that the uncertainty reported from the fit provides a reliable estimate of the fitting error,
 * potentially scaled by a scalefactor (such that the standard deviation of the Gaussian is scaleFactor * reported error).
 * The upperlimit value is present since occasionally the fitting routine may enter an unstable region produce erratic values.
 * Points that have reported errors over upperLimit will be rejected.
 */
class PALMBitmapImageDeviationCalculator_FitUncertainty : public PALMBitmapImageDeviationCalculator {
public:
	PALMBitmapImageDeviationCalculator_FitUncertainty(double scaleFactor_rhs, double upperLimit_rhs) {scaleFactor = scaleFactor_rhs; upperLimit = upperLimit_rhs;}
	~PALMBitmapImageDeviationCalculator_FitUncertainty() {;}
	
	double getDeviation(std::shared_ptr<LocalizedPositionsContainer> positions, size_t index) {return ((positions->getXPositionDeviation(index) < upperLimit) ? (positions->getXPositionDeviation(index) * scaleFactor) : -1);}
	
private:
	double scaleFactor;
	double upperLimit;
};

/**
 * @brief Assumes that the positional error is more or less the same for every emitter and is given as an
 * argument to the constructor.
 */
class PALMBitmapImageDeviationCalculator_Constant : public PALMBitmapImageDeviationCalculator {
public:
	PALMBitmapImageDeviationCalculator_Constant(double deviation_rhs) {deviation = deviation_rhs;}
	PALMBitmapImageDeviationCalculator_Constant() {;}
	
	double getDeviation(std::shared_ptr<LocalizedPositionsContainer> positions, size_t index) {return deviation;}
	
private:
	double deviation;
};

/**
 * @brief Makes use of the the equation for the variance proposed by Mortensen et al in nmeth.1447. This requires
 * approximate knowledge of the detected number of photons for a given position, which we estimate
 * by taking the multiplication factor of the camera as an argument.
 */
class PALMBitmapImageDeviationCalculator_GaussianMask : public PALMBitmapImageDeviationCalculator {
public:
	PALMBitmapImageDeviationCalculator_GaussianMask(double PSFWidth_rhs, double cameraOffset_rhs, double cameraMultiplicationFactor_rhs) : 
		PSFWidth(PSFWidth_rhs), 
		cameraMultiplicationFactor(cameraMultiplicationFactor_rhs),
		cameraOffset(cameraOffset_rhs)
		{}
	
	~PALMBitmapImageDeviationCalculator_GaussianMask() {;}
	
	double getDeviation(std::shared_ptr<LocalizedPositionsContainer> positions, size_t index);
	
private:
	double PSFWidth;
	double cameraMultiplicationFactor;
	double cameraOffset;
};

#endif
