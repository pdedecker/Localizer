/*
 *  PALM_analysis_PALMImages.h
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 01/10/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef PALM_ANALYSIS_PALMIMAGES_H
#define PALM_ANALYSIS_PALMIMAGES_H

#include "PALM_analysis.h"
#include "PALM_analysis_defines.h"
#include "boost/smart_ptr.hpp"
#include "boost/thread.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include "PALM_analysis_storage.h"
#include <gsl/gsl_cdf.h>

namespace ublas = boost::numeric::ublas;

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
	
	double getNormalCDF(double x, double sigma) {	// this function is defined here to make it inline
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
	
protected:
	boost::shared_array<double> cdfTable;
};

/**
 * @brief A class responsible for calculating a bitmap PALM output image by adding an individual Gaussian for every position to the image.
 * The output image is returned as a ublas::matrix<double>.
 * 
 * PALMBitmapImageCalculator takes a set of positions in a LocalizedPositionsContainer and creates an output image (2D) with sizes (xSize, ySize).
 * For every position it adds a Gaussian to the output image with an integrated contribution that is either equal to to integrated intensity
 * of the fitted emission, or equal to one, depending on the value of emitterWeighingMethod_rhs. The standard deviation of the Gaussian
 * is determined by which of the PALMBitmapImageDeviationCalculator is passed into the constructor.
 */
class PALMBitmapImageCalculator {
public:
	PALMBitmapImageCalculator(boost::shared_ptr<PALMBitmapImageDeviationCalculator> devationCalculator_rhs, int emitterWeighingMethod_rhs,
							  boost::shared_ptr<PALMAnalysisProgressReporter> progressReporter_rhs) {
		devationCalculator = devationCalculator_rhs; emitterWeighingMethod = emitterWeighingMethod_rhs; progressReporter = progressReporter_rhs;
	}
	~PALMBitmapImageCalculator() {;}
	
	boost::shared_ptr<ublas::matrix<double> > CalculateImage(boost::shared_ptr<LocalizedPositionsContainer> positions, size_t xSize, 
														 size_t ySize, size_t imageWidth, size_t imageHeight);
	
	
protected:
	boost::shared_ptr<PALMBitmapImageDeviationCalculator> devationCalculator;
	boost::shared_ptr<PALMAnalysisProgressReporter> progressReporter;
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
	
	virtual double getDeviation(boost::shared_ptr<LocalizedPositionsContainer> positions, size_t index) = 0;
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
	
	double getDeviation(boost::shared_ptr<LocalizedPositionsContainer> positions, size_t index) {return ((positions->getXPositionDeviation(index) < upperLimit) ? (positions->getXPositionDeviation(index) * scaleFactor) : -1);}
	
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
	
	double getDeviation(boost::shared_ptr<LocalizedPositionsContainer> positions, size_t index) {return deviation;}
	
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
	PALMBitmapImageDeviationCalculator_GaussianMask(double PSFWidth_rhs, double cameraMultiplicationFactor_rhs) : PSFWidth(PSFWidth_rhs), cameraMultiplicationFactor(cameraMultiplicationFactor_rhs) {}
	~PALMBitmapImageDeviationCalculator_GaussianMask() {;}
	
	double getDeviation(boost::shared_ptr<LocalizedPositionsContainer> positions, size_t index) {
		if (positions->getIntegral(index) != 0.0) {
			return sqrt(PSFWidth * PSFWidth / (positions->getIntegral(index) / cameraMultiplicationFactor) * (16.0 / 9.0 + 8.0 * M_PI * PSFWidth * PSFWidth / (positions->getIntegral(index) / cameraMultiplicationFactor)));
		} else {
			throw std::runtime_error("the used positions do not provide a PSF width estimate (use a different localization algorithm)");
		}
	}
	
private:
	double PSFWidth;
	double cameraMultiplicationFactor;
};

#endif
