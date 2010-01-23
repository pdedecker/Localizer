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
#include "PALM_analysis_storage.h"
#include <gsl/gsl_cdf.h>

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
	boost::shared_array<double> cdfTable;
};

/**
 * @brief A class responsible for calculating a bitmap PALM output image by adding an individual Gaussian for every position to the image.
 * The output image is returned as a PALMMatrix<float>.
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
	
	boost::shared_ptr<PALMMatrix<float> > CalculateImage(boost::shared_ptr<LocalizedPositionsContainer> positions, size_t xSize, 
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
 * @brief Assumes that the positional error can be estimated as the width of the point-spread function divided
 * by some scalefactor times the square root of the fitted intensity.
 */
class PALMBitmapImageDeviationCalculator_IntegralSquareRoot : public PALMBitmapImageDeviationCalculator {
public:
	PALMBitmapImageDeviationCalculator_IntegralSquareRoot(double PSFWidth_rhs, double scaleFactor_rhs) {PSFWidth = PSFWidth_rhs; scaleFactor = scaleFactor_rhs;}
	PALMBitmapImageDeviationCalculator_IntegralSquareRoot() {;}
	
	double getDeviation(boost::shared_ptr<LocalizedPositionsContainer> positions, size_t index) {return PSFWidth / (scaleFactor * sqrt(positions->getIntegral(index)));}
	
private:
	double PSFWidth;
	double scaleFactor;
};

#endif
