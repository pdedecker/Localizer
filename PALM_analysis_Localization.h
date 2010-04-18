/*
 *  PALM_analysis_Localization.h
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 27/09/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef PALM_ANALYSIS_LOCALIZATION_H
#define PALM_ANALYSIS_LOCALIZATION_H

#define GSL_RANGE_CHECK_OFF	// this is not required since ublas::matrix<double> does range checks

#include "PALM_analysis.h"
#include "PALM_analysis_storage.h"
#include "boost/smart_ptr.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_blas.h>
#include <stdexcept>

namespace ublas = boost::numeric::ublas;

class LocalizedPosition;
class LocalizedPositionsContainer;

/**
 * @brief An abstract base class from which all other 'FitPositions' classes must derive
 * 
 * This family of classes contains algorithms that accept a single CCD image and a list of positions as input and return a list of subpixel-localized positions.
 *
 * This class is abstract, meaning that it can not be instantiated. However, it provides the virtual function 'fit_positions()', 
 * which is where derived classes should do their work. Any necessary parameters should be passed in the constructors of derived classes.
 */
class FitPositions {
public:
	FitPositions() {;}
	virtual ~FitPositions() {;}
	
	boost::shared_ptr<LocalizedPositionsContainer> fit_positions(const boost::shared_ptr<ublas::matrix<double> > image, boost::shared_ptr<std::list<position> > positions);
	virtual boost::shared_ptr<LocalizedPositionsContainer> fit_positions(const boost::shared_ptr<ublas::matrix<double> > image, boost::shared_ptr<std::list<position> > positions, std::list<position>::iterator startIt, std::list<position>::iterator endIt) = 0;
	// the second function fits the positions between startPos and endPos (indices in the array passed in positions)
	// it's mainly provided to help with multithreading
};

/**
 * @brief Localizes the particles using nonlinear least-squares Levenberg-Marquardt fitting of a symmetric 2D Gaussian.
 */
class FitPositions_SymmetricGaussian : public FitPositions {
public:
	FitPositions_SymmetricGaussian(double initialPSFWidth_rhs, double sigma_rhs) {initialPSFWidth = initialPSFWidth_rhs; sigma = sigma_rhs; cutoff_radius = std::ceil(initialPSFWidth_rhs * 4.0);}
	~FitPositions_SymmetricGaussian() {;}
	
	boost::shared_ptr<LocalizedPositionsContainer> fit_positions(const boost::shared_ptr<ublas::matrix<double> > image, boost::shared_ptr<std::list<position> > positions, std::list<position>::iterator startIt, std::list<position>::iterator endIt);
	
protected:
	double sigma;
	size_t cutoff_radius;
	double initialPSFWidth;
};

/**
 * @brief Localizes the particles using nonlinear least-squares Levenberg-Marquardt fitting of a symmetric 2D Gaussian with a fixed standard deviation
 */
class FitPositions_FixedWidthGaussian : public FitPositions {
public:
	FitPositions_FixedWidthGaussian(double initialPSFWidth_rhs, double sigma_rhs) {initialPSFWidth = initialPSFWidth_rhs; sigma = sigma_rhs; cutoff_radius = std::ceil(initialPSFWidth_rhs * 4.0);}
	~FitPositions_FixedWidthGaussian() {;}
	
	boost::shared_ptr<LocalizedPositionsContainer> fit_positions(const boost::shared_ptr<ublas::matrix<double> > image, boost::shared_ptr<std::list<position> > positions, std::list<position>::iterator startIt, std::list<position>::iterator endIt);
	
protected:
	double sigma;
	size_t cutoff_radius;
	double initialPSFWidth;
};

class FitPositions_EllipsoidalGaussian : public FitPositions {
public:
	FitPositions_EllipsoidalGaussian(double initialPSFWidth_rhs, double sigma_rhs) {initialPSFWidth = initialPSFWidth_rhs; sigma = sigma_rhs; cutoff_radius = std::ceil(initialPSFWidth_rhs * 4.0);}
	~FitPositions_EllipsoidalGaussian() {;}
	
	boost::shared_ptr<LocalizedPositionsContainer> fit_positions(const boost::shared_ptr<ublas::matrix<double> > image, boost::shared_ptr<std::list<position> > positions, std::list<position>::iterator startIt, std::list<position>::iterator endIt);
	
protected:
	double sigma;
	size_t cutoff_radius;
	double initialPSFWidth;
};

/**
 * @brief fits the positions by doing an interative multiplication of the data with a Gaussian at the current best-guess position
 * if this converges then we assume that we have found the actual position. A description is given in Thompson, Biophys J 82:2775 2002.
 */
class FitPositionsMultiplication : public FitPositions {
public:
	FitPositionsMultiplication(double initialPSFWidth_rhs, double convergence_rhs) {initialPSFWidth = initialPSFWidth_rhs; convergence_threshold = convergence_rhs; cutoff_radius = std::ceil(initialPSFWidth_rhs * 4.0);}
	~FitPositionsMultiplication() {;}
	
	// initialPSFWidth should be the standard deviation of the Gaussian
	boost::shared_ptr<LocalizedPositionsContainer> fit_positions(const boost::shared_ptr<ublas::matrix<double> > image, boost::shared_ptr<std::list<position> > positions, std::list<position>::iterator startIt, std::list<position>::iterator endIt);
	
protected:
	int multiply_with_gaussian(boost::shared_ptr<ublas::matrix<double> > original_image, boost::shared_ptr<ublas::matrix<double> > masked_image, double x, double y,
							   double std_dev, double background, double amplitude);
	// masked_image should be provided with the same dimensions as original_image. It will be overwritten with the contents of the multiplication
	int determine_x_y_position(boost::shared_ptr<ublas::matrix<double> > masked_image, double &x, double &y);
	
	double convergence_threshold;
	size_t cutoff_radius;
	double initialPSFWidth;
};

/**
 * @brief fits the positions by calculating weighted averages for x and y using the intensity of the emission at each point as the weight.
 */
class FitPositionsCentroid : public FitPositions {
	// fits the positions by calculating a centroid for the pixel values
	
public:
	FitPositionsCentroid(double initialPSFWidth_rhs) {cutoff_radius = std::ceil(initialPSFWidth_rhs * 4.0);}
	~FitPositionsCentroid() {;}
	
	// initialPSFWidth should be the standard deviation of the Gaussian
	boost::shared_ptr<LocalizedPositionsContainer> fit_positions(const boost::shared_ptr<ublas::matrix<double> > image, boost::shared_ptr<std::list<position> > positions, std::list<position>::iterator startIt, std::list<position>::iterator endIt);
	
protected:
	size_t cutoff_radius;
};

/**
 * @brief A higher-level positions fitter that implements deflation, but needs another FitPositions to do the actual fitting work
 *
 * Deflation is the idea that fitted positions can be effectively subtracted from the image to reveal other underlying emitters that were obscured
 * at first. The validity of this approach is open to debate, however, some papers in the literature make use of this approach.
 * This class does no fitting by itself, it merely takes another FitPositions as well as all the classes required for segmentation.
 * The process is iterative: the class will keep running localize-deflate circles until no more positions are recovered.
 * The returned positions will contain all positions localized, both originally and after deflation.
 */
class FitPositionsDeflate : public FitPositions {
public:
	FitPositionsDeflate(boost::shared_ptr <ThresholdImage_Preprocessor> preprocessor_rhs, boost::shared_ptr <ThresholdImage_Postprocessor> postprocessor_rhs,
						boost::shared_ptr<ThresholdImage> thresholder_rhs, boost::shared_ptr<ParticleFinder> particleFinder_rhs, 
						boost::shared_ptr<FitPositions> positionsFitter_rhs) {
		preprocessor = preprocessor_rhs; postprocessor = postprocessor_rhs; thresholder = thresholder_rhs; positionsFitter = positionsFitter_rhs; positionsFitter = positionsFitter_rhs;}
	
	~FitPositionsDeflate() {;}
	
	boost::shared_ptr<LocalizedPositionsContainer> fit_positions(const boost::shared_ptr<ublas::matrix<double> > image, boost::shared_ptr<std::list<position> > positions, std::list<position>::iterator startIt, std::list<position>::iterator endIt);
	
protected:
	boost::shared_ptr<ublas::matrix<double> > subtractLocalizedPositions(boost::shared_ptr<ublas::matrix<double> > image, boost::shared_ptr<LocalizedPositionsContainer> positions);
	
	boost::shared_ptr <ThresholdImage_Preprocessor> preprocessor;
	boost::shared_ptr <ThresholdImage_Postprocessor> postprocessor;
	boost::shared_ptr<ThresholdImage> thresholder;
	boost::shared_ptr<ParticleFinder> particleFinder;
	boost::shared_ptr<FitPositions> positionsFitter;
};

// the routines below are used in the least-squares fitting of a Gaussian to the spots

int FitFunction_SymmetricGaussian(const gsl_vector *params, void *measured_intensities, gsl_vector *deviations);
int FitFunction_FixedWidthGaussian(const gsl_vector *params, void *measured_intensities, gsl_vector *deviations);
int FitFunction_EllipsoidalGaussian(const gsl_vector *params, void *measured_intensities, gsl_vector *deviations);


int Jacobian_SymmetricGaussian(const gsl_vector *params, void *measured_intensities, gsl_matrix *jacobian);
int Jacobian_FixedWidthGaussian(const gsl_vector *params, void *measured_intensities, gsl_matrix *jacobian);
int Jacobian_EllipsoidalGaussian(const gsl_vector *params, void *measured_intensities, gsl_matrix *jacobian);

int FitFunctionAndJacobian_SymmetricGaussian(const gsl_vector *params, void *measured_intensities_struct, gsl_vector *model_values, gsl_matrix *jacobian);
int FitFunctionAndJacobian_FixedWidthGaussian(const gsl_vector *params, void *measured_intensities_struct, gsl_vector *model_values, gsl_matrix *jacobian);
int FitFunctionAndJacobian_EllipsoidalGaussian(const gsl_vector *params, void *measured_intensities_struct, gsl_vector *model_values, gsl_matrix *jacobian);

class measured_data_Gauss_fits {
public:	
	measured_data_Gauss_fits() {;}
	~measured_data_Gauss_fits() {;}
	
	double xOffset;
	double yOffset;
	double sigma;
	double width;
	boost::shared_ptr<ublas::matrix<double> > imageSubset;
};


#endif
