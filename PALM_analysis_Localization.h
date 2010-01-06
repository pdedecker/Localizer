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

#define GSL_RANGE_CHECK_OFF	// this is not required since PALMMatrix<double> does range checks

#include "PALM_analysis.h"
#include "PALM_analysis_storage.h"
#include "boost/smart_ptr.hpp"
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_blas.h>
#include <stdexcept>

using namespace std;
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
	
	boost::shared_ptr<LocalizedPositionsContainer> fit_positions(const boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<std::vector<position> > positions);
	virtual boost::shared_ptr<LocalizedPositionsContainer> fit_positions(const boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<std::vector<position> >, size_t startPos, size_t endPos) = 0;
	// the second function fits the positions between startPos and endPos (indices in the array passed in positions)
	// it's mainly provided to help with multithreading
};

/**
 * @brief Localizes the particles using nonlinear least-squares Levenberg-Marquardt fitting of a symmetric 2D Gaussian.
 */
class FitPositionsGaussian : public FitPositions {
public:
	FitPositionsGaussian(size_t cutoff_radius_rhs, double r_initial_rhs, double sigma_rhs) {cutoff_radius = cutoff_radius_rhs; r_initial = r_initial_rhs; sigma = sigma_rhs;}
	~FitPositionsGaussian() {;}
	
	boost::shared_ptr<LocalizedPositionsContainer> fit_positions(const boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<std::vector<position> > positions, size_t startPos, size_t endPos);
	
protected:
	double sigma;
	size_t cutoff_radius;
	double r_initial;
};

/**
 * @brief Localizes the particles using nonlinear least-squares Levenberg-Marquardt fitting of a symmetric 2D Gaussian with a fixed standard deviation
 */
class FitPositionsGaussian_FixedWidth : public FitPositions {
public:
	FitPositionsGaussian_FixedWidth(size_t cutoff_radius_rhs, double r_initial_rhs, double sigma_rhs) {cutoff_radius = cutoff_radius_rhs; r_initial = r_initial_rhs; sigma = sigma_rhs;}
	~FitPositionsGaussian_FixedWidth() {;}
	
	boost::shared_ptr<LocalizedPositionsContainer> fit_positions(const boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<std::vector<position> > positions, size_t startPos, size_t endPos);
	
protected:
	double sigma;
	size_t cutoff_radius;
	double r_initial;
};

/*class FitPositionsEllipsoidalGaussian : public FitPositions {
public:
	FitPositionsEllipsoidalGaussian(size_t cutoff_radius_rhs, double r_initial_rhs, double sigma_rhs) {cutoff_radius = cutoff_radius_rhs; r_initial = r_initial_rhs; sigma = sigma_rhs;}
	~FitPositionsEllipsoidalGaussian() {;}
	
	boost::shared_ptr<std::vector<LocalizedPosition> > fit_positions(const boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<PALMMatrix<double> > positions, size_t startPos, size_t endPos);
	
protected:
	double sigma;
	size_t cutoff_radius;
	double r_initial;
};*/

/**
 * @brief fits the positions by doing an interative multiplication of the data with a Gaussian at the current best-guess position
 * if this converges then we assume that we have found the actual position. A description is given in Thompson, Biophys J 82:2775 2002.
 */
class FitPositionsMultiplication : public FitPositions {
public:
	FitPositionsMultiplication(size_t cutoff_radius_rhs, double r_initial_rhs, double convergence_rhs) {cutoff_radius = cutoff_radius_rhs; r_initial = r_initial_rhs; convergence_threshold = convergence_rhs;}
	~FitPositionsMultiplication() {;}
	
	// r_initial should be the standard deviation of the Gaussian
	boost::shared_ptr<LocalizedPositionsContainer> fit_positions(const boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<std::vector<position> > positions, size_t startPos, size_t endPos);
	
protected:
	int multiply_with_gaussian(boost::shared_ptr<PALMMatrix<double> > original_image, boost::shared_ptr<PALMMatrix<double> > masked_image, double x, double y,
							   double std_dev, double background, double amplitude);
	// masked_image should be provided with the same dimensions as original_image. It will be overwritten with the contents of the multiplication
	int determine_x_y_position(boost::shared_ptr<PALMMatrix<double> > masked_image, double &x, double &y);
	
	double convergence_threshold;
	size_t cutoff_radius;
	double r_initial;
};

/**
 * @brief fits the positions by calculating weighted averages for x and y using the intensity of the emission at each point as the weight.
 */
class FitPositionsCentroid : public FitPositions {
	// fits the positions by calculating a centroid for the pixel values
	
public:
	FitPositionsCentroid(size_t cutoff_radius_rhs) {cutoff_radius = cutoff_radius_rhs;}
	~FitPositionsCentroid() {;}
	
	// r_initial should be the standard deviation of the Gaussian
	boost::shared_ptr<LocalizedPositionsContainer> fit_positions(const boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<std::vector<position> > positions, size_t startPos, size_t endPos);
	
protected:
	size_t cutoff_radius;
};





// the routines below are used in the least-squares fitting of a Gaussian to the spots

int Gauss_2D_fit_function(const gsl_vector *params, void *measured_intensities, gsl_vector *deviations);
int Gauss_2D_fit_function_FixedWidth(const gsl_vector *params, void *measured_intensities, gsl_vector *deviations);
int EllipsoidalGauss_2D_fit_function(const gsl_vector *params, void *measured_intensities, gsl_vector *deviations);


int Gauss_2D_fit_function_Jacobian(const gsl_vector *params, void *measured_intensities, gsl_matrix *jacobian);
int Gauss_2D_fit_function_FixedWidth(const gsl_vector *params, void *measured_intensities, gsl_matrix *jacobian);
int Elliposoidal_Gauss_2D_fit_function_Jacobian(const gsl_vector *params, void *measured_intensities, gsl_matrix *jacobian);

int Gauss_2D_fit_function_and_Jacobian(const gsl_vector *params, void *measured_intensities_struct, gsl_vector *model_values, gsl_matrix *jacobian);
int Gauss_2D_fit_function_and_Jacobian_FixedWidth(const gsl_vector *params, void *measured_intensities_struct, gsl_vector *model_values, gsl_matrix *jacobian);
int EllipsoidalGauss_2D_fit_function_and_Jacobian(const gsl_vector *params, void *measured_intensities_struct, gsl_vector *model_values, gsl_matrix *jacobian);

// the routines below are simply adapted versions of the least-squares routines above, but have been 'tweaked' to approximate Poissonian instead of Gaussian error distributions

// int Gauss_2D_Poissonian_fit_function(const gsl_vector *params, void *measured_intensities, gsl_vector *model_values);

// int Gauss_2D_Poissonian_fit_function_Jacobian(const gsl_vector *params, void *measured_intensities, gsl_matrix *jacobian);

// int Gauss_2D_Poissonian_fit_function_and_Jacobian(const gsl_vector *params, void *measured_intensities_struct, gsl_vector *model_values, gsl_matrix *jacobian);

class measured_data_Gauss_fits {
public:	
	measured_data_Gauss_fits() {;}
	~measured_data_Gauss_fits() {;}
	
	double xOffset;
	double yOffset;
	double sigma;
	double width;
	boost::shared_ptr<PALMMatrix<double> > imageSubset;
};


#endif
