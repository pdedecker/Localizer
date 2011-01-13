/*
 *  PALM_analysis_Localization.cpp
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 27/09/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "PALM_analysis_Localization.h"

boost::shared_ptr<LocalizedPositionsContainer> FitPositions_SymmetricGaussian::fit_positions(const boost::shared_ptr<Eigen::MatrixXd> image,
																							 boost::shared_ptr<std::list<position> > positions) {
																			
	// some safety checks
	if (positions->size() == 0) {
		// if no positions were found then there is no reason to run the analysis
		// we need to catch this here and not upstream since then we can return an appropriate
		// instance of LocalizedPositionsContainer
		return boost::shared_ptr<LocalizedPositionsContainer_2DGauss> (new LocalizedPositionsContainer_2DGauss());
	}
	
	size_t size_of_subset = 2 * cutoff_radius + 1;
	double x_offset, y_offset, x_max, y_max;
	size_t number_of_intensities = size_of_subset * size_of_subset;
	size_t xSize = image->size1();
	size_t ySize = image->size2();
	
	double x0_initial, y0_initial, amplitude, background;
	double chi, degreesOfFreedom, c;
	long iterations = 0;
	int status;
	
	double relativeAmplitudeError, relativeWidthError;
	
	boost::shared_ptr<Eigen::MatrixXd> image_subset;
	boost::shared_ptr<LocalizedPositionsContainer_2DGauss> fitted_positions (new LocalizedPositionsContainer_2DGauss());
	boost::shared_ptr<LocalizedPosition_2DGauss> localizationResult (new LocalizedPosition_2DGauss());
	
	image_subset = boost::shared_ptr<Eigen::MatrixXd> (new ublas::matrix<double>(size_of_subset, size_of_subset));
	
	// initialize the solver
	const gsl_multifit_fdfsolver_type *solver;
	gsl_multifit_fdfsolver *fit_iterator;
	
	solver = gsl_multifit_fdfsolver_lmsder;
	measured_data_Gauss_fits fitData;
	gsl_multifit_function_fdf f;
	
	gsl_vector *fit_parameters = gsl_vector_alloc(5);
	if (fit_parameters == NULL) {
		throw std::bad_alloc();
	}
	
	fit_iterator = gsl_multifit_fdfsolver_alloc(solver, number_of_intensities, 5);
	if (fit_iterator == NULL) {
		gsl_vector_free(fit_parameters);
		throw std::bad_alloc();
	}
	
	gsl_matrix *covarianceMatrix = gsl_matrix_alloc(5, 5);
	if (covarianceMatrix == NULL) {
		gsl_vector_free(fit_parameters);
		gsl_multifit_fdfsolver_free(fit_iterator);
		throw std::bad_alloc();
	}
	
	f.f = &FitFunction_SymmetricGaussian;
	f.df = &Jacobian_SymmetricGaussian;
	f.fdf = &FitFunctionAndJacobian_SymmetricGaussian;
	f.n = number_of_intensities;
	f.p = 5;
	f.params = (void *)&fitData;
	
	
	// iterate over all the determined positions
	std::list<position>::iterator it = positions->begin();
	while (it != positions->end()) {
		iterations = 0;
		
		amplitude = (*it).get_intensity();
		x0_initial = (*it).get_x();
		y0_initial = (*it).get_y();
		background = (*it).get_background();
		
		x_offset = floor(x0_initial - (double)cutoff_radius);
		y_offset = floor(y0_initial - (double)cutoff_radius);
		x_max = floor(x0_initial + (double)cutoff_radius);
		y_max = floor(y0_initial + (double)cutoff_radius);
		
		if ((x_offset < 0) || (x_max > (xSize - 1)) || (y_offset < 0) || (y_max > (ySize - 1))) {	// this position is too close to the edge of the image
			// we cannot include it
			it = positions->erase(it);
			continue;
		}
		
		for (size_t j = x_offset; j <= x_max; j++) {
			for (size_t k = y_offset; k <= y_max; k++) {
				(*image_subset)(j - x_offset, k - y_offset) = (*image)(j, k);
			}
		}
		
		fitData.xOffset = x_offset;
		fitData.yOffset = y_offset;
		fitData.imageSubset = image_subset;
		fitData.sigma = sigma;
		
		// provide the initial parameters
		gsl_vector_set(fit_parameters, 0, amplitude);
		gsl_vector_set(fit_parameters, 1, initialPSFWidth);
		gsl_vector_set(fit_parameters, 2, x0_initial);
		gsl_vector_set(fit_parameters, 3, y0_initial);
		gsl_vector_set(fit_parameters, 4, background);
		
		// set the solver
		gsl_multifit_fdfsolver_set(fit_iterator, &f, fit_parameters);
		
		// run the iterations
		do {
			iterations++;
			status = gsl_multifit_fdfsolver_iterate(fit_iterator);
			if (status != 0)
				break;
			
			status = gsl_multifit_test_delta(fit_iterator->dx, fit_iterator->x, 1e-4, 1e-4);
		} while ((status == GSL_CONTINUE) && (iterations < 200));
		
		if (gsl_vector_get(fit_iterator->x, 0) <= 0) {	// reject fits that have negative amplitudes
			it = positions->erase(it);
			continue;
		}
		
		// are the reported positions within the window?
		if ((gsl_vector_get(fit_iterator->x, 2) < x_offset) || (gsl_vector_get(fit_iterator->x, 2) > x_max) || (gsl_vector_get(fit_iterator->x, 3) < y_offset) || (gsl_vector_get(fit_iterator->x, 3) > y_max)) {
			// the reported positions are not within the window, we should reject them
			it = positions->erase(it);
			continue;
		}
		
		// are the fitted coordinates close enough to the initial guess?
		if ((fabs(gsl_vector_get(fit_iterator->x, 2) - x0_initial) > 2.0 * initialPSFWidth) || (fabs(gsl_vector_get(fit_iterator->x, 3) - y0_initial) > 2.0 * initialPSFWidth)) {
			it = positions->erase(it);
			continue;
		}
		
		// is the amplitude close enough to the initial value to be trusted?
		if ((gsl_vector_get(fit_iterator->x, 0) < amplitude / 2.0) || (gsl_vector_get(fit_iterator->x, 0) > amplitude * 1.5)) {
			// the output fit amplitude is more than a factor of two different from the initial value, drop this point
			it = positions->erase(it);
			continue;
		}
		
		if ((gsl_vector_get(fit_iterator->x, 1) < initialPSFWidth / 2.0) || (gsl_vector_get(fit_iterator->x, 1) > initialPSFWidth * 1.5)) {
			// the output fit width is more than a factor of two different from the initial value, drop this point
			it = positions->erase(it);
			continue;
		}
		
		// calculate the covariance matrix
		gsl_multifit_covar(fit_iterator->J, 0.0, covarianceMatrix);
		chi = gsl_blas_dnrm2(fit_iterator->f);
		degreesOfFreedom = (2 * cutoff_radius - 1) * (2 * cutoff_radius - 1) - 5;
		c = GSL_MAX_DBL(1, chi / sqrt(degreesOfFreedom));
		
		// store the fitted parameters
		
		localizationResult->width = gsl_vector_get(fit_iterator->x, 1);
		localizationResult->xPosition = gsl_vector_get(fit_iterator->x, 2);
		localizationResult->yPosition = gsl_vector_get(fit_iterator->x, 3);
		localizationResult->background = gsl_vector_get(fit_iterator->x, 4);
		
		localizationResult->widthDeviation = c * sqrt(gsl_matrix_get(covarianceMatrix, 1, 1));
		localizationResult->xPositionDeviation = c * sqrt(gsl_matrix_get(covarianceMatrix, 2, 2));
		localizationResult->yPositionDeviation = c * sqrt(gsl_matrix_get(covarianceMatrix, 3, 3));
		localizationResult->backgroundDeviation = c * sqrt(gsl_matrix_get(covarianceMatrix, 4, 4));
		
		localizationResult->integral = 2 * PI * localizationResult->width * localizationResult->width * gsl_vector_get(fit_iterator->x, 0);
		
		// use the rules for error propagation to calculate the error on the integrated intensity
		// \sigma_I	= 2 \pi I \sqrt{\left(\frac{\sigma_A}{A} \right)^2 + 2 \left(\frac{\sigma_r}{r} \right)^2}
		relativeAmplitudeError = c * sqrt(gsl_matrix_get(covarianceMatrix, 0, 0)) / gsl_vector_get(fit_iterator->x, 0); // the relative error on the amplitude
		relativeWidthError = localizationResult->widthDeviation / localizationResult->width;
		localizationResult->integralDeviation = 2 * PI * localizationResult->integral * sqrt(relativeAmplitudeError * relativeAmplitudeError + 2 * relativeWidthError * relativeWidthError);
		
		fitted_positions->addPosition(localizationResult);
		++it;
	}
	
	gsl_multifit_fdfsolver_free(fit_iterator);
	gsl_vector_free(fit_parameters);
	gsl_matrix_free(covarianceMatrix);
	
	return fitted_positions;
	
}

boost::shared_ptr<LocalizedPositionsContainer> FitPositions_FixedWidthGaussian::fit_positions(const boost::shared_ptr<Eigen::MatrixXd> image,
																							  boost::shared_ptr<std::list<position> > positions) {
	
	// some safety checks
	if (positions->size() == 0) {
		// if no positions were found then there is no reason to run the analysis
		// we need to catch this here and not upstream since then we can return an appropriate
		// instance of LocalizedPositionsContainer
		return boost::shared_ptr<LocalizedPositionsContainer_2DGaussFixedWidth> (new LocalizedPositionsContainer_2DGaussFixedWidth());
	}
	
	size_t size_of_subset = 2 * cutoff_radius + 1;
	double x_offset, y_offset, x_max, y_max;
	size_t number_of_intensities = size_of_subset * size_of_subset;
	size_t xSize = image->size1();
	size_t ySize = image->size2();
	
	double x0_initial, y0_initial, amplitude, background;
	double chi, degreesOfFreedom, c;
	long iterations = 0;
	int status;
	
	double relativeAmplitudeError;
	
	boost::shared_ptr<Eigen::MatrixXd> image_subset;
	boost::shared_ptr<LocalizedPositionsContainer_2DGaussFixedWidth> fitted_positions (new LocalizedPositionsContainer_2DGaussFixedWidth());
	boost::shared_ptr<LocalizedPosition_2DGaussFixedWidth> localizationResult (new LocalizedPosition_2DGaussFixedWidth());
	
	image_subset = boost::shared_ptr<Eigen::MatrixXd> (new ublas::matrix<double>(size_of_subset, size_of_subset));
	
	// initialize the solver
	const gsl_multifit_fdfsolver_type *solver;
	gsl_multifit_fdfsolver *fit_iterator;
	
	solver = gsl_multifit_fdfsolver_lmsder;
	measured_data_Gauss_fits fitData;
	gsl_multifit_function_fdf f;
	
	gsl_vector *fit_parameters = gsl_vector_alloc(4);
	if (fit_parameters == NULL) {
		throw std::bad_alloc();
	}
	
	fit_iterator = gsl_multifit_fdfsolver_alloc(solver, number_of_intensities, 4);
	if (fit_iterator == NULL) {
		gsl_vector_free(fit_parameters);
		throw std::bad_alloc();
	}
	
	gsl_matrix *covarianceMatrix = gsl_matrix_alloc(4, 4);
	if (covarianceMatrix == NULL) {
		gsl_vector_free(fit_parameters);
		gsl_multifit_fdfsolver_free(fit_iterator);
		throw std::bad_alloc();
	}
	
	f.f = &FitFunction_FixedWidthGaussian;
	f.df = &Jacobian_FixedWidthGaussian;
	f.fdf = &FitFunctionAndJacobian_FixedWidthGaussian;
	f.n = number_of_intensities;
	f.p = 4;
	f.params = (void *)&fitData;
	
	
	// iterate over all the determined positions
	std::list<position>::iterator it = positions->begin();
	while (it != positions->end()) {
		iterations = 0;
		
		amplitude = (*it).get_intensity();
		x0_initial = (*it).get_x();
		y0_initial = (*it).get_y();
		background = (*it).get_background();
		
		x_offset = floor(x0_initial - (double)cutoff_radius);
		y_offset = floor(y0_initial - (double)cutoff_radius);
		x_max = floor(x0_initial + (double)cutoff_radius);
		y_max = floor(y0_initial + (double)cutoff_radius);
		
		if ((x_offset < 0) || (x_max > (xSize - 1)) || (y_offset < 0) || (y_max > (ySize - 1))) {	// this position is too close to the edge of the image
			// we cannot include it
			it = positions->erase(it);
			continue;
		}
		
		for (size_t j = x_offset; j <= x_max; j++) {
			for (size_t k = y_offset; k <= y_max; k++) {
				(*image_subset)(j - x_offset, k - y_offset) = (*image)(j, k);
			}
		}
		
		fitData.xOffset = x_offset;
		fitData.yOffset = y_offset;
		fitData.imageSubset = image_subset;
		fitData.sigma = sigma;
		fitData.width = initialPSFWidth;
		
		// provide the initial parameters
		gsl_vector_set(fit_parameters, 0, amplitude);
		gsl_vector_set(fit_parameters, 1, x0_initial);
		gsl_vector_set(fit_parameters, 2, y0_initial);
		gsl_vector_set(fit_parameters, 3, background);
		
		// set the solver
		gsl_multifit_fdfsolver_set(fit_iterator, &f, fit_parameters);
		
		// run the iterations
		do {
			iterations++;
			status = gsl_multifit_fdfsolver_iterate(fit_iterator);
			if (status != 0)
				break;
			
			status = gsl_multifit_test_delta(fit_iterator->dx, fit_iterator->x, 1e-4, 1e-4);
		} while ((status == GSL_CONTINUE) && (iterations < 200));
		
		if (gsl_vector_get(fit_iterator->x, 0) <= 0) {	// reject fits that have negative amplitudes
			it = positions->erase(it);
			continue;
		}
		
		// are the reported positions within the window?
		if ((gsl_vector_get(fit_iterator->x, 1) < x_offset) || (gsl_vector_get(fit_iterator->x, 1) > x_max) || (gsl_vector_get(fit_iterator->x, 2) < y_offset) || (gsl_vector_get(fit_iterator->x, 2) > y_max)) {
			// the reported positions are not within the window, we should reject them
			it = positions->erase(it);
			continue;
		}
		
		// are the fitted coordinates close enough to the initial guess?
		if ((fabs(gsl_vector_get(fit_iterator->x, 1) - x0_initial) > 2.0 * initialPSFWidth) || (fabs(gsl_vector_get(fit_iterator->x, 2) - y0_initial) > 2.0 * initialPSFWidth)) {
			it = positions->erase(it);
			continue;
		}
		
		// is the fitted amplitude close enough to the initial value to be trusted?
		if ((gsl_vector_get(fit_iterator->x, 0) < amplitude / 2.0) || (gsl_vector_get(fit_iterator->x, 0) > amplitude * 2.0)) {
			// the output fit amplitude is more than a factor of two different from the initial value, drop this point
			it = positions->erase(it);
			continue;
		}
		
		// calculate the covariance matrix
		gsl_multifit_covar(fit_iterator->J, 0.0, covarianceMatrix);
		chi = gsl_blas_dnrm2(fit_iterator->f);
		degreesOfFreedom = (2 * cutoff_radius - 1) * (2 * cutoff_radius - 1) - 4;
		c = GSL_MAX_DBL(1, chi / sqrt(degreesOfFreedom));
		
		// store the data
		// store the fitted parameters
		localizationResult->xPosition = gsl_vector_get(fit_iterator->x, 1);
		localizationResult->yPosition = gsl_vector_get(fit_iterator->x, 2);
		localizationResult->background = gsl_vector_get(fit_iterator->x, 3);
		
		localizationResult->integral = 2 * PI * initialPSFWidth * initialPSFWidth * gsl_vector_get(fit_iterator->x, 0);
		
		
		localizationResult->xPositionDeviation = c * sqrt(gsl_matrix_get(covarianceMatrix, 1, 1));
		localizationResult->yPositionDeviation = c * sqrt(gsl_matrix_get(covarianceMatrix, 2, 2));
		localizationResult->backgroundDeviation = c * sqrt(gsl_matrix_get(covarianceMatrix, 3, 3));
		
		// calculate the uncertainty on the integrated intensity using error propagation
		relativeAmplitudeError = c * sqrt(gsl_matrix_get(covarianceMatrix, 0, 0)) / gsl_vector_get(fit_iterator->x, 0); // the relative error on the amplitude
		localizationResult->integralDeviation = 2 * PI * localizationResult->integral * initialPSFWidth * initialPSFWidth * sqrt(relativeAmplitudeError * relativeAmplitudeError);
		
		
		fitted_positions->addPosition(localizationResult);
		++it;
	}
	
	gsl_multifit_fdfsolver_free(fit_iterator);
	gsl_vector_free(fit_parameters);
	gsl_matrix_free(covarianceMatrix);
	
	return fitted_positions;
	
}


boost::shared_ptr<LocalizedPositionsContainer> FitPositions_EllipsoidalGaussian::fit_positions(const boost::shared_ptr<Eigen::MatrixXd> image,
																							   boost::shared_ptr<std::list<position> > positions) {
	
	// some safety checks
	if (positions->size() == 0) {
		// if no positions were found then there is no reason to run the analysis
		// we need to catch this here and not upstream since then we can return an appropriate
		// instance of LocalizedPositionsContainer
		return boost::shared_ptr<LocalizedPositionsContainer> (new LocalizedPositionsContainer_Ellipsoidal2DGaussian());
	}
	
	size_t size_of_subset = 2 * cutoff_radius + 1;
	double x_offset, y_offset, x_max, y_max;
	size_t number_of_intensities = size_of_subset * size_of_subset;
	size_t xSize = image->size1();
	size_t ySize = image->size2();
	
	double x0_initial, y0_initial, amplitude, background, correlation_initial;
	double chi, degreesOfFreedom, c;
	long iterations = 0;
	int status;
	
	boost::shared_ptr<Eigen::MatrixXd> image_subset;
	boost::shared_ptr<LocalizedPositionsContainer_Ellipsoidal2DGaussian> fitted_positions (new LocalizedPositionsContainer_Ellipsoidal2DGaussian());
	boost::shared_ptr<LocalizedPosition_Ellipsoidal2DGauss> localizationResult (new LocalizedPosition_Ellipsoidal2DGauss());
	
	image_subset = boost::shared_ptr<Eigen::MatrixXd> (new ublas::matrix<double>(size_of_subset, size_of_subset));
	
	// initialize the solver
	const gsl_multifit_fdfsolver_type *solver;
	gsl_multifit_fdfsolver *fit_iterator;
	
	solver = gsl_multifit_fdfsolver_lmsder;
	measured_data_Gauss_fits fitData;
	gsl_multifit_function_fdf f;
	
	gsl_vector *fit_parameters = gsl_vector_alloc(7);
	if (fit_parameters == NULL) {
		throw std::bad_alloc();
	}
	
	fit_iterator = gsl_multifit_fdfsolver_alloc(solver, number_of_intensities, 7);
	if (fit_iterator == NULL) {
		gsl_vector_free(fit_parameters);
		throw std::bad_alloc();
	}
	
	gsl_matrix *covarianceMatrix = gsl_matrix_alloc(7, 7);
	if (covarianceMatrix == NULL) {
		gsl_vector_free(fit_parameters);
		gsl_multifit_fdfsolver_free(fit_iterator);
		throw std::bad_alloc();
	}
	
	f.f = &FitFunction_EllipsoidalGaussian;
	f.df = &Jacobian_EllipsoidalGaussian;
	f.fdf = &FitFunctionAndJacobian_EllipsoidalGaussian;
	f.n = number_of_intensities;
	f.p = 7;
	f.params = (void *)&fitData;
	
	
	// iterate over all the determined positions
	std::list<position>::iterator it = positions->begin();
	while (it != positions->end()) {
		iterations = 0;
		
		amplitude = (*it).get_intensity();
		x0_initial = (*it).get_x();
		y0_initial = (*it).get_y();
		correlation_initial = 0.0;
		background = (*it).get_background();
		
		x_offset = floor(x0_initial - (double)cutoff_radius);
		y_offset = floor(y0_initial - (double)cutoff_radius);
		x_max = floor(x0_initial + (double)cutoff_radius);
		y_max = floor(y0_initial + (double)cutoff_radius);
		
		if ((x_offset < 0) || (x_max > (xSize - 1)) || (y_offset < 0) || (y_max > (ySize - 1))) {	// this position is too close to the edge of the image
			// we cannot include it
			it = positions->erase(it);
			continue;
		}
		
		for (size_t j = x_offset; j <= x_max; j++) {
			for (size_t k = y_offset; k <= y_max; k++) {
				(*image_subset)(j - x_offset, k - y_offset) = (*image)(j, k);
			}
		}
		
		fitData.xOffset = x_offset;
		fitData.yOffset = y_offset;
		fitData.imageSubset = image_subset;
		fitData.sigma = sigma;
		
		// provide the initial parameters
		gsl_vector_set(fit_parameters, 0, amplitude);
		gsl_vector_set(fit_parameters, 1, this->initialPSFWidth);
		gsl_vector_set(fit_parameters, 2, this->initialPSFWidth);
		gsl_vector_set(fit_parameters, 3, x0_initial);
		gsl_vector_set(fit_parameters, 4, y0_initial);
		gsl_vector_set(fit_parameters, 5, correlation_initial);
		gsl_vector_set(fit_parameters, 6, background);
		
		// set the solver
		gsl_multifit_fdfsolver_set(fit_iterator, &f, fit_parameters);
		
		// run the iterations
		do {
			iterations++;
			status = gsl_multifit_fdfsolver_iterate(fit_iterator);
			if (status != 0)
				break;
			
			status = gsl_multifit_test_delta(fit_iterator->dx, fit_iterator->x, 1e-4, 1e-4);
		} while ((status == GSL_CONTINUE) && (iterations < 200));
		
		if (gsl_vector_get(fit_iterator->x, 0) <= 0) {	// reject fits that have negative amplitudes
			it = positions->erase(it);
			continue;
		}
		
		// are the reported positions within the window?
		if ((gsl_vector_get(fit_iterator->x, 3) < x_offset) || (gsl_vector_get(fit_iterator->x, 3) > x_max) || (gsl_vector_get(fit_iterator->x, 4) < y_offset) || (gsl_vector_get(fit_iterator->x, 4) > y_max)) {
			// the reported positions are not within the window, we should reject them
			it = positions->erase(it);
			continue;
		}
		
		// are the fitted coordinates close enough to the initial guess?
		if ((fabs(gsl_vector_get(fit_iterator->x, 3) - x0_initial) > 2.0 * initialPSFWidth) || (fabs(gsl_vector_get(fit_iterator->x, 4) - y0_initial) > 2.0 * initialPSFWidth)) {
			it = positions->erase(it);
			continue;
		}
		
		// is the fitted amplitude close enough to the initial value to be trusted?
		if ((gsl_vector_get(fit_iterator->x, 0) < amplitude / 2.0) || (gsl_vector_get(fit_iterator->x, 0) > amplitude * 1.5)) {
			// the output fit amplitude is more than a factor of two different from the initial value, drop this point
			it = positions->erase(it);
			continue;
		}
		
		// are the standard deviations in the x and y direction close enough?
		if ((0.75 * gsl_vector_get(fit_iterator->x, 1) > gsl_vector_get(fit_iterator->x, 2)) || (1.25 * gsl_vector_get(fit_iterator->x, 1) < gsl_vector_get(fit_iterator->x, 2))) {
			// the standard deviations in x and y differ by more than 25%, drop this point
			it = positions->erase(it);
			continue;
		}
		
		if ((gsl_vector_get(fit_iterator->x, 1) < initialPSFWidth / 2.0) || (gsl_vector_get(fit_iterator->x, 1) > initialPSFWidth * 1.5)) {
			// the output fit width is more than a factor of two different from the initial value, drop this point
			it = positions->erase(it);
			continue;
		}
		
		// calculate the covariance matrix
		gsl_multifit_covar(fit_iterator->J, 0.0, covarianceMatrix);
		chi = gsl_blas_dnrm2(fit_iterator->f);
		degreesOfFreedom = (2 * cutoff_radius - 1) * (2 * cutoff_radius - 1) - 5;
		c = GSL_MAX_DBL(1, chi / sqrt(degreesOfFreedom));
		
		// store the fitted parameters
		
		localizationResult->xWidth = gsl_vector_get(fit_iterator->x, 1);
		localizationResult->yWidth = gsl_vector_get(fit_iterator->x, 2);
		localizationResult->xPosition = gsl_vector_get(fit_iterator->x, 3);
		localizationResult->yPosition = gsl_vector_get(fit_iterator->x, 4);
		localizationResult->correlation = gsl_vector_get(fit_iterator->x, 5);
		localizationResult->background = gsl_vector_get(fit_iterator->x, 6);
		
		localizationResult->xWidthDeviation = c * sqrt(gsl_matrix_get(covarianceMatrix, 1, 1));
		localizationResult->yWidthDeviation = c * sqrt(gsl_matrix_get(covarianceMatrix, 2, 2));
		localizationResult->xPositionDeviation = c * sqrt(gsl_matrix_get(covarianceMatrix, 3, 3));
		localizationResult->yPositionDeviation = c * sqrt(gsl_matrix_get(covarianceMatrix, 4, 4));
		localizationResult->correlationDeviation = c * sqrt(gsl_matrix_get(covarianceMatrix, 5, 5));
		localizationResult->backgroundDeviation = c * sqrt(gsl_matrix_get(covarianceMatrix, 6, 6));
		
		localizationResult->integral = 2 * PI * gsl_vector_get(fit_iterator->x, 0) * sqrt(1 - localizationResult->correlation 
										* localizationResult->correlation) * localizationResult->xWidth * localizationResult->yWidth;
		localizationResult->integralDeviation = 2 * M_PI * gsl_vector_get(fit_iterator->x, 0) * sqrt(1 - localizationResult->correlation * localizationResult->correlation)
												* localizationResult->xWidth * localizationResult->yWidth * sqrt(c * sqrt(gsl_matrix_get(covarianceMatrix, 0, 0))
												* c * sqrt(gsl_matrix_get(covarianceMatrix, 0, 0)) / gsl_vector_get(fit_iterator->x, 0) / gsl_vector_get(fit_iterator->x, 0)
												+ localizationResult->correlation * localizationResult->correlation * localizationResult->correlationDeviation
												* localizationResult->correlationDeviation / (1 - localizationResult->correlation * localizationResult->correlation)
												/ (1 - localizationResult->correlation * localizationResult->correlation) + localizationResult->xWidthDeviation
												* localizationResult->xWidthDeviation / localizationResult->xWidth / localizationResult->xWidth
												+ localizationResult->yWidthDeviation * localizationResult->yWidthDeviation / localizationResult->yWidth
												/ localizationResult->yWidth);
		
		fitted_positions->addPosition(localizationResult);
		++it;
	}
	
	gsl_multifit_fdfsolver_free(fit_iterator);
	gsl_vector_free(fit_parameters);
	gsl_matrix_free(covarianceMatrix);
	
	return fitted_positions;
	
}

boost::shared_ptr<LocalizedPositionsContainer> FitPositions_MLEwG::fit_positions(const boost::shared_ptr<Eigen::MatrixXd> image, boost::shared_ptr<std::list<position> > positions) {
	
	// some safety checks
	if (positions->size() == 0) {
		// if no positions were found then there is no reason to run the analysis
		// we need to catch this here and not upstream since then we can return an appropriate
		// instance of LocalizedPositionsContainer
		return boost::shared_ptr<LocalizedPositionsContainer_2DGauss> (new LocalizedPositionsContainer_2DGauss());
	}
	
	size_t size_of_subset = 2 * cutoff_radius + 1;
	double x_offset, y_offset, x_max, y_max;
	size_t xSize = image->size1();
	size_t ySize = image->size2();
	
	double x0_initial, y0_initial, integral, background;
	long iterations = 0;
	int status;
	double size;
	double dummy;
		
	boost::shared_ptr<Eigen::MatrixXd> image_subset;
	boost::shared_ptr<LocalizedPositionsContainer_MLEwG> fitted_positions (new LocalizedPositionsContainer_MLEwG());
	boost::shared_ptr<LocalizedPosition_MLEwG> localizationResult (new LocalizedPosition_MLEwG());
	
	image_subset = boost::shared_ptr<Eigen::MatrixXd> (new ublas::matrix<double>(size_of_subset, size_of_subset));
	
	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer *fit_iterator = NULL;
	gsl_multimin_function minex_func;
	
	measured_data_Gauss_fits fitData;
	
	gsl_vector *fit_parameters = gsl_vector_alloc(5);
	if (fit_parameters == NULL) {
		throw std::bad_alloc();
	}
	gsl_vector *stepSizes = gsl_vector_alloc(5);
	if (stepSizes == NULL) {
		gsl_vector_free(fit_parameters);
		throw std::bad_alloc();
	}
	
	minex_func.n = 5;
	minex_func.f = MinimizationFunction_MLEwG;
	minex_func.params = (void *)&fitData;
	
	fit_iterator = gsl_multimin_fminimizer_alloc (T, 5);
	if (fit_iterator == NULL) {
		gsl_vector_free(fit_parameters);
		gsl_vector_free(stepSizes);
		throw std::bad_alloc();
	}
	
	// iterate over all the determined positions
	std::list<position>::iterator it = positions->begin();
	while (it != positions->end()) {
		iterations = 0;
		
		integral = (*it).get_intensity() * 2.0 * M_PI * this->initialPSFWidth * this->initialPSFWidth;;
		x0_initial = (*it).get_x();
		y0_initial = (*it).get_y();
		background = (*it).get_background();
		
		x_offset = floor(x0_initial - (double)cutoff_radius);
		y_offset = floor(y0_initial - (double)cutoff_radius);
		x_max = floor(x0_initial + (double)cutoff_radius);
		y_max = floor(y0_initial + (double)cutoff_radius);
		
		if ((x_offset < 0) || (x_max > (xSize - 1)) || (y_offset < 0) || (y_max > (ySize - 1))) {	// this position is too close to the edge of the image
			// we cannot include it
			it = positions->erase(it);
			continue;
		}
		
		for (size_t j = x_offset; j <= x_max; j++) {
			for (size_t k = y_offset; k <= y_max; k++) {
				(*image_subset)(j - x_offset, k - y_offset) = (*image)(j, k);
			}
		}
		
		fitData.xOffset = x_offset;
		fitData.yOffset = y_offset;
		fitData.imageSubset = image_subset;
		
		// provide the initial parameters
		gsl_vector_set(fit_parameters, 0, x0_initial);
		gsl_vector_set(fit_parameters, 1, y0_initial);
		gsl_vector_set(fit_parameters, 2, this->initialPSFWidth);
		gsl_vector_set(fit_parameters, 3, background);
		gsl_vector_set(fit_parameters, 4, integral);
		
		// set the step sizes to 1
		gsl_vector_set_all(stepSizes, 1.0);
		// set the solver
		gsl_multimin_fminimizer_set(fit_iterator, &minex_func, fit_parameters, stepSizes);
		
		// iterate
		do {
			iterations++;
			status = gsl_multimin_fminimizer_iterate(fit_iterator);
			if (status != 0)
				break;
			
			size = gsl_multimin_fminimizer_size (fit_iterator);
			status = gsl_multimin_test_size (size, 1e-2);
		} while ((status == GSL_CONTINUE) && (iterations < 200));
		
		if (gsl_vector_get(fit_iterator->x, 4) <= 0) {	// reject fits that have negative amplitudes
			it = positions->erase(it);
			continue;
		}
		
		// are the reported positions within the window?
		if ((gsl_vector_get(fit_iterator->x, 0) < x_offset) || (gsl_vector_get(fit_iterator->x, 0) > x_max) || (gsl_vector_get(fit_iterator->x, 1) < y_offset) || (gsl_vector_get(fit_iterator->x, 1) > y_max)) {
			// the reported positions are not within the window, we should reject them
			it = positions->erase(it);
			continue;
		}
		
		// are the fitted coordinates close enough to the initial guess?
		if ((fabs(gsl_vector_get(fit_iterator->x, 0) - x0_initial) > 2.0 * initialPSFWidth) || (fabs(gsl_vector_get(fit_iterator->x, 1) - y0_initial) > 2.0 * initialPSFWidth)) {
			it = positions->erase(it);
			continue;
		}
		
		// is the amplitude close enough to the initial value to be trusted?
		if ((gsl_vector_get(fit_iterator->x, 4) < integral / 2.0) || (gsl_vector_get(fit_iterator->x, 4) > integral * 1.5)) {
			// the output fit integral is more than a factor of two different from the initial value, drop this point
			it = positions->erase(it);
			continue;
		}
		
		if ((gsl_vector_get(fit_iterator->x, 2) < initialPSFWidth / 2.0) || (gsl_vector_get(fit_iterator->x, 2) > initialPSFWidth * 1.5)) {
			// the output fit width is more than a factor of two different from the initial value, drop this point
			dummy = gsl_vector_get(fit_iterator->x, 2);
			it = positions->erase(it);
			continue;
		}
		
		// store the fitted parameters
		localizationResult->xPosition = gsl_vector_get(fit_iterator->x, 0);
		localizationResult->yPosition = gsl_vector_get(fit_iterator->x, 1);
		localizationResult->width = gsl_vector_get(fit_iterator->x, 2);
		localizationResult->background = gsl_vector_get(fit_iterator->x, 3);
		localizationResult->integral = gsl_vector_get(fit_iterator->x, 4);
		localizationResult->positionDeviation = sqrt(CalculateMLEwGVariance(initialPSFWidth, localizationResult->integral, localizationResult->background));
		
		fitted_positions->addPosition(localizationResult);
		++it;
	}
	
	gsl_multimin_fminimizer_free(fit_iterator);
	gsl_vector_free(fit_parameters);
	gsl_vector_free(stepSizes);
	
	return fitted_positions;
	
}

boost::shared_ptr<LocalizedPositionsContainer> FitPositionsMultiplication::fit_positions(const boost::shared_ptr<Eigen::MatrixXd> image,
																						 boost::shared_ptr<std::list<position> > positions) {
	
	// some safety checks
	if (positions->size() == 0) {
		// if no positions were found then there is no reason to run the analysis
		// we need to catch this here and not upstream since then we can return an appropriate
		// instance of LocalizedPositionsContainer
		return boost::shared_ptr<LocalizedPositionsContainer_Multiplication> (new LocalizedPositionsContainer_Multiplication());
	}
	
	size_t size_of_subset = 2 * cutoff_radius + 1;
	size_t xSize = image->size1();
	size_t ySize = image->size2();
	double x_offset, y_offset, x_max, y_max;
	
	double x0_initial, y0_initial, amplitude, background;
	size_t iterations = 0;
	int converged;
	
	double convergence_treshold_squared = convergence_threshold * convergence_threshold;
	double delta_squared = 10 * convergence_treshold_squared;	// this test the convergence of the position determined by the iteration
	// it is the distance between (xn-1, yn-1) and (xn, yn)
	// we initialize it to a value well over the treshold so that we will run at least two iterations
	double previous_position_x, previous_position_y;
	double current_x, current_y;
	
	boost::shared_ptr<Eigen::MatrixXd> image_subset;
	boost::shared_ptr<Eigen::MatrixXd> image_subset_mask;
	boost::shared_ptr<LocalizedPositionsContainer_Multiplication> fitted_positions (new LocalizedPositionsContainer_Multiplication());
	boost::shared_ptr<LocalizedPosition_Multiplication> localizationResult (new LocalizedPosition_Multiplication());
	
	image_subset = boost::shared_ptr<Eigen::MatrixXd>(new ublas::matrix<double>(size_of_subset, size_of_subset));
	image_subset_mask = boost::shared_ptr<Eigen::MatrixXd> (new ublas::matrix<double>(size_of_subset, size_of_subset));
	
	std::list<position>::iterator it = positions->begin();
	while (it != positions->end()) {
		amplitude = (*it).get_intensity();
		x0_initial = (*it).get_x();
		y0_initial = (*it).get_y();
		background = (*it).get_background();
		
		x_offset = floor(x0_initial - (double)cutoff_radius);
		y_offset = floor(y0_initial - (double)cutoff_radius);
		x_max = floor(x0_initial + (double)cutoff_radius);
		y_max = floor(y0_initial + (double)cutoff_radius);
		
		if ((x_offset < 0) || (x_max > (xSize - 1)) || (y_offset < 0) || (y_max > (ySize - 1))) {	// this position is too close to the edge of the image, we cannot include it
			it = positions->erase(it);
			continue;
		}
		
		for (size_t j = x_offset; j <= x_max; ++j) {
			for (size_t k = y_offset; k <= y_max; ++k) {
				(*image_subset)(j - x_offset, k - y_offset) = (*image)(j, k);
			}
		}
		
		iterations = 0;
		
		current_x = x0_initial - x_offset;	// correct the x- and y-values for the fact that we analyze in a subset of the image rather than the complete frame
		current_y = y0_initial - y_offset;
		
		converged = 1;
		while (delta_squared > convergence_treshold_squared) {
			previous_position_x = current_x;
			previous_position_y = current_y;
			
			++iterations;
			
			if (iterations > 100) {	// the multiplication is not converging, we should stop
				converged = 0;
				break;
			}
			
			multiply_with_gaussian(image_subset, image_subset_mask, current_x, current_y, initialPSFWidth, background, amplitude);
			determine_x_y_position(image_subset_mask, current_x, current_y);
			
			if (iterations == 1)	// this is the first iteration, we should not check for termination
				continue;
			
			delta_squared = (current_x - previous_position_x) * (current_x - previous_position_x) + (current_y - previous_position_y) * (current_y - previous_position_y);
		}
		
		if (converged == 0)
			continue;
		
		delta_squared = 10 * convergence_treshold_squared;
		
		localizationResult->width = initialPSFWidth;
		localizationResult->xPosition = (double)current_x + (double)x_offset;
		localizationResult->yPosition = (double)current_y + (double)y_offset;
		
		fitted_positions->addPosition(localizationResult);
		++it;
	}
	
	return fitted_positions;
}


int FitPositionsMultiplication::multiply_with_gaussian(boost::shared_ptr<Eigen::MatrixXd> original_image, boost::shared_ptr<Eigen::MatrixXd> masked_image, double x, double y, 
													   double std_dev, double background, double amplitude) {
	// we will replace the contents of masked_image with the multiplication of original_image and a gaussian centered at position (x,y)
	
	size_t x_size = masked_image->size1();
	size_t y_size = masked_image->size2();
	
	double gaussian_value, distance_squared;
	
	if ((original_image->size1() != x_size) || (original_image->size2() != y_size)) {
		throw DIMENSIONS_SHOULD_BE_EQUAL(std::string("Matrix dimensions are not equal in FitPositionsMultiplication::multiply_with_gaussian"));
	}
	
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {
			distance_squared = (x - (double)i) * (x - (double)i) + (y - (double)j) * (y - (double)j);
			
			gaussian_value = amplitude * exp(- distance_squared / (2 * std_dev * std_dev)) + background;
			
			(*masked_image)(i, j) = gaussian_value * (*original_image)(i, j);
		}
	}
	
	return 0;
}


int FitPositionsMultiplication::determine_x_y_position(boost::shared_ptr<Eigen::MatrixXd> masked_image, double &x, double &y) {
	// based on eq (3) in Thompson Biophys J 2002
	
	size_t x_size = (size_t)masked_image->size1();
	size_t y_size = (size_t)masked_image->size2();
	
	double numerator_x = 0, denominator = 0;
	double numerator_y = 0;
	
	// start with determining the x-position
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {
			numerator_x += (double)i * (*masked_image)(i, j);
			numerator_y += (double)j * (*masked_image)(i, j);
			denominator += (*masked_image)(i, j);
		}
	}
	
	x = numerator_x / denominator;
	
	y = numerator_y / denominator;
	
	return 0;
}


boost::shared_ptr<LocalizedPositionsContainer> FitPositionsCentroid::fit_positions(const boost::shared_ptr<Eigen::MatrixXd> image,
																				   boost::shared_ptr<std::list<position> > positions) {
	
	// some safety checks
	if (positions->size() == 0) {
		// if no positions were found then there is no reason to run the analysis
		// we need to catch this here and not upstream since then we can return an appropriate
		// instance of LocalizedPositionsContainer
		return boost::shared_ptr<LocalizedPositionsContainer_Centroid> (new LocalizedPositionsContainer_Centroid());
	}
	
	size_t xSize = image->size1();
	size_t ySize = image->size2();
	double x_offset, y_offset, x_max, y_max;
	
	size_t x0_initial, y0_initial;
	double current_x, current_y;
	double denominator;
	
	boost::shared_ptr<LocalizedPositionsContainer_Centroid> fitted_positions (new LocalizedPositionsContainer_Centroid());
	boost::shared_ptr<LocalizedPosition_Centroid> localizationResult (new LocalizedPosition_Centroid());
	
	std::list<position>::iterator it = positions->begin();
	while (it != positions->end()) {
		x0_initial = (*it).get_x();
		y0_initial = (*it).get_y();
		current_x = 0;
		current_y = 0;
		denominator = 0;
		
		x_offset = floor(x0_initial - (double)cutoff_radius);
		y_offset = floor(y0_initial - (double)cutoff_radius);
		x_max = floor(x0_initial + (double)cutoff_radius);
		y_max = floor(y0_initial + (double)cutoff_radius);
		
		if ((x_offset < 0) || (x_max > (xSize - 1)) || (y_offset < 0) || (y_max > (ySize - 1))) {	// the point is too close to the edge
			it = positions->erase(it);
			continue;
		}
		
		for (size_t j = x_offset; j <= x_max; ++j) {
			for (size_t k = y_offset; k <= y_max; ++k) {
				current_x += (double)j * (*image)(j, k);
				current_y += (double)k * (*image)(j, k);
				denominator += (*image)(j, k);
			}
		}
		
		current_x /= denominator;
		current_y /= denominator;		
		
		localizationResult->xPosition = current_x;
		localizationResult->yPosition = current_y;
		
		fitted_positions->addPosition(localizationResult);
		++it;
	}
	
	return fitted_positions;
}

int FitFunction_SymmetricGaussian(const gsl_vector *params, void *fitData_rhs, gsl_vector *deviations) {
	// params contains the current values of the parameters - amplitude, width, etc.
	// fitData is an object that contains the experimental data
	
	
	measured_data_Gauss_fits *fitDataLocal = (measured_data_Gauss_fits *)fitData_rhs;
	boost::shared_ptr<Eigen::MatrixXd> imageSubset = fitDataLocal->imageSubset;
	
	size_t xSize = imageSubset->size1();
	size_t ySize = imageSubset->size2();
	size_t arrayOffset = 0;
	double xOffset = fitDataLocal->xOffset;
	double yOffset = fitDataLocal->yOffset;
	double sigma = fitDataLocal->sigma;
	
	double amplitude = gsl_vector_get(params, 0);
	double r = gsl_vector_get(params, 1);
	double x0 = gsl_vector_get(params, 2);
	double y0 = gsl_vector_get(params, 3);
	double offset = gsl_vector_get(params, 4);
	
	double x,y, function_value, square_deviation;
	
	if (r == 0) {
		return GSL_FAILURE;
	}
	
	for (size_t i = 0; i < xSize; ++i) {
		for (size_t j = 0; j < ySize; ++j) {
			x = xOffset + (double)i;
			y = yOffset + (double)j;
			function_value = offset + amplitude * exp(- (((x0 - x)/ (SQRT2 * r)) * ((x0 - x)/ (SQRT2 * r)) + ((y0 - y) / (SQRT2 * r)) * ((y0 - y) / (SQRT2 * r))));
			square_deviation = (function_value - (*imageSubset)(i, j)) / sigma;
			gsl_vector_set(deviations, arrayOffset, square_deviation);
			++arrayOffset;
		}
	}
	return GSL_SUCCESS;
}

int FitFunction_FixedWidthGaussian(const gsl_vector *params, void *fitData_rhs, gsl_vector *deviations) {
	// params contains the current values of the parameters - amplitude, width, etc.
	// fitData is an object that contains the experimental data
	
	
	measured_data_Gauss_fits *fitDataLocal = (measured_data_Gauss_fits *)fitData_rhs;
	boost::shared_ptr<Eigen::MatrixXd> imageSubset = fitDataLocal->imageSubset;
	
	size_t xSize = imageSubset->size1();
	size_t ySize = imageSubset->size2();
	size_t arrayOffset = 0;
	double xOffset = fitDataLocal->xOffset;
	double yOffset = fitDataLocal->yOffset;
	double sigma = fitDataLocal->sigma;
	double r = fitDataLocal->width;
	
	double amplitude = gsl_vector_get(params, 0);
	double x0 = gsl_vector_get(params, 1);
	double y0 = gsl_vector_get(params, 2);
	double offset = gsl_vector_get(params, 3);
	
	double x,y, function_value, square_deviation;
	
	if (r == 0) {
		return GSL_FAILURE;
	}
	
	for (size_t i = 0; i < xSize; ++i) {
		for (size_t j = 0; j < ySize; ++j) {
			x = xOffset + (double)i;
			y = yOffset + (double)j;
			function_value = offset + amplitude * exp(- (((x0 - x)/ (SQRT2 * r)) * ((x0 - x)/ (SQRT2 * r)) + ((y0 - y) / (SQRT2 * r)) * ((y0 - y) / (SQRT2 * r))));
			square_deviation = (function_value - (*imageSubset)(i, j)) / sigma;
			gsl_vector_set(deviations, arrayOffset, square_deviation);
			++arrayOffset;
		}
	}
	return GSL_SUCCESS;
}

int FitFunction_EllipsoidalGaussian(const gsl_vector *params, void *fitData_rhs, gsl_vector *deviations) {
	// params contains the current values of the parameters - amplitude, width, etc.
	// fitData is an object that contains the experimental data
	
	
	measured_data_Gauss_fits *fitDataLocal = (measured_data_Gauss_fits *)fitData_rhs;
	boost::shared_ptr<Eigen::MatrixXd> imageSubset = fitDataLocal->imageSubset;
	
	size_t xSize = imageSubset->size1();
	size_t ySize = imageSubset->size2();
	size_t arrayOffset = 0;
	double xOffset = fitDataLocal->xOffset;
	double yOffset = fitDataLocal->yOffset;
	double sigma = fitDataLocal->sigma;
	
	double amplitude = gsl_vector_get(params, 0);
	double sigmaX = gsl_vector_get(params, 1);
	double sigmaY = gsl_vector_get(params, 2);
	double x0 = gsl_vector_get(params, 3);
	double y0 = gsl_vector_get(params, 4);
	double corr = gsl_vector_get(params, 5);
	double offset = gsl_vector_get(params, 6);
	
	if ((sigmaX == 0) || (sigmaY == 0)) {
		return GSL_FAILURE;
	}
	
	double x,y, function_value, square_deviation;
	
	for (size_t i = 0; i < xSize; ++i) {
		for (size_t j = 0; j < ySize; ++j) {
			x = xOffset + (double)i;
			y = yOffset + (double)j;
			function_value = amplitude * exp(- 1.0 / (2.0 * (1.0 - corr * corr)) * ((x - x0) * (x - x0) / sigmaX / sigmaX 
																					+ (y - y0) * (y - y0) / sigmaY / sigmaY 
																					- 2.0 * corr * (x - x0) * (y - y0) / sigmaX / sigmaY)) + offset;
			square_deviation = (function_value - (*imageSubset)(i, j)) / sigma;
			gsl_vector_set(deviations, arrayOffset, square_deviation);
			++arrayOffset;
		}
	}
	return GSL_SUCCESS;
}

int Jacobian_SymmetricGaussian(const gsl_vector *params, void *fitData_rhs, gsl_matrix *jacobian) {
	measured_data_Gauss_fits *fitDataLocal = (measured_data_Gauss_fits *)fitData_rhs;
	boost::shared_ptr<Eigen::MatrixXd> imageSubset = fitDataLocal->imageSubset;
	
	size_t xSize = imageSubset->size1();
	size_t ySize = imageSubset->size2();
	size_t arrayOffset = 0;
	double xOffset = fitDataLocal->xOffset;
	double yOffset = fitDataLocal->yOffset;
	double sigma = fitDataLocal->sigma;
	
	double amplitude = gsl_vector_get(params, 0);
	double r = gsl_vector_get(params, 1);
	double x0 = gsl_vector_get(params, 2);
	double y0 = gsl_vector_get(params, 3);
	// double offset = gsl_vector_get(params, 4);
	
	double x,y, exp_factor;
	double dfdA, dfdr, dfdx0, dfdy0,dfdoffset;
	
	if (r == 0) {
		return GSL_FAILURE;
	}
	
	for (size_t i = 0; i < xSize; ++i) {
		for (size_t j = 0; j < ySize; ++j) {
			x = xOffset + (double)i;
			y = yOffset + (double)j;
			
			exp_factor = exp(- ((x0 - x)/ (SQRT2 * r)) * ((x0 - x)/ (SQRT2 * r)) - ((y0 - y) / (SQRT2 * r)) * ((y0 - y) / (SQRT2 * r)));
			
			dfdA = exp_factor / sigma;
			dfdr = (2 * (y - y0) * (y - y0) / (SQRT2 * r) / (SQRT2 * r) / (SQRT2 * r) + 2 * (x - x0) * (x - x0) / (SQRT2 * r) / (SQRT2 * r) / (SQRT2 * r)) * exp_factor * amplitude / sigma;
			dfdx0 = (2 * (x - x0) * exp_factor * amplitude) / ((SQRT2 * r) * (SQRT2 * r) * sigma);
			dfdy0 = (2 * (y - y0) * exp_factor * amplitude) / ((SQRT2 * r) * (SQRT2 * r) * sigma);
			dfdoffset = 1/sigma;
			
			gsl_matrix_set(jacobian, arrayOffset, 0, dfdA);
			gsl_matrix_set(jacobian, arrayOffset, 1, dfdr);
			gsl_matrix_set(jacobian, arrayOffset, 2, dfdx0);
			gsl_matrix_set(jacobian, arrayOffset, 3, dfdy0);
			gsl_matrix_set(jacobian, arrayOffset, 4, dfdoffset);
			++arrayOffset;
		}
	}
	
	return GSL_SUCCESS;
	
}

int Jacobian_FixedWidthGaussian(const gsl_vector *params, void *fitData_rhs, gsl_matrix *jacobian) {
	measured_data_Gauss_fits *fitDataLocal = (measured_data_Gauss_fits *)fitData_rhs;
	boost::shared_ptr<Eigen::MatrixXd> imageSubset = fitDataLocal->imageSubset;
	
	size_t xSize = imageSubset->size1();
	size_t ySize = imageSubset->size2();
	size_t arrayOffset = 0;
	double xOffset = fitDataLocal->xOffset;
	double yOffset = fitDataLocal->yOffset;
	double sigma = fitDataLocal->sigma;
	double r = fitDataLocal->width;
	
	double amplitude = gsl_vector_get(params, 0);
	double x0 = gsl_vector_get(params, 1);
	double y0 = gsl_vector_get(params, 2);
	
	double x,y, exp_factor;
	double dfdA, dfdx0, dfdy0,dfdoffset;
	
	if (r == 0) {
		return GSL_FAILURE;
	}
	
	for (size_t i = 0; i < xSize; ++i) {
		for (size_t j = 0; j < ySize; ++j) {
			x = xOffset + (double)i;
			y = yOffset + (double)j;
			
			exp_factor = exp(- ((x0 - x)/ (SQRT2 * r)) * ((x0 - x)/ (SQRT2 * r)) - ((y0 - y) / (SQRT2 * r)) * ((y0 - y) / (SQRT2 * r)));
			
			dfdA = exp_factor / sigma;
			dfdx0 = (2 * (x - x0) * exp_factor * amplitude) / ((SQRT2 * r) * (SQRT2 * r) * sigma);
			dfdy0 = (2 * (y - y0) * exp_factor * amplitude) / ((SQRT2 * r) * (SQRT2 * r) * sigma);
			dfdoffset = 1/sigma;
			
			gsl_matrix_set(jacobian, arrayOffset, 0, dfdA);
			gsl_matrix_set(jacobian, arrayOffset, 1, dfdx0);
			gsl_matrix_set(jacobian, arrayOffset, 2, dfdy0);
			gsl_matrix_set(jacobian, arrayOffset, 3, dfdoffset);
			++arrayOffset;
		}
	}
	
	return GSL_SUCCESS;
	
}

int Jacobian_EllipsoidalGaussian(const gsl_vector *params, void *fitData_rhs, gsl_matrix *jacobian) {
	measured_data_Gauss_fits *fitDataLocal = (measured_data_Gauss_fits *)fitData_rhs;
	boost::shared_ptr<Eigen::MatrixXd> imageSubset = fitDataLocal->imageSubset;
	
	size_t xSize = imageSubset->size1();
	size_t ySize = imageSubset->size2();
	size_t arrayOffset = 0;
	double xOffset = fitDataLocal->xOffset;
	double yOffset = fitDataLocal->yOffset;
	double sigma = fitDataLocal->sigma;
	
	double amplitude = gsl_vector_get(params, 0);
	double sigmaX = gsl_vector_get(params, 1);
	double sigmaY = gsl_vector_get(params, 2);
	double x0 = gsl_vector_get(params, 3);
	double y0 = gsl_vector_get(params, 4);
	double corr = gsl_vector_get(params, 5);
	
	double x,y, exp_factor, denominator;
	double dfdA, dfdsigmaX, dfdsigmaY, dfdx0, dfdy0, dfdcorr,dfdoffset;
	
	if ((sigmaX == 0) || (sigmaY == 0)) {
		return GSL_FAILURE;
	}
	
	denominator = 2 * (1 - corr * corr) * sigma;
	
	for (size_t i = 0; i < xSize; ++i) {
		for (size_t j = 0; j < ySize; ++j) {
			x = xOffset + (double)i;
			y = yOffset + (double)j;
			
			exp_factor = exp(- 1.0 / (2.0 * (1.0 - corr * corr)) * ((x - x0) * (x - x0) / sigmaX / sigmaX 
																	+ (y - y0) * (y - y0) / sigmaY / sigmaY
																	- 2.0 * corr * (x - x0) * (y - y0) / sigmaX / sigmaY));
			
			dfdA = exp_factor / sigma;
			dfdsigmaX = - amplitude * (2.0 * corr * (x - x0) * (y - y0) / sigmaX / sigmaX / sigmaY - 2.0 * (x - x0) * (x - x0) / sigmaX / sigmaX / sigmaX) * exp_factor / denominator;
			dfdsigmaY = - amplitude * (2.0 * corr * (x - x0) * (y - y0) / sigmaX / sigmaY / sigmaY - 2.0 * (y - y0) * (y - y0) / sigmaY / sigmaY / sigmaY) * exp_factor / denominator;
			dfdx0 = - amplitude * (2 * corr * (y - y0) / sigmaX / sigmaY - 2 * (x - x0) / sigmaX / sigmaX) * exp_factor / denominator;
			dfdy0 = - amplitude * (2 * corr * (x - x0) / sigmaX / sigmaY - 2 * (y - y0) / sigmaY / sigmaY) * exp_factor / denominator;
			dfdcorr = amplitude * ((x - x0) * (y - y0) / (1 - corr * corr) / sigmaX / sigmaY - corr * (- 2 * corr * (x - x0) * (y - y0) / sigmaX / sigmaY 
																									   + (y - y0) * (y - y0) / sigmaY / sigmaY + (x - x0) * (x - x0) / sigmaX / sigmaX) / (1 - corr * corr) / (1 - corr * corr)) * exp_factor / sigma;
			dfdoffset = 1.0 / sigma;
			
			gsl_matrix_set(jacobian, arrayOffset, 0, dfdA);
			gsl_matrix_set(jacobian, arrayOffset, 1, dfdsigmaX);
			gsl_matrix_set(jacobian, arrayOffset, 2, dfdsigmaY);
			gsl_matrix_set(jacobian, arrayOffset, 3, dfdx0);
			gsl_matrix_set(jacobian, arrayOffset, 4, dfdy0);
			gsl_matrix_set(jacobian, arrayOffset, 5, dfdcorr);
			gsl_matrix_set(jacobian, arrayOffset, 6, dfdoffset);
			++arrayOffset;
		}
	}
	
	return GSL_SUCCESS;
	
}

int FitFunctionAndJacobian_SymmetricGaussian(const gsl_vector *params, void *fitData_rhs, gsl_vector *model_values, gsl_matrix *jacobian) {
	measured_data_Gauss_fits *fitDataLocal = (measured_data_Gauss_fits *)fitData_rhs;
	boost::shared_ptr<Eigen::MatrixXd> imageSubset = fitDataLocal->imageSubset;
	
	size_t xSize = imageSubset->size1();
	size_t ySize = imageSubset->size2();
	size_t arrayOffset = 0;
	double xOffset = fitDataLocal->xOffset;
	double yOffset = fitDataLocal->yOffset;
	double sigma = fitDataLocal->sigma;
	
	double amplitude = gsl_vector_get(params, 0);
	double r = gsl_vector_get(params, 1);
	double x0 = gsl_vector_get(params, 2);
	double y0 = gsl_vector_get(params, 3);
	double offset = gsl_vector_get(params, 4);
	
	double x,y, function_value, square_deviation, exp_factor;
	double dfdA, dfdr, dfdx0, dfdy0, dfdoffset;
	
	if (r == 0) {
		return GSL_FAILURE;
	}
	
	for (size_t i = 0; i < xSize; ++i) {
		for (size_t j = 0; j < ySize; ++j) {
			x = xOffset + (double)i;
			y = yOffset + (double)j;
			
			exp_factor = exp(- ((x0 - x)/ (SQRT2 * r)) * ((x0 - x)/ (SQRT2 * r)) - ((y0 - y) / (SQRT2 * r)) * ((y0 - y) / (SQRT2 * r)));
			
			function_value = offset + amplitude * exp_factor;
			square_deviation = (function_value - (*imageSubset)(i, j)) / sigma;
			
			dfdA = exp_factor / sigma;
			dfdr = (2 * (y - y0) * (y - y0) / (SQRT2 * r) / (SQRT2 * r) / (SQRT2 * r) + 2 * (x - x0) * (x - x0) / (SQRT2 * r) / (SQRT2 * r) / (SQRT2 * r)) * exp_factor * amplitude / sigma;
			dfdx0 = (2 * (x - x0) * exp_factor * amplitude) / ((SQRT2 * r) * (SQRT2 * r) * sigma);
			dfdy0 = (2 * (y - y0) * exp_factor * amplitude) / ((SQRT2 * r) * (SQRT2 * r) * sigma);
			dfdoffset = 1/sigma;
			
			gsl_vector_set(model_values, arrayOffset, square_deviation);
			
			gsl_matrix_set(jacobian, arrayOffset, 0, dfdA);
			gsl_matrix_set(jacobian, arrayOffset, 1, dfdr);
			gsl_matrix_set(jacobian, arrayOffset, 2, dfdx0);
			gsl_matrix_set(jacobian, arrayOffset, 3, dfdy0);
			gsl_matrix_set(jacobian, arrayOffset, 4, dfdoffset);
			
			++arrayOffset;
		}
	}
	return GSL_SUCCESS;
}

int FitFunctionAndJacobian_FixedWidthGaussian(const gsl_vector *params, void *measured_intensities_struct, gsl_vector *model_values, gsl_matrix *jacobian) {
	int result;
	result = FitFunction_FixedWidthGaussian(params, measured_intensities_struct, model_values);
	if (result != GSL_SUCCESS)
		return result;
	result = Jacobian_FixedWidthGaussian(params, measured_intensities_struct, jacobian);
	if (result != GSL_SUCCESS)
		return result;
	
	return GSL_SUCCESS;
}

int FitFunctionAndJacobian_EllipsoidalGaussian(const gsl_vector *params, void *measured_intensities_struct, gsl_vector *model_values, gsl_matrix *jacobian) {
	int result;
	result = FitFunction_EllipsoidalGaussian(params, measured_intensities_struct, model_values);
	if (result != GSL_SUCCESS)
		return result;
	result = Jacobian_EllipsoidalGaussian(params, measured_intensities_struct, jacobian);
	if (result != GSL_SUCCESS)
		return result;
	
	return GSL_SUCCESS;
}

double MinimizationFunction_MLEwG(const gsl_vector *fittedParams, void *fitData_rhs) {
	measured_data_Gauss_fits *fitDataLocal = (measured_data_Gauss_fits *)fitData_rhs;
	boost::shared_ptr<Eigen::MatrixXd> imageSubset = fitDataLocal->imageSubset;
	
	size_t xSize = imageSubset->size1();
	size_t ySize = imageSubset->size2();
	double xOffset = fitDataLocal->xOffset;
	double yOffset = fitDataLocal->yOffset;
	
	double x0 = gsl_vector_get(fittedParams, 0);
	double y0 = gsl_vector_get(fittedParams, 1);
	double stdDev = gsl_vector_get(fittedParams, 2);
	double background = gsl_vector_get(fittedParams, 3);
	double nPhotons = gsl_vector_get(fittedParams, 4);
	
	double expectationValue, recordedSignal, summedLikelihood, x, y;
	
	summedLikelihood = 0.0;
	for (size_t i = 0; i < xSize; ++i) {
		for (size_t j = 0; j < ySize; ++j) {
			x = xOffset + (double)i;
			y = yOffset + (double)j;
			recordedSignal = (*imageSubset)(i, j);
			expectationValue = nPhotons / (2 * PI * stdDev * stdDev) * exp(-((x - x0) * (x - x0) + (y - y0) * (y - y0)) / stdDev / stdDev) + background;
			if (recordedSignal >= 0.0)
				summedLikelihood += - expectationValue + recordedSignal * log(expectationValue) - gsl_sf_lngamma(recordedSignal + 1.0);
		}
	}
	return (-1.0 * summedLikelihood);	// make the number negative since what we really want is maximization
}

double CalculateMLEwGVariance(double PSFWidth, double nPhotons, double background) {
	// based on eq 6 in the Mortensen paper since eq 5 does not seem very stable
	double variance, sigmaA;
	
	sigmaA = PSFWidth * PSFWidth + 1.0 / 12.0;
	
	variance = sigmaA * sigmaA / nPhotons * (16.0 / 9.0 + 8 * M_PI * PSFWidth * PSFWidth * background / nPhotons);
	
	return variance;
}

double MLEwGIntegrand(double t, void *params_rhs) {
	double *params = (double *)params_rhs;
	
	double sigmaStar = params[0];
	double nPhotons = params[1];
	double background = params[2];
	
	return (log(t) / (1 + t / (2 * M_PI * sigmaStar * sigmaStar * background / nPhotons)));
}




