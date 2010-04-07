/*
 *  PALM_analysis_Localization.cpp
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 27/09/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "PALM_analysis_Localization.h"

int Gauss_2D_fit_function(const gsl_vector *params, void *fitData_rhs, gsl_vector *deviations) {
	// params contains the current values of the parameters - amplitude, width, etc.
	// fitData is an object that contains the experimental data
	
	
	measured_data_Gauss_fits *fitDataLocal = (measured_data_Gauss_fits *)fitData_rhs;
	boost::shared_ptr<PALMMatrix<double> > imageSubset = fitDataLocal->imageSubset;
	
	size_t xSize = imageSubset->getXSize();
	size_t ySize = imageSubset->getYSize();
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
			square_deviation = (function_value - imageSubset->get(i, j)) / sigma;
			gsl_vector_set(deviations, arrayOffset, square_deviation);
			++arrayOffset;
		}
	}
	return GSL_SUCCESS;
}

int Gauss_2D_fit_function_FixedWidth(const gsl_vector *params, void *fitData_rhs, gsl_vector *deviations) {
	// params contains the current values of the parameters - amplitude, width, etc.
	// fitData is an object that contains the experimental data
	
	
	measured_data_Gauss_fits *fitDataLocal = (measured_data_Gauss_fits *)fitData_rhs;
	boost::shared_ptr<PALMMatrix<double> > imageSubset = fitDataLocal->imageSubset;
	
	size_t xSize = imageSubset->getXSize();
	size_t ySize = imageSubset->getYSize();
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
			square_deviation = (function_value - imageSubset->get(i, j)) / sigma;
			gsl_vector_set(deviations, arrayOffset, square_deviation);
			++arrayOffset;
		}
	}
	return GSL_SUCCESS;
}

int Ellipsoidal_Gauss_2D_fit_function(const gsl_vector *params, void *fitData_rhs, gsl_vector *deviations) {
	// params contains the current values of the parameters - amplitude, width, etc.
	// fitData is an object that contains the experimental data
	
	
	measured_data_Gauss_fits *fitDataLocal = (measured_data_Gauss_fits *)fitData_rhs;
	boost::shared_ptr<PALMMatrix<double> > imageSubset = fitDataLocal->imageSubset;
	
	size_t xSize = imageSubset->getXSize();
	size_t ySize = imageSubset->getYSize();
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
			square_deviation = (function_value - imageSubset->get(i, j)) / sigma;
			gsl_vector_set(deviations, arrayOffset, square_deviation);
			++arrayOffset;
		}
	}
	return GSL_SUCCESS;
}

int Gauss_2D_fit_function_Jacobian(const gsl_vector *params, void *fitData_rhs, gsl_matrix *jacobian) {
	measured_data_Gauss_fits *fitDataLocal = (measured_data_Gauss_fits *)fitData_rhs;
	boost::shared_ptr<PALMMatrix<double> > imageSubset = fitDataLocal->imageSubset;
	
	size_t xSize = imageSubset->getXSize();
	size_t ySize = imageSubset->getXSize();
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

int Gauss_2D_fit_function_Jacobian_FixedWidth(const gsl_vector *params, void *fitData_rhs, gsl_matrix *jacobian) {
	measured_data_Gauss_fits *fitDataLocal = (measured_data_Gauss_fits *)fitData_rhs;
	boost::shared_ptr<PALMMatrix<double> > imageSubset = fitDataLocal->imageSubset;
	
	size_t xSize = imageSubset->getXSize();
	size_t ySize = imageSubset->getXSize();
	size_t arrayOffset = 0;
	double xOffset = fitDataLocal->xOffset;
	double yOffset = fitDataLocal->yOffset;
	double sigma = fitDataLocal->sigma;
	double r = fitDataLocal->width;
	
	double amplitude = gsl_vector_get(params, 0);
	double x0 = gsl_vector_get(params, 2);
	double y0 = gsl_vector_get(params, 3);
	// double offset = gsl_vector_get(params, 4);
	
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
			// dfdr = (2 * (y - y0) * (y - y0) / r / r / r + 2 * (x - x0) * (x - x0) / r / r / r) * exp_factor * amplitude / sigma;
			dfdx0 = (2 * (x - x0) * exp_factor * amplitude) / ((SQRT2 * r) * (SQRT2 * r) * sigma);
			dfdy0 = (2 * (y - y0) * exp_factor * amplitude) / ((SQRT2 * r) * (SQRT2 * r) * sigma);
			dfdoffset = 1/sigma;
			
			gsl_matrix_set(jacobian, arrayOffset, 0, dfdA);
			// gsl_matrix_set(jacobian, arrayOffset, 1, dfdr);
			gsl_matrix_set(jacobian, arrayOffset, 1, dfdx0);
			gsl_matrix_set(jacobian, arrayOffset, 2, dfdy0);
			gsl_matrix_set(jacobian, arrayOffset, 3, dfdoffset);
			++arrayOffset;
		}
	}
	
	return GSL_SUCCESS;
	
}

int Ellipsoidal_Gauss_2D_fit_function_Jacobian(const gsl_vector *params, void *fitData_rhs, gsl_matrix *jacobian) {
	measured_data_Gauss_fits *fitDataLocal = (measured_data_Gauss_fits *)fitData_rhs;
	boost::shared_ptr<PALMMatrix<double> > imageSubset = fitDataLocal->imageSubset;
	
	size_t xSize = imageSubset->getXSize();
	size_t ySize = imageSubset->getXSize();
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

int Gauss_2D_fit_function_and_Jacobian(const gsl_vector *params, void *measured_intensities_struct, gsl_vector *model_values, gsl_matrix *jacobian) {
	int result;
	result = Gauss_2D_fit_function(params, measured_intensities_struct, model_values);
	if (result != GSL_SUCCESS)
		return result;
	result = Gauss_2D_fit_function_Jacobian(params, measured_intensities_struct, jacobian);
	if (result != GSL_SUCCESS)
		return result;
	return GSL_SUCCESS;
}

int Gauss_2D_fit_function_and_Jacobian_FixedWidth(const gsl_vector *params, void *measured_intensities_struct, gsl_vector *model_values, gsl_matrix *jacobian) {
	int result;
	result = Gauss_2D_fit_function_FixedWidth(params, measured_intensities_struct, model_values);
	if (result != GSL_SUCCESS)
		return result;
	result = Gauss_2D_fit_function_Jacobian_FixedWidth(params, measured_intensities_struct, jacobian);
	if (result != GSL_SUCCESS)
		return result;
	
	return GSL_SUCCESS;
}

int Ellipsoidal_Gauss_2D_fit_function_and_Jacobian(const gsl_vector *params, void *measured_intensities_struct, gsl_vector *model_values, gsl_matrix *jacobian) {
	int result;
	result = Ellipsoidal_Gauss_2D_fit_function(params, measured_intensities_struct, model_values);
	if (result != GSL_SUCCESS)
		return result;
	result = Ellipsoidal_Gauss_2D_fit_function_Jacobian(params, measured_intensities_struct, jacobian);
	if (result != GSL_SUCCESS)
		return result;
	
	return GSL_SUCCESS;
}

boost::shared_ptr<LocalizedPositionsContainer> FitPositions::fit_positions(const boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<std::vector<position> > positions) {
	size_t startPosition, endPosition;
	
	startPosition = 0;
	endPosition = positions->size() - 1;
	
	return fit_positions(image, positions, startPosition, endPosition);
}

boost::shared_ptr<LocalizedPositionsContainer> FitPositionsGaussian::fit_positions(const boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<std::vector<position> > positions, 
																					   size_t startPos, size_t endPos) {
	
	// some safety checks
	if (positions->size() == 0) {
		// if no positions were found then there is no reason to run the analysis
		// we need to catch this here and not upstream since then we can return an appropriate
		// instance of LocalizedPositionsContainer
		return boost::shared_ptr<LocalizedPositionsContainer_2DGauss> (new LocalizedPositionsContainer_2DGauss());
	}
	
	if ((endPos >= positions->size()) || (startPos >= positions->size())) {
		std::string error;
		error = "Requested start and end positions are outside the range of positions supplied in FitPositionsGaussian::fit_positions";
		throw std::range_error(error);
	}
	
	if (startPos > endPos) {
		std::string error;
		error = "Start is beyond end in FitPositionsGaussian::fit_positions";
		throw std::range_error(error);
	}
	
	size_t size_of_subset = 2 * cutoff_radius + 1;
	size_t x_offset, y_offset, x_max, y_max;
	size_t number_of_intensities = size_of_subset * size_of_subset;
	size_t xSize = image->getXSize();
	size_t ySize = image->getYSize();
	
	double x0_initial, y0_initial, amplitude, background;
	double chi, degreesOfFreedom, c;
	long iterations = 0;
	int status;
	
	double relativeAmplitudeError, relativeWidthError;
	
	boost::shared_ptr<PALMMatrix<double> > image_subset;
	boost::shared_ptr<LocalizedPositionsContainer_2DGauss> fitted_positions (new LocalizedPositionsContainer_2DGauss());
	boost::shared_ptr<LocalizedPosition_2DGauss> localizationResult (new LocalizedPosition_2DGauss());
	
	image_subset = boost::shared_ptr<PALMMatrix<double> > (new PALMMatrix<double>(size_of_subset, size_of_subset));
	
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
	
	f.f = &Gauss_2D_fit_function;
	f.df = &Gauss_2D_fit_function_Jacobian;
	f.fdf = &Gauss_2D_fit_function_and_Jacobian;
	f.n = number_of_intensities;
	f.p = 5;
	f.params = (void *)&fitData;
	
	
	// iterate over all the determined positions
	for (size_t i = startPos; i <= endPos; i++) {
		iterations = 0;
		
		amplitude = (*positions)[i].get_intensity();
		x0_initial = (*positions)[i].get_x();
		y0_initial = (*positions)[i].get_y();
		background = (*positions)[i].get_background();
		
		x_offset = x0_initial - cutoff_radius;
		y_offset = y0_initial - cutoff_radius;
		x_max = x0_initial + cutoff_radius;
		y_max = y0_initial + cutoff_radius;
		
		if ((x_offset < 0) || (x_max > (xSize - 1)) || (y_offset < 0) || (y_max > (ySize - 1))) {	// this position is too close to the edge of the image
			// we cannot include it
			continue;
		}
		
		for (size_t j = x_offset; j <= x_max; j++) {
			for (size_t k = y_offset; k <= y_max; k++) {
				image_subset->set(j - x_offset, k - y_offset, image->get(j, k));
			}
		}
		
		fitData.xOffset = (double)x_offset;
		fitData.yOffset = (double)y_offset;
		fitData.imageSubset = image_subset;
		fitData.sigma = sigma;
		
		// provide the initial parameters
		gsl_vector_set(fit_parameters, 0, amplitude);
		gsl_vector_set(fit_parameters, 1, r_initial);
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
			continue;
		}
		
		// are the reported positions within the window?
		if ((gsl_vector_get(fit_iterator->x, 2) < x_offset) || (gsl_vector_get(fit_iterator->x, 2) > x_max) || (gsl_vector_get(fit_iterator->x, 3) < y_offset) || (gsl_vector_get(fit_iterator->x, 3) > y_max)) {
			// the reported positions are not within the window, we should reject them
			continue;
		}
		
		// are the fit results close enough to the initial values to be trusted?
		if ((gsl_vector_get(fit_iterator->x, 0) < amplitude / 2.0) || (gsl_vector_get(fit_iterator->x, 0) > amplitude * 1.5)) {
			// the output fit amplitude is more than a factor of two different from the initial value, drop this point
			continue;
		}
		
		if ((gsl_vector_get(fit_iterator->x, 1) < r_initial / 2.0) || (gsl_vector_get(fit_iterator->x, 1) > r_initial * 1.5)) {
			// the output fit width is more than a factor of two different from the initial value, drop this point
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
	}
	
	gsl_multifit_fdfsolver_free(fit_iterator);
	gsl_vector_free(fit_parameters);
	gsl_matrix_free(covarianceMatrix);
	
	return fitted_positions;
	
}

boost::shared_ptr<LocalizedPositionsContainer> FitPositionsGaussian_FixedWidth::fit_positions(const boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<std::vector<position> > positions, 
																								  size_t startPos, size_t endPos) {
	
	// some safety checks
	if (positions->size() == 0) {
		// if no positions were found then there is no reason to run the analysis
		// we need to catch this here and not upstream since then we can return an appropriate
		// instance of LocalizedPositionsContainer
		return boost::shared_ptr<LocalizedPositionsContainer_2DGaussFixedWidth> (new LocalizedPositionsContainer_2DGaussFixedWidth());
	}
	
	if ((endPos >= positions->size()) || (startPos >= positions->size())) {
		std::string error;
		error = "Requested start and end positions are outside the range of positions supplied in FitPositionsGaussian::fit_positions";
		throw std::range_error(error);
	}
	
	if (startPos > endPos) {
		std::string error;
		error = "Start is beyond end in FitPositionsGaussian::fit_positions";
		throw std::range_error(error);
	}
	
	size_t size_of_subset = 2 * cutoff_radius + 1;
	size_t x_offset, y_offset, x_max, y_max;
	size_t number_of_intensities = size_of_subset * size_of_subset;
	size_t xSize = image->getXSize();
	size_t ySize = image->getYSize();
	
	double x0_initial, y0_initial, amplitude, background;
	double chi, degreesOfFreedom, c;
	long iterations = 0;
	int status;
	
	double relativeAmplitudeError;
	
	boost::shared_ptr<PALMMatrix<double> > image_subset;
	boost::shared_ptr<LocalizedPositionsContainer_2DGaussFixedWidth> fitted_positions (new LocalizedPositionsContainer_2DGaussFixedWidth());
	boost::shared_ptr<LocalizedPosition_2DGaussFixedWidth> localizationResult (new LocalizedPosition_2DGaussFixedWidth());
	
	image_subset = boost::shared_ptr<PALMMatrix<double> > (new PALMMatrix<double>(size_of_subset, size_of_subset));
	
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
	
	f.f = &Gauss_2D_fit_function_FixedWidth;
	f.df = &Gauss_2D_fit_function_Jacobian_FixedWidth;
	f.fdf = &Gauss_2D_fit_function_and_Jacobian_FixedWidth;
	f.n = number_of_intensities;
	f.p = 4;
	f.params = (void *)&fitData;
	
	
	// iterate over all the determined positions
	for (size_t i = startPos; i <= endPos; i++) {
		iterations = 0;
		
		amplitude = (*positions)[i].get_intensity();
		x0_initial = (*positions)[i].get_x();
		y0_initial = (*positions)[i].get_y();
		background = (*positions)[i].get_background();
		
		x_offset = x0_initial - cutoff_radius;
		y_offset = y0_initial - cutoff_radius;
		x_max = x0_initial + cutoff_radius;
		y_max = y0_initial + cutoff_radius;
		
		if ((x_offset < 0) || (x_max > (xSize - 1)) || (y_offset < 0) || (y_max > (ySize - 1))) {	// this position is too close to the edge of the image
			// we cannot include it
			continue;
		}
		
		for (size_t j = x_offset; j <= x_max; j++) {
			for (size_t k = y_offset; k <= y_max; k++) {
				image_subset->set(j - x_offset, k - y_offset, image->get(j, k));
			}
		}
		
		fitData.xOffset = (double)x_offset;
		fitData.yOffset = (double)y_offset;
		fitData.imageSubset = image_subset;
		fitData.sigma = sigma;
		fitData.width = r_initial;
		
		// provide the initial parameters
		gsl_vector_set(fit_parameters, 0, amplitude);
		// gsl_vector_set(fit_parameters, 1, r_initial * 1.414213562373095);	// because the fitting function is of the form 1/r^2, but standard deviation is 1/(2 r^2), we have to correct by sqrt(2)
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
			continue;
		}
		
		// are the reported positions within the window?
		if ((gsl_vector_get(fit_iterator->x, 1) < x_offset) || (gsl_vector_get(fit_iterator->x, 1) > x_max) || (gsl_vector_get(fit_iterator->x, 2) < y_offset) || (gsl_vector_get(fit_iterator->x, 2) > y_max)) {
			// the reported positions are not within the window, we should reject them
			continue;
		}
		
		// are the fit results close enough to the initial values to be trusted?
		if ((gsl_vector_get(fit_iterator->x, 0) < amplitude / 2.0) || (gsl_vector_get(fit_iterator->x, 0) > amplitude * 2.0)) {
			// the output fit amplitude is more than a factor of two different from the initial value, drop this point
			continue;
		}
		
		// calculate the covariance matrix
		gsl_multifit_covar(fit_iterator->J, 0.0, covarianceMatrix);
		chi = gsl_blas_dnrm2(fit_iterator->f);
		degreesOfFreedom = (2 * cutoff_radius - 1) * (2 * cutoff_radius - 1) - 4;
		c = GSL_MAX_DBL(1, chi / sqrt(degreesOfFreedom));
		
		// store the data
		// store the fitted parameters
		
		// the width returned by the fit function is not equal to the standard deviation (a factor of sqrt 2 is missing)
		// so we correct for that
		localizationResult->width = r_initial;
		localizationResult->xPosition = gsl_vector_get(fit_iterator->x, 1);
		localizationResult->yPosition = gsl_vector_get(fit_iterator->x, 2);
		localizationResult->background = gsl_vector_get(fit_iterator->x, 3);
		
		localizationResult->integral = 2 * PI * localizationResult->width * localizationResult->width * gsl_vector_get(fit_iterator->x, 0);
		
		
		localizationResult->xPositionDeviation = c * sqrt(gsl_matrix_get(covarianceMatrix, 1, 1));
		localizationResult->yPositionDeviation = c * sqrt(gsl_matrix_get(covarianceMatrix, 2, 2));
		localizationResult->backgroundDeviation = c * sqrt(gsl_matrix_get(covarianceMatrix, 3, 3));
		
		// calculate the uncertainty on the integrated intensity using error propagation
		relativeAmplitudeError = c * sqrt(gsl_matrix_get(covarianceMatrix, 0, 0)) / gsl_vector_get(fit_iterator->x, 0); // the relative error on the amplitude
		localizationResult->integralDeviation = 2 * PI * localizationResult->integral * localizationResult->width * localizationResult->width * sqrt(relativeAmplitudeError * relativeAmplitudeError);
		
		
		fitted_positions->addPosition(localizationResult);
	}
	
	gsl_multifit_fdfsolver_free(fit_iterator);
	gsl_vector_free(fit_parameters);
	gsl_matrix_free(covarianceMatrix);
	
	return fitted_positions;
	
}


boost::shared_ptr<LocalizedPositionsContainer> FitPositionsEllipsoidalGaussian::fit_positions(const boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<std::vector<position> > positions,
															 size_t startPos, size_t endPos) {
	
	// some safety checks
	if (positions->size() == 0) {
		// if no positions were found then there is no reason to run the analysis
		// we need to catch this here and not upstream since then we can return an appropriate
		// instance of LocalizedPositionsContainer
		return boost::shared_ptr<LocalizedPositionsContainer> (new LocalizedPositionsContainer_Ellipsoidal2DGaussian());
	}
	
	if ((endPos >= positions->size()) || (startPos >= positions->size())) {
		std::string error;
		error = "Requested start and end positions are outside the range of positions supplied in FitPositionsGaussian::fit_positions";
		throw std::range_error(error);
	}
	
	if (startPos > endPos) {
		std::string error;
		error = "Start is beyond end in FitPositionsGaussian::fit_positions";
		throw std::range_error(error);
	}
	
	size_t size_of_subset = 2 * cutoff_radius + 1;
	size_t x_offset, y_offset, x_max, y_max;
	size_t number_of_intensities = size_of_subset * size_of_subset;
	size_t xSize = image->getXSize();
	size_t ySize = image->getYSize();
	
	double x0_initial, y0_initial, amplitude, background, correlation_initial;
	double chi, degreesOfFreedom, c;
	long iterations = 0;
	int status;
	
	boost::shared_ptr<PALMMatrix<double> > image_subset;
	boost::shared_ptr<LocalizedPositionsContainer_Ellipsoidal2DGaussian> fitted_positions (new LocalizedPositionsContainer_Ellipsoidal2DGaussian());
	boost::shared_ptr<LocalizedPosition_Ellipsoidal2DGauss> localizationResult (new LocalizedPosition_Ellipsoidal2DGauss());
	
	image_subset = boost::shared_ptr<PALMMatrix<double> > (new PALMMatrix<double>(size_of_subset, size_of_subset));
	
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
	
	f.f = &Ellipsoidal_Gauss_2D_fit_function;
	f.df = &Ellipsoidal_Gauss_2D_fit_function_Jacobian;
	f.fdf = &Ellipsoidal_Gauss_2D_fit_function_and_Jacobian;
	f.n = number_of_intensities;
	f.p = 7;
	f.params = (void *)&fitData;
	
	
	// iterate over all the determined positions
	for (size_t i = startPos; i <= endPos; i++) {
		iterations = 0;
		
		amplitude = (*positions)[i].get_intensity();
		x0_initial = (*positions)[i].get_x();
		y0_initial = (*positions)[i].get_y();
		correlation_initial = 0.0;
		background = (*positions)[i].get_background();
		
		x_offset = x0_initial - cutoff_radius;
		y_offset = y0_initial - cutoff_radius;
		x_max = x0_initial + cutoff_radius;
		y_max = y0_initial + cutoff_radius;
		
		if ((x_offset < 0) || (x_max > (xSize - 1)) || (y_offset < 0) || (y_max > (ySize - 1))) {	// this position is too close to the edge of the image
			// we cannot include it
			continue;
		}
		
		for (size_t j = x_offset; j <= x_max; j++) {
			for (size_t k = y_offset; k <= y_max; k++) {
				image_subset->set(j - x_offset, k - y_offset, image->get(j, k));
			}
		}
		
		fitData.xOffset = (double)x_offset;
		fitData.yOffset = (double)y_offset;
		fitData.imageSubset = image_subset;
		fitData.sigma = sigma;
		
		// provide the initial parameters
		gsl_vector_set(fit_parameters, 0, amplitude);
		gsl_vector_set(fit_parameters, 1, this->r_initial);
		gsl_vector_set(fit_parameters, 2, this->r_initial);
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
			status = gsl_multifit_test_delta(fit_iterator->dx, fit_iterator->x, 1e-4, 1e-4);
		} while ((status == GSL_CONTINUE) && (iterations < 200));
		
		if (gsl_vector_get(fit_iterator->x, 0) <= 0) {	// reject fits that have negative amplitudes
			continue;
		}
		
		// are the reported positions within the window?
		if ((gsl_vector_get(fit_iterator->x, 3) < x_offset) || (gsl_vector_get(fit_iterator->x, 3) > x_max) || (gsl_vector_get(fit_iterator->x, 4) < y_offset) || (gsl_vector_get(fit_iterator->x, 4) > y_max)) {
			// the reported positions are not within the window, we should reject them
			continue;
		}
		
		// are the fit results close enough to the initial values to be trusted?
		if ((gsl_vector_get(fit_iterator->x, 0) < amplitude / 2.0) || (gsl_vector_get(fit_iterator->x, 0) > amplitude * 1.5)) {
			// the output fit amplitude is more than a factor of two different from the initial value, drop this point
			continue;
		}
		
		// are the standard deviations in the x and y direction close enough?
		if ((0.75 * gsl_vector_get(fit_iterator->x, 1) > gsl_vector_get(fit_iterator->x, 2)) || (1.25 * gsl_vector_get(fit_iterator->x, 1) < gsl_vector_get(fit_iterator->x, 2))) {
			// the standard deviations in x and y differ by more than 25%, drop this point
			continue;
		}
		
		if ((gsl_vector_get(fit_iterator->x, 1) < r_initial / 2.0) || (gsl_vector_get(fit_iterator->x, 1) > r_initial * 1.5)) {
			// the output fit width is more than a factor of two different from the initial value, drop this point
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
	}
	
	gsl_multifit_fdfsolver_free(fit_iterator);
	gsl_vector_free(fit_parameters);
	gsl_matrix_free(covarianceMatrix);
	
	return fitted_positions;
	
}


boost::shared_ptr<LocalizedPositionsContainer> FitPositionsMultiplication::fit_positions(const boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<std::vector<position> > positions, 
																							 size_t startPos, size_t endPos) {
	
	// some safety checks
	if (positions->size() == 0) {
		// if no positions were found then there is no reason to run the analysis
		// we need to catch this here and not upstream since then we can return an appropriate
		// instance of LocalizedPositionsContainer
		return boost::shared_ptr<LocalizedPositionsContainer_Multiplication> (new LocalizedPositionsContainer_Multiplication());
	}
	
	if ((endPos >= positions->size()) || (startPos >= positions->size())) {
		std::string error;
		error = "Requested start and end positions are outside the range of positions supplied in FitPositionsGaussian::fit_positions";
		throw std::range_error(error);
	}
	
	if (startPos > endPos) {
		std::string error;
		error = "Start is beyond end in FitPositionsGaussian::fit_positions";
		throw std::range_error(error);
	}
	
	size_t size_of_subset = 2 * cutoff_radius + 1;
	size_t xSize = image->getXSize();
	size_t ySize = image->getYSize();
	size_t x_offset, y_offset, x_max, y_max;
	
	double x0_initial, y0_initial, amplitude, background;
	size_t iterations = 0;
	int converged;
	
	double convergence_treshold_squared = convergence_threshold * convergence_threshold;
	double delta_squared = 10 * convergence_treshold_squared;	// this test the convergence of the position determined by the iteration
	// it is the distance between (xn-1, yn-1) and (xn, yn)
	// we initialize it to a value well over the treshold so that we will run at least two iterations
	double previous_position_x, previous_position_y;
	double current_x, current_y;
	
	boost::shared_ptr<PALMMatrix<double> > image_subset;
	boost::shared_ptr<PALMMatrix<double> > image_subset_mask;
	boost::shared_ptr<LocalizedPositionsContainer_Multiplication> fitted_positions (new LocalizedPositionsContainer_Multiplication());
	boost::shared_ptr<LocalizedPosition_Multiplication> localizationResult (new LocalizedPosition_Multiplication());
	
	image_subset = boost::shared_ptr<PALMMatrix<double> >(new PALMMatrix<double>(size_of_subset, size_of_subset));
	image_subset_mask = boost::shared_ptr<PALMMatrix<double> > (new PALMMatrix<double>(size_of_subset, size_of_subset));
	
	for (size_t i = startPos; i <= endPos; ++i) {
		amplitude = (*positions)[i].get_intensity();
		x0_initial = (*positions)[i].get_x();
		y0_initial = (*positions)[i].get_y();
		background = (*positions)[i].get_background();
		
		x_offset = (size_t)x0_initial - cutoff_radius;
		y_offset = (size_t)y0_initial - cutoff_radius;
		x_max = (size_t)x0_initial + cutoff_radius;
		y_max = (size_t)y0_initial + cutoff_radius;
		
		if ((x_offset < 0) || (x_max > (xSize - 1)) || (y_offset < 0) || (y_max > (ySize - 1))) {	// this position is too close to the edge of the image, we cannot include it
			continue;
		}
		
		for (size_t k = y_offset; k <= y_max; ++k) {
			for (size_t j = x_offset; j <= x_max; ++j) {
				image_subset->set(j - x_offset, k - y_offset, image->get(j, k));
			}
		}
		
		iterations = 0;
		
		current_x = x0_initial - (double)x_offset;	// correct the x- and y-values for the fact that we analyze in a subset of the image rather than the complete frame
		current_y = y0_initial - (double)y_offset;
		
		converged = 1;
		while (delta_squared > convergence_treshold_squared) {
			previous_position_x = current_x;
			previous_position_y = current_y;
			
			++iterations;
			
			if (iterations > 100) {	// the multiplication is not converging, we should stop
				converged = 0;
				break;
			}
			
			multiply_with_gaussian(image_subset, image_subset_mask, current_x, current_y, r_initial, background, amplitude);
			determine_x_y_position(image_subset_mask, current_x, current_y);
			
			if (iterations == 1)	// this is the first iteration, we should not check for termination
				continue;
			
			delta_squared = (current_x - previous_position_x) * (current_x - previous_position_x) + (current_y - previous_position_y) * (current_y - previous_position_y);
		}
		
		if (converged == 0)
			continue;
		
		delta_squared = 10 * convergence_treshold_squared;
		
		localizationResult->width = r_initial;
		localizationResult->xPosition = (double)current_x + (double)x_offset;
		localizationResult->yPosition = (double)current_y + (double)y_offset;
		
		fitted_positions->addPosition(localizationResult);
	}
	
	return fitted_positions;
}


int FitPositionsMultiplication::multiply_with_gaussian(boost::shared_ptr<PALMMatrix<double> > original_image, boost::shared_ptr<PALMMatrix<double> > masked_image, double x, double y, 
													   double std_dev, double background, double amplitude) {
	// we will replace the contents of masked_image with the multiplication of original_image and a gaussian centered at position (x,y)
	
	size_t x_size = masked_image->getXSize();
	size_t y_size = masked_image->getYSize();
	
	double gaussian_value, distance_squared;
	
	if ((original_image->getXSize() != x_size) || (original_image->getYSize() != y_size)) {
		throw DIMENSIONS_SHOULD_BE_EQUAL(std::string("Matrix dimensions are not equal in FitPositionsMultiplication::multiply_with_gaussian"));
	}
	
	for (size_t j = 0; j < y_size; j++) {
		for (size_t i = 0; i < x_size; i++) {
			distance_squared = (x - (double)i) * (x - (double)i) + (y - (double)j) * (y - (double)j);
			
			gaussian_value = amplitude * exp(- distance_squared / (2 * std_dev * std_dev)) + background;
			
			masked_image->set(i, j, gaussian_value * original_image->get(i, j));
		}
	}
	
	return 0;
}


int FitPositionsMultiplication::determine_x_y_position(boost::shared_ptr<PALMMatrix<double> > masked_image, double &x, double &y) {
	// based on eq (3) in Thompson Biophys J 2002
	
	size_t x_size = (size_t)masked_image->getXSize();
	size_t y_size = (size_t)masked_image->getYSize();
	
	double numerator_x = 0, denominator = 0;
	double numerator_y = 0;
	
	// start with determining the x-position
	for (size_t j = 0; j < y_size; j++) {
		for (size_t i = 0; i < x_size; i++) {
			numerator_x += (double)i * masked_image->get(i, j);
			numerator_y += (double)j * masked_image->get(i, j);
			denominator += masked_image->get(i, j);
		}
	}
	
	x = numerator_x / denominator;
	
	y = numerator_y / denominator;
	
	return 0;
}


boost::shared_ptr<LocalizedPositionsContainer> FitPositionsCentroid::fit_positions(const boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<std::vector<position> > positions, 
																					   size_t startPos, size_t endPos) {
	
	// some safety checks
	if (positions->size() == 0) {
		// if no positions were found then there is no reason to run the analysis
		// we need to catch this here and not upstream since then we can return an appropriate
		// instance of LocalizedPositionsContainer
		return boost::shared_ptr<LocalizedPositionsContainer_Centroid> (new LocalizedPositionsContainer_Centroid());
	}
	
	if ((endPos >= positions->size()) || (startPos >= positions->size())) {
		std::string error;
		error = "Requested start and end positions are outside the range of positions supplied in FitPositionsCentroid::fit_positions";
		throw std::range_error(error);
	}
	
	if (startPos > endPos) {
		std::string error;
		error = "Start is beyond end in FitPositionsCentroid::fit_positions";
		throw std::range_error(error);
	}
	
	size_t xSize = image->getXSize();
	size_t ySize = image->getYSize();
	size_t x_offset, y_offset, x_max, y_max;
	
	size_t x0_initial, y0_initial;
	double current_x, current_y;
	double denominator;
	
	boost::shared_ptr<LocalizedPositionsContainer_Centroid> fitted_positions (new LocalizedPositionsContainer_Centroid());
	boost::shared_ptr<LocalizedPosition_Centroid> localizationResult (new LocalizedPosition_Centroid());
	
	for (size_t i = startPos; i <= endPos; ++i) {
		x0_initial = (*positions)[i].get_x();
		y0_initial = (*positions)[i].get_y();
		current_x = 0;
		current_y = 0;
		denominator = 0;
		
		x_offset = x0_initial - cutoff_radius;
		y_offset = y0_initial - cutoff_radius;
		x_max = x0_initial + cutoff_radius;
		y_max = y0_initial + cutoff_radius;
		
		if ((x_offset > xSize) || (x_max > (xSize - 1)) || (y_offset > ySize) || (y_max > (ySize - 1))) {	// the point is too close to the edge
			// because all the variables are unsigned, a negative value will
			// actually end up being larger than xSize or ySize
			continue;
		}
		
		for (size_t j = x_offset; j <= x_max; ++j) {
			for (size_t k = y_offset; k <= y_max; ++k) {
				current_x += (double)j * image->get(j, k);
				current_y += (double)k * image->get(j, k);
				denominator += image->get(j, k);
			}
		}
		
		current_x /= denominator;
		current_y /= denominator;		
		
		localizationResult->xPosition = current_x;
		localizationResult->yPosition = current_y;
		
		fitted_positions->addPosition(localizationResult);
	}
	
	return fitted_positions;
}

boost::shared_ptr<LocalizedPositionsContainer> FitPositionsDeflate::fit_positions(const boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<std::vector<position> > positions, size_t startPos, size_t endPos) {
	// TODO: for now we ignore the starting positions and ending position provided as arguments
	
	boost::shared_ptr<LocalizedPositionsContainer> positionsFittedThusFar;
	boost::shared_ptr<LocalizedPositionsContainer> positionsLocalizedThisFrame;
	boost::shared_ptr<PALMMatrix<double> > subtractedImage;
	boost::shared_ptr<PALMMatrix <unsigned char> > segmentedImage;
	boost::shared_ptr<std::vector<position> > locatedParticles;
	
	// get an initial set of positions to start from
	positionsFittedThusFar = this->positionsFitter->fit_positions(image, positions);
	
	// if no positions were localized then don't bother
	if (positionsFittedThusFar->getNPositions() == 0)
		return positionsFittedThusFar;
	
	// check if the type of FitPositions provided returns the required information to reconstruct
	// a Gaussian emitter
	if ((positionsFittedThusFar->getIntegral(0) == 0) || (positionsFittedThusFar->getXWidth(0) == 0) || (positionsFittedThusFar->getYWidth(0) == 0) || (positionsFittedThusFar->getBackground(0) == 0))
		throw std::runtime_error("The selected fitting algorithm does not provide sufficient information for deflation analysis");
	
	// set up for iteration
	positionsLocalizedThisFrame = positionsFittedThusFar;
	subtractedImage = image;
	
	while (positionsLocalizedThisFrame->getNPositions() != 0) {
		// continue while new positions are still being localized
		subtractedImage = this->subtractLocalizedPositions(subtractedImage, positionsLocalizedThisFrame);
		
		segmentedImage = do_processing_and_thresholding(subtractedImage, this->preprocessor, 
																					  this->thresholder, this->postprocessor);
		locatedParticles = this->particleFinder->findPositions(subtractedImage, segmentedImage);
		positionsLocalizedThisFrame = this->positionsFitter->fit_positions(subtractedImage, locatedParticles);
		
		// if we found new particles then append them
		if (positionsLocalizedThisFrame->getNPositions() != 0) {
			positionsFittedThusFar->addPositions(positionsLocalizedThisFrame);
		}
	}
	
	return positionsFittedThusFar;
	
}


boost::shared_ptr<PALMMatrix<double> > FitPositionsDeflate::subtractLocalizedPositions(boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<LocalizedPositionsContainer> positions) {
	double fittedXPos, fittedYPos, fittedIntegral, fittedXWidth, fittedYWidth;
	double centerX, centerY, calculatedAmplitude;
	size_t startX, endX, startY, endY;
	double distanceXSquared, distanceYSquared, currentIntensity;
	
	size_t nPositions = positions->getNPositions();
	size_t xSize = image->getXSize();
	size_t ySize = image->getYSize();
	boost::shared_ptr<PALMMatrix<double> > outputImage(new PALMMatrix<double> (xSize, ySize));
	outputImage = image;
	
	for (size_t n = 0; n < nPositions; ++n) {
		fittedIntegral = positions->getIntegral(n);
		fittedXWidth = positions->getXWidth(n);
		fittedYWidth = positions->getYWidth(n);
		fittedXPos = positions->getXPosition(n);
		fittedYPos = positions->getYPosition(n);
		
		calculatedAmplitude = fittedIntegral / (2 * PI * fittedXWidth * fittedYWidth);
		
		centerX = fittedXPos;
		centerY = fittedYPos;
		
		startX = floor((double)centerX - 4.0 * fittedXWidth);	// only run the calculation over a subset of the image surrounding the position
		startY = floor((double)centerY - 4.0 * fittedYWidth);
		endX = ceil((double)centerX + 4.0 * fittedXWidth);
		endY = ceil((double)centerY + 4.0 * fittedYWidth);
		
		if (startX < 0)
			startX = 0;
		if (endX >= xSize)
			endX = xSize - 1;
		if (startY < 0)
			startY = 0;
		if (endY >= ySize)
			endY = ySize - 1;
		
		for (size_t i = startX; i <= endX; ++i) {
			for (size_t j = startY; j < endY; ++j) {
				distanceXSquared = ((double)i - (double)centerX) * ((double)i - (double)centerX);
				distanceYSquared = ((double)j - (double)centerY) * ((double)j - (double)centerY);
				currentIntensity = calculatedAmplitude * exp(- (distanceXSquared + distanceYSquared) / (2 * fittedXWidth * fittedYWidth));
				
				(*outputImage)(i, j) = (*outputImage)(i, j) - currentIntensity;
			}
		}
	}
	return outputImage;
}


