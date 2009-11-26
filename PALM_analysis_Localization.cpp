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
	size_t ySize = imageSubset->getXSize();
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
			function_value = offset + amplitude * exp(- (((x0 - x)/ r) * ((x0 - x)/ r) + ((y0 - y) / r) * ((y0 - y) / r)));
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
	size_t ySize = imageSubset->getXSize();
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
			function_value = offset + amplitude * exp(- (((x0 - x)/ r) * ((x0 - x)/ r) + ((y0 - y) / r) * ((y0 - y) / r)));
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
			
			exp_factor = exp(- ((x0 - x)/ r) * ((x0 - x)/ r) - ((y0 - y) / r) * ((y0 - y) / r));
			
			dfdA = exp_factor / sigma;
			dfdr = (2 * (y - y0) * (y - y0) / r / r / r + 2 * (x - x0) * (x - x0) / r / r / r) * exp_factor * amplitude / sigma;
			dfdx0 = (2 * (x - x0) * exp_factor * amplitude) / (r * r *sigma);
			dfdy0 = (2 * (y - y0) * exp_factor * amplitude) / (r * r *sigma);
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
	double dfdA, dfdr, dfdx0, dfdy0,dfdoffset;
	
	if (r == 0) {
		return GSL_FAILURE;
	}
	
	for (size_t i = 0; i < xSize; ++i) {
		for (size_t j = 0; j < ySize; ++j) {
			x = xOffset + (double)i;
			y = yOffset + (double)j;
			
			exp_factor = exp(- ((x0 - x)/ r) * ((x0 - x)/ r) - ((y0 - y) / r) * ((y0 - y) / r));
			
			dfdA = exp_factor / sigma;
			// dfdr = (2 * (y - y0) * (y - y0) / r / r / r + 2 * (x - x0) * (x - x0) / r / r / r) * exp_factor * amplitude / sigma;
			dfdx0 = (2 * (x - x0) * exp_factor * amplitude) / (r * r *sigma);
			dfdy0 = (2 * (y - y0) * exp_factor * amplitude) / (r * r *sigma);
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

int Gauss_2D_fit_function_and_Jacobian(const gsl_vector *params, void *measured_intensities_struct, gsl_vector *model_values, gsl_matrix *jacobian) {
	
	Gauss_2D_fit_function(params, measured_intensities_struct, model_values);
	Gauss_2D_fit_function_Jacobian(params, measured_intensities_struct, jacobian);
	
	return GSL_SUCCESS;
}

int Gauss_2D_fit_function_and_Jacobian_FixedWidth(const gsl_vector *params, void *measured_intensities_struct, gsl_vector *model_values, gsl_matrix *jacobian) {
	
	Gauss_2D_fit_function_FixedWidth(params, measured_intensities_struct, model_values);
	Gauss_2D_fit_function_Jacobian_FixedWidth(params, measured_intensities_struct, jacobian);
	
	return GSL_SUCCESS;
}

/*	Maxima code to generate the required derivatives for a 2D ellipsoidal Gaussian
	g: A * %e^(-(a * (x - x0)**2 + 2*b*(x - x0)*(y-y0) + c * (y - y0)**2)) + offset;
	\begin{align}
		g &= e^{-c\,\left(y-{\it y_0}\right)^2-2\,b\,\left(x-{\it x_0}\right)\,\left(y-{\it y_0}\right)-a\,\left(x-{\it x_0}\right)^2}\,A+z_0 \\
		\frac{\partial g}{\partial A} &= e^{-c\,\left(y-{\it y_0}\right)^2-2\,b\,\left(x-{\it x_0}\right)\,\left(y-{\it y_0}\right)-a\,\left(x-{\it x_0}\right)^2} \\
		\frac{\partial g}{\partial x_0} &= e^{-c\,\left(y-{\it y_0}\right)^2-2\,b\,\left(x-{\it x_0}\right)\,\left(y-{\it y_0}\right)-a\,\left(x-{\it x_0}\right)^2}\, \left(2\,b\,\left(y-{\it y_0}\right)+2\,a\,\left(x-{\it x_0}\right) \right)\,A \\
		\frac{\partial g}{\partial y_0} &= e^{-c\,\left(y-{\it y_0}\right)^2-2\,b\,\left(x-{\it x_0}\right)\,\left(y-{\it y_0}\right)-a\,\left(x-{\it x_0}\right)^2}\,\left(2\,c\,\left(y-{\it y_0}\right)+2\,b\,\left(x-{\it x_0}\right)\right)\,A \\
		\frac{\partial g}{\partial a} &= -e^{-c\,\left(y-{\it y_0}\right)^2-2\,b\,\left(x-{\it x_0}\right)\,\left(y-{\it y_0}\right)-a\,\left(x-{\it x_0}\right)^2}\,\left(x-{\it x_0}\right)^2\,A \\
		\frac{\partial g}{\partial b} &= e^{-c\,\left(y-{\it y_0}\right)^2-2\,b\,\left(x-{\it x_0}\right)\,\left(y-{\it y_0}\right)-a\,\left(x-{\it x_0}\right)^2}\,\left(x-{\it x_0}\right)\,\left(y-{\it y_0}\right)\,A \\
		\frac{\partial g}{\partial c} &= -e^{-c\,\left(y-{\it y_0}\right)^2-2\,b\,\left(x-{\it x_0}\right)\,\left(y-{\it y_0}\right)-a\,\left(x-{\it x_0}\right)^2}\,\left(y-{\it y_0}\right)^2\,A \\
		\frac{\partial g}{\partial z_0} &= 1
 \end{align}
 */


/*int Gauss_2D_Poissonian_fit_function(const gsl_vector *params, void *measured_intensities_struct, gsl_vector *model_values) {
 measured_data_Gauss_fits *intensities_local = (measured_data_Gauss_fits *)measured_intensities_struct;
 
 size_t number_of_intensities = intensities_local->get_number_of_intensities();
 double *measured_intensities = intensities_local->get_intensities();
 // double *sigma = intensities_local->get_sigma();
 size_t x_size = intensities_local->getXSize();
 size_t y_size = intensities_local->getYSize();
 double x_offset = intensities_local->get_x_offset();
 double y_offset = intensities_local->get_y_offset();
 
 double amplitude = gsl_vector_get(params, 0);
 double r = gsl_vector_get(params, 1);
 double x0 = gsl_vector_get(params, 2);
 double y0 = gsl_vector_get(params, 3);
 double offset = gsl_vector_get(params, 4);
 
 double x,y, function_value, square_deviation;
 
 if (r == 0) {
 return GSL_FAILURE;
 }
 
 for (size_t i = 0; i < number_of_intensities; i++) {
 return_xy_from_array(i, x_size, y_size, x_offset, y_offset, x, y);
 
 function_value = offset + amplitude * exp(- (((x0 - x)/ r) * ((x0 - x)/ r) + ((y0 - y) / r) * ((y0 - y) / r)));
 // square_deviation = (function_value - measured_intensities[i]) / sigma[i];
 square_deviation = sqrt(2 * measured_intensities[i] * log(measured_intensities[i] / function_value));
 gsl_vector_set(model_values, i, square_deviation);
 }
 return GSL_SUCCESS;
 }
 
 int Gauss_2D_Poissonian_fit_function_Jacobian(const gsl_vector *params, void *measured_intensities, gsl_matrix *jacobian) {
 
 measured_data_Gauss_fits *intensities_local = (measured_data_Gauss_fits *)measured_intensities;
 
 size_t number_of_intensities = intensities_local->get_number_of_intensities();
 // double *measured_intensities = intensities_local->get_intensities();
 // double *sigma = intensities_local->get_sigma();
 size_t x_size = intensities_local->getXSize();
 size_t y_size = intensities_local->getYSize();
 double x_offset = intensities_local->get_x_offset();
 double y_offset = intensities_local->get_y_offset();
 double *measured_intensities_array = intensities_local->get_intensities();
 
 double amplitude = gsl_vector_get(params, 0);
 double r = gsl_vector_get(params, 1);
 double x0 = gsl_vector_get(params, 2);
 double y0 = gsl_vector_get(params, 3);
 double offset = gsl_vector_get(params, 4);
 
 double x,y, exp_factor;
 double dfdA, dfdr, dfdx0, dfdy0,dfdoffset;
 
 if (r == 0) {
 return GSL_FAILURE;
 }
 
 // the maxima code to get the expression to derive from:
 // sqrt(2 * yi * log(yi / (offset + A * exp(-(((x0 - x) / r)^2 + ((y0 - y) / r)^2)))))
 
 double sqrt_2 = 1.414213562373095;
 double measured_intensity;
 double denominator;
 
 for (size_t i = 0; i < number_of_intensities; i++) {
 measured_intensity = measured_intensities_array[i];
 return_xy_from_array(i, x_size, y_size, x_offset, y_offset, x, y);
 
 exp_factor = exp(- ((x0 - x)/ r) * ((x0 - x)/ r) - ((y0 - y) / r) * ((y0 - y) / r));
 denominator = (exp_factor * amplitude  + offset) * sqrt(measured_intensity * log(measured_intensity / (exp_factor * amplitude + offset)));
 
 dfdA = - (sqrt_2 * exp_factor * measured_intensity) / (2.0 * denominator);
 dfdr = - sqrt_2 * (2.0 * (y - y0) * (y - y0) / r / r / r + 2 * (x - x0) * (x - x0) / r / r / r) * exp_factor * amplitude * measured_intensity / (2 * denominator);
 dfdx0 = (sqrt_2 * (x0 - x) * exp_factor * amplitude * measured_intensity) / (r * r * denominator);
 dfdy0 = (sqrt_2 * (y0 - y) * exp_factor * amplitude * measured_intensity) / (r * r * denominator);
 dfdoffset = (- sqrt_2 * measured_intensity) / (2.0 * denominator);
 
 gsl_matrix_set(jacobian, i, 0, dfdA);
 gsl_matrix_set(jacobian, i, 1, dfdr);
 gsl_matrix_set(jacobian, i, 2, dfdx0);
 gsl_matrix_set(jacobian, i, 3, dfdy0);
 gsl_matrix_set(jacobian, i, 4, dfdoffset);
 }
 
 return GSL_SUCCESS;
 
 }
 
 int Gauss_2D_Poissonian_fit_function_and_Jacobian(const gsl_vector *params, void *measured_intensities_struct, gsl_vector *model_values, gsl_matrix *jacobian) {
 
 Gauss_2D_Poissonian_fit_function(params, measured_intensities_struct, model_values);
 Gauss_2D_Poissonian_fit_function_Jacobian(params, measured_intensities_struct, jacobian);
 
 return GSL_SUCCESS;
 } */

boost::shared_ptr<std::vector<LocalizedPosition> > FitPositions::fit_positions(const boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<PALMMatrix<double> > positions) {
	size_t startPosition, endPosition;
	boost::shared_ptr<std::vector<LocalizedPosition> > fittedPositions;
	
	// if no positions were found, return a NULL fitted positions
	if (NULL == positions.get()) {
		fittedPositions = boost::shared_ptr<std::vector<LocalizedPosition> > (new std::vector<LocalizedPosition>);
		return fittedPositions;
	}
	
	startPosition = 0;
	endPosition = positions->getXSize() - 1;
	
	fittedPositions = fit_positions(image, positions, startPosition, endPosition);
	
	return fittedPositions;
}

boost::shared_ptr<std::vector<LocalizedPosition> > FitPositionsGaussian::fit_positions(const boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<PALMMatrix<double> > positions, 
																		   size_t startPos, size_t endPos) {
	
	// some safety checks
	if ((endPos >= positions->getXSize()) || (startPos >= positions->getXSize())) {
		string error;
		error = "Requested start and end positions are outside the range of positions supplied in FitPositionsGaussian::fit_positions";
		throw std::range_error(error);
	}
	
	if (startPos > endPos) {
		string error;
		error = "Start is beyond end in FitPositionsGaussian::fit_positions";
		throw std::range_error(error);
	}
	
	size_t number_of_positions = endPos - startPos + 1;
	size_t size_of_subset = 2 * cutoff_radius + 1;
	size_t x_offset, y_offset, x_max, y_max;
	size_t number_of_intensities = size_of_subset * size_of_subset;
	size_t xSize = image->getXSize();
	size_t ySize = image->getYSize();
	
	double x0_initial, y0_initial, amplitude, background;
	double chi, degreesOfFreedom, c;
	long iterations = 0;
	int status;
	
	boost::shared_ptr<PALMMatrix<double> > image_subset;
	boost::shared_ptr<std::vector<LocalizedPosition> > fitted_positions (new std::vector<LocalizedPosition>);
	LocalizedPosition localizationResult;
	
	image_subset = boost::shared_ptr<PALMMatrix<double> > (new PALMMatrix<double>(size_of_subset, size_of_subset));
	
	fitted_positions->reserve(number_of_positions);
	
	// initialize the solver
	const gsl_multifit_fdfsolver_type *solver;
	gsl_multifit_fdfsolver *fit_iterator;
	
	solver = gsl_multifit_fdfsolver_lmsder;
	measured_data_Gauss_fits fitData;
	gsl_multifit_function_fdf f;
	
	gsl_vector *fit_parameters = gsl_vector_alloc(5);
	if (fit_parameters == NULL) {
		string error;
		error = "unable to allocate fit_parameters in FitPositionsGaussian::fit_positions()\r";
		throw OUT_OF_MEMORY(error);
	}
	
	fit_iterator = gsl_multifit_fdfsolver_alloc(solver, number_of_intensities, 5);
	if (fit_iterator == NULL) {
		gsl_vector_free(fit_parameters);
		string error;
		error = "unable to allocate fit_iterator in FitPositionsGaussian::fit_positions()\r";
		throw OUT_OF_MEMORY(error);
	}
	
	gsl_matrix *covarianceMatrix = gsl_matrix_alloc(5, 5);
	if (covarianceMatrix == NULL) {
		gsl_vector_free(fit_parameters);
		gsl_multifit_fdfsolver_free(fit_iterator);
		string error;
		error = "unable to allocate covarianceMatrix in FitPositionsGaussian::fit_positions()\r";
		throw OUT_OF_MEMORY(error);
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
		
		amplitude = positions->get(i, 0);
		x0_initial = positions->get(i, 1);
		y0_initial = positions->get(i, 2);
		background = positions->get(i, 3);
		
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
		gsl_vector_set(fit_parameters, 1, r_initial * SQRT2);	// because the fitting function is of the form 1/r^2, but standard deviation is 1/(2 r^2), we have to correct by sqrt(2)
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
			
			status = gsl_multifit_test_delta(fit_iterator->dx, fit_iterator->x, 10, 10);
		} while ((status = GSL_CONTINUE) && (iterations < 200));
		
		if (gsl_vector_get(fit_iterator->x, 0) <= 0) {	// reject fits that have negative amplitudes
			continue;
		}
		
		// are the reported positions within the window?
		if ((gsl_vector_get(fit_iterator->x, 2) < x_offset) || (gsl_vector_get(fit_iterator->x, 2) > x_max) || (gsl_vector_get(fit_iterator->x, 3) < y_offset) || (gsl_vector_get(fit_iterator->x, 3) > y_max)) {
			// the reported positions are not within the window, we should reject them
			continue;
		}
		
		// are the fit results close enough to the initial values to be trusted?
		if ((gsl_vector_get(fit_iterator->x, 0) < amplitude / 2.0) || (gsl_vector_get(fit_iterator->x, 0) > amplitude * 2.0)) {
			// the output fit amplitude is more than a factor of two different from the initial value, drop this point
			continue;
		}
		
		if ((gsl_vector_get(fit_iterator->x, 1) < r_initial * SQRT2 / 2.0) || (gsl_vector_get(fit_iterator->x, 1) > r_initial * SQRT2 * 2.0)) {
			// the output fit width is more than a factor of two different from the initial value, drop this point
			continue;
		}
		
		// calculate the covariance matrix
		gsl_multifit_covar(fit_iterator->J, 0.0, covarianceMatrix);
		chi = gsl_blas_dnrm2(fit_iterator->f);
		degreesOfFreedom = (2 * cutoff_radius - 1) * (2 * cutoff_radius - 1) - 5;
		c = GSL_MAX_DBL(1, chi / sqrt(degreesOfFreedom));
		
		// store the fitted parameters
		localizationResult.amplitude = gsl_vector_get(fit_iterator->x, 0);
		localizationResult.width = gsl_vector_get(fit_iterator->x, 1);
		localizationResult.xPos = gsl_vector_get(fit_iterator->x, 2);
		localizationResult.yPos = gsl_vector_get(fit_iterator->x, 3);
		localizationResult.offset = gsl_vector_get(fit_iterator->x, 4);
		
		localizationResult.amplitudeError = c * sqrt(gsl_matrix_get(covarianceMatrix, 0, 0));
		localizationResult.widthError = c * sqrt(gsl_matrix_get(covarianceMatrix, 1, 1));
		localizationResult.xPosError = c * sqrt(gsl_matrix_get(covarianceMatrix, 2, 2));
		localizationResult.yPosError = c * sqrt(gsl_matrix_get(covarianceMatrix, 3, 3));
		localizationResult.offsetError = c * sqrt(gsl_matrix_get(covarianceMatrix, 4, 4));
		
		// store the number of iterations
		localizationResult.nIterations = iterations;
		
		// the width returned by the fit function is not equal to the standard deviation (a factor of sqrt 2 is missing)
		// so we correct for that
		localizationResult.width = localizationResult.width / SQRT2;
		localizationResult.widthError = localizationResult.widthError / SQRT2;	// the same for the error
		
		fitted_positions->push_back(localizationResult);
	}
	
	gsl_multifit_fdfsolver_free(fit_iterator);
	gsl_vector_free(fit_parameters);
	
	return fitted_positions;
	
}

boost::shared_ptr<std::vector<LocalizedPosition> > FitPositionsGaussian_FixedWidth::fit_positions(const boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<PALMMatrix<double> > positions, 
																					  size_t startPos, size_t endPos) {
	
	// some safety checks
	if ((endPos >= positions->getXSize()) || (startPos >= positions->getXSize())) {
		string error;
		error = "Requested start and end positions are outside the range of positions supplied in FitPositionsGaussian::fit_positions";
		throw std::range_error(error);
	}
	
	if (startPos > endPos) {
		string error;
		error = "Start is beyond end in FitPositionsGaussian::fit_positions";
		throw std::range_error(error);
	}
	
	size_t number_of_positions = endPos - startPos + 1;
	size_t size_of_subset = 2 * cutoff_radius + 1;
	size_t x_offset, y_offset, x_max, y_max;
	size_t number_of_intensities = size_of_subset * size_of_subset;
	size_t xSize = image->getXSize();
	size_t ySize = image->getYSize();
	
	double x0_initial, y0_initial, amplitude, background;
	double chi, degreesOfFreedom, c;
	long iterations = 0;
	int status;
	
	boost::shared_ptr<PALMMatrix<double> > image_subset;
	boost::shared_ptr<std::vector<LocalizedPosition> > fitted_positions (new std::vector<LocalizedPosition>);
	LocalizedPosition localizationResult;
	
	image_subset = boost::shared_ptr<PALMMatrix<double> > (new PALMMatrix<double>(size_of_subset, size_of_subset));
	
	fitted_positions->reserve(number_of_positions);
	
	// initialize the solver
	const gsl_multifit_fdfsolver_type *solver;
	gsl_multifit_fdfsolver *fit_iterator;
	
	solver = gsl_multifit_fdfsolver_lmsder;
	measured_data_Gauss_fits fitData;
	gsl_multifit_function_fdf f;
	
	gsl_vector *fit_parameters = gsl_vector_alloc(4);
	if (fit_parameters == NULL) {
		string error;
		error = "unable to allocate fit_parameters in FitPositionsGaussian::fit_positions()\r";
		throw OUT_OF_MEMORY(error);
	}
	
	fit_iterator = gsl_multifit_fdfsolver_alloc(solver, number_of_intensities, 4);
	if (fit_iterator == NULL) {
		gsl_vector_free(fit_parameters);
		string error;
		error = "unable to allocate fit_iterator in FitPositionsGaussian::fit_positions()\r";
		throw OUT_OF_MEMORY(error);
	}
	
	gsl_matrix *covarianceMatrix = gsl_matrix_alloc(4, 4);
	if (covarianceMatrix == NULL) {
		gsl_vector_free(fit_parameters);
		gsl_multifit_fdfsolver_free(fit_iterator);
		string error;
		error = "unable to allocate covarianceMatrix in FitPositionsGaussian::fit_positions()\r";
		throw OUT_OF_MEMORY(error);
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
		
		amplitude = positions->get(i, 0);
		x0_initial = positions->get(i, 1);
		y0_initial = positions->get(i, 2);
		background = positions->get(i, 3);
		
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
		fitData.width = r_initial * 1.414213562373095;
		
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
			
			status = gsl_multifit_test_delta(fit_iterator->dx, fit_iterator->x, 10, 10);
		} while ((status = GSL_CONTINUE) && (iterations < 200));
		
		if (gsl_vector_get(fit_iterator->x, 0) <= 0) {	// reject fits that have negative amplitudes
			continue;
		}
		
		// are the reported positions within the window?
		if ((gsl_vector_get(fit_iterator->x, 2) < x_offset) || (gsl_vector_get(fit_iterator->x, 2) > x_max) || (gsl_vector_get(fit_iterator->x, 3) < y_offset) || (gsl_vector_get(fit_iterator->x, 3) > y_max)) {
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
		localizationResult.amplitude = gsl_vector_get(fit_iterator->x, 0);
		localizationResult.width = r_initial;
		localizationResult.xPos = gsl_vector_get(fit_iterator->x, 1);
		localizationResult.yPos = gsl_vector_get(fit_iterator->x, 2);
		localizationResult.offset = gsl_vector_get(fit_iterator->x, 3);
		
		localizationResult.amplitudeError = c * sqrt(gsl_matrix_get(covarianceMatrix, 0, 0));
		localizationResult.widthError = 0;
		localizationResult.xPosError = c * sqrt(gsl_matrix_get(covarianceMatrix, 1, 1));
		localizationResult.yPosError = c * sqrt(gsl_matrix_get(covarianceMatrix, 2, 2));
		localizationResult.offsetError = c * sqrt(gsl_matrix_get(covarianceMatrix, 3, 3));
		
		// store the number of iterations
		localizationResult.nIterations = iterations;
		
		fitted_positions->push_back(localizationResult);
	}
	
	gsl_multifit_fdfsolver_free(fit_iterator);
	gsl_vector_free(fit_parameters);
	
	return fitted_positions;
	
}


boost::shared_ptr<std::vector<LocalizedPosition> > FitPositionsMultiplication::fit_positions(const boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<PALMMatrix<double> > positions, 
																				 size_t startPos, size_t endPos) {
	
	// some safety checks
	if ((endPos >= positions->getXSize()) || (startPos >= positions->getXSize())) {
		string error;
		error = "Requested start and end positions are outside the range of positions supplied in FitPositionsGaussian::fit_positions";
		throw std::range_error(error);
	}
	
	if (startPos > endPos) {
		string error;
		error = "Start is beyond end in FitPositionsGaussian::fit_positions";
		throw std::range_error(error);
	}
	
	size_t number_of_positions = endPos - startPos + 1;
	size_t size_of_subset = 2 * cutoff_radius + 1;
	size_t xSize = image->getXSize();
	size_t ySize = image->getYSize();
	size_t x_offset, y_offset, x_max, y_max;
	
	double x0_initial, y0_initial, amplitude, background;
	size_t iterations = 0;
	
	double convergence_treshold_squared = convergence_threshold * convergence_threshold;
	double delta_squared = 10 * convergence_treshold_squared;	// this test the convergence of the position determined by the iteration
	// it is the distance between (xn-1, yn-1) and (xn, yn)
	// we initialize it to a value well over the treshold so that we will run at least two iterations
	double previous_position_x, previous_position_y;
	double current_x, current_y;
	
	boost::shared_ptr<PALMMatrix<double> > image_subset;
	boost::shared_ptr<PALMMatrix<double> > image_subset_mask;
	boost::shared_ptr<std::vector<LocalizedPosition> > fitted_positions (new std::vector<LocalizedPosition>);
	LocalizedPosition localizationResult;
	
	image_subset = boost::shared_ptr<PALMMatrix<double> >(new PALMMatrix<double>(size_of_subset, size_of_subset));
	image_subset_mask = boost::shared_ptr<PALMMatrix<double> > (new PALMMatrix<double>(size_of_subset, size_of_subset));
	
	fitted_positions->reserve(number_of_positions);
	
	for (size_t i = startPos; i <= endPos; ++i) {
		amplitude = positions->get(i, 0);
		x0_initial = positions->get(i, 1);
		y0_initial = positions->get(i, 2);
		background = positions->get(i, 3);
		
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
		
		while (delta_squared > convergence_treshold_squared) {
			previous_position_x = current_x;
			previous_position_y = current_y;
			
			++iterations;
			
			if (iterations > 100) {	// the multiplication is not converging, we should stop
				SetNaN64(&current_x);
				SetNaN64(&current_y);
				break;
			}
			
			multiply_with_gaussian(image_subset, image_subset_mask, current_x, current_y, r_initial, background, amplitude);
			determine_x_y_position(image_subset_mask, current_x, current_y);
			
			if (iterations == 1)	// this is the first iteration, we should not check for termination
				continue;
			
			delta_squared = (current_x - previous_position_x) * (current_x - previous_position_x) + (current_y - previous_position_y) * (current_y - previous_position_y);
		}
		
		delta_squared = 10 * convergence_treshold_squared;
		
		localizationResult.width = r_initial;
		localizationResult.xPos = (double)current_x + (double)x_offset;
		localizationResult.yPos = (double)current_y + (double)y_offset;
		localizationResult.nIterations = iterations;
		
		fitted_positions->push_back(localizationResult);
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
		throw DIMENSIONS_SHOULD_BE_EQUAL();
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


boost::shared_ptr<std::vector<LocalizedPosition> > FitPositionsCentroid::fit_positions(const boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<PALMMatrix<double> > positions, 
																		   size_t startPos, size_t endPos) {
	
	// some safety checks
	if ((endPos >= positions->getXSize()) || (startPos >= positions->getXSize())) {
		string error;
		error = "Requested start and end positions are outside the range of positions supplied in FitPositionsCentroid::fit_positions";
		throw std::range_error(error);
	}
	
	if (startPos > endPos) {
		string error;
		error = "Start is beyond end in FitPositionsCentroid::fit_positions";
		throw std::range_error(error);
	}
	
	size_t number_of_positions = endPos - startPos + 1;
	size_t xSize = image->getXSize();
	size_t ySize = image->getYSize();
	size_t x_offset, y_offset, x_max, y_max;
	
	size_t x0_initial, y0_initial;
	double current_x, current_y;
	double denominator;
	
	boost::shared_ptr<std::vector<LocalizedPosition> > fitted_positions (new std::vector<LocalizedPosition>);
	LocalizedPosition localizationResult;
	
	fitted_positions->reserve(number_of_positions);
	
	for (size_t i = startPos; i <= endPos; ++i) {
		x0_initial = positions->get(i, 1);
		y0_initial = positions->get(i, 2);
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
		
		localizationResult.xPos = current_x;
		localizationResult.yPos = current_y;
		
		fitted_positions->push_back(localizationResult);
	}
	
	return fitted_positions;
}
