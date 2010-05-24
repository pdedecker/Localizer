/*
 *  PALM_analysis_segmentation.cpp
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 27/09/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "PALM_analysis_segmentation.h"

boost::shared_ptr<ublas::matrix <unsigned char> > ThresholdImage_Direct::do_thresholding(boost::shared_ptr<ublas::matrix<double> > image) {
	
	size_t x_size, y_size;
	double current_value;
	
	x_size = image->size1();
	y_size = image->size2();
	
	boost::shared_ptr<ublas::matrix <unsigned char> > thresholded_image(new ublas::matrix<unsigned char>(x_size, y_size));
	
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {
			current_value = (*image)(i, j);
			if (current_value >= threshold) {
				(*thresholded_image)(i, j) = 255;
			} else {
				(*thresholded_image)(i, j) = 0;
			}
		}
	}
	
	return thresholded_image;
	
}

#ifdef WITH_IGOR
boost::shared_ptr<ublas::matrix <unsigned char> > ThresholdImage_Igor_Iterative::do_thresholding(boost::shared_ptr<ublas::matrix<double> > image) {
	// we use a built-in Igor operation to handle this case
	// the approach is quite clumsy: we create a temporary wave in Igor that has the same values as the frame we are interested in
	// the igor routine then runs on that frame
	// and next we copy the data back out into a new image
	
	// this is wasteful, but in this way the interface is the same as the other classes
	size_t x_size, y_size;
	long dimensionSizes[MAX_DIMENSIONS+1];
	long indices[MAX_DIMENSIONS];
	double value[2];
	int result;
	unsigned char threshold_result;
	
	waveHndl tmp_storage_wave;
	
	x_size = image->size1();
	y_size = image->size2();
	
	dimensionSizes[0] = (long)x_size;
	dimensionSizes[1] = (long)y_size;
	dimensionSizes[2] = 0;
	
	/****** this function is not reentrant, only one thread can be past this stage. Make sure that this is so by locking the mutex ******/
	threadMutex.lock();
	
	
	// try to make the wave for the temporary storage
	result = MDMakeWave(&tmp_storage_wave, "tmp_thresh_storage", NULL, dimensionSizes, NT_FP64, 1);
	if (result != 0) {
		throw result;
	}
	
	boost::shared_ptr<ublas::matrix <unsigned char> > thresholded_image = boost::shared_ptr<ublas::matrix <unsigned char> >(new ublas::matrix<unsigned char>(x_size, y_size));
	
	// now copy the data
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {
			indices[0] = i;
			indices[1] = j;
			
			value[0] = (*image)(i, j);
			
			result = MDSetNumericWavePointValue(tmp_storage_wave, indices, value) ;
			if (result != 0) {
				throw result;
			}
		}
	}
	
	// now let Igor handle the calculation of the threshold
	result = XOPSilentCommand("ImageThreshold /Q /O /M=1 tmp_thresh_storage");
	if (result != 0) {
		throw result;
	}
	
	// now copy the data back to the output matrix
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {
			indices[0] = i;
			indices[1] = j;
			
			result = MDGetNumericWavePointValue(tmp_storage_wave, indices, value) ;
			if (result != 0) {
				throw result;
			}
			
			threshold_result = (unsigned char)(value[0] + 0.5);
			
			(*thresholded_image)(i, j) = threshold_result;
		}
	}
	
	// now delete the temporary storage wave
	result = KillWave(tmp_storage_wave);
	if (result != 0) {
		throw result;
	}
	
	/****** release the mutex ******/
	threadMutex.unlock();
	
	return thresholded_image;
}

boost::shared_ptr<ublas::matrix <unsigned char> > ThresholdImage_Igor_Bimodal::do_thresholding(boost::shared_ptr<ublas::matrix<double> > image) {
	// we use a built-in Igor operation to handle this case
	// the approach is quite clumsy: we create a temporary wave in Igor that has the same values as the frame we are interested in
	// the igor routine then runs on that frame
	// and next we copy the data back out into a new image
	
	// this is wasteful, but the interface is the same as the other classes
	size_t x_size, y_size;
	long dimensionSizes[MAX_DIMENSIONS+1];
	long indices[MAX_DIMENSIONS];
	double value[2];
	int result;
	unsigned char threshold_result;
	
	waveHndl tmp_storage_wave;
	
	x_size = image->size1();
	y_size = image->size2();
	
	dimensionSizes[0] = (long)x_size;
	dimensionSizes[1] = (long)y_size;
	dimensionSizes[2] = 0;
	
	/****** this function is not reentrant, only one thread can be past this stage. Make sure that this is so by locking the mutex ******/
	threadMutex.lock();
	
	// try to make the wave for the temporary storage
	result = MDMakeWave(&tmp_storage_wave, "tmp_thresh_storage", NULL, dimensionSizes, NT_FP64, 1);
	if (result != 0) {
		throw result;
	}
	
	boost::shared_ptr<ublas::matrix <unsigned char> > thresholded_image = boost::shared_ptr<ublas::matrix <unsigned char> >(new ublas::matrix<unsigned char>(x_size, y_size));
	
	// now copy the data
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {
			indices[0] = i;
			indices[1] = j;
			
			value[0] = (*image)(i, j);
			
			result = MDSetNumericWavePointValue(tmp_storage_wave, indices, value) ;
			if (result != 0) {
				throw result;
			}
		}
	}
	
	// now let Igor handle the calculation of the threshold
	result = XOPSilentCommand("ImageThreshold /Q /O /M=2 tmp_thresh_storage");
	if (result != 0) {
		throw result;
	}
	
	// now copy the data back to the output matrix
	for (size_t i = 0; i < x_size; ++i) {
		for (size_t j = 0; j < y_size; ++j) {
			indices[0] = i;
			indices[1] = j;
			
			result = MDGetNumericWavePointValue(tmp_storage_wave, indices, value) ;
			if (result != 0) {
				throw result;
			}
			
			threshold_result = (unsigned char)(value[0] + 0.5);
			
			(*thresholded_image)(i, j) = threshold_result;
		}
	}
	
	// now delete the temporary storage wave
	result = KillWave(tmp_storage_wave);
	if (result != 0) {
		throw result;
	}
	
	/****** release the mutex ******/
	threadMutex.unlock();
	
	return thresholded_image;
}

boost::shared_ptr<ublas::matrix <unsigned char> > ThresholdImage_Igor_Adaptive::do_thresholding(boost::shared_ptr<ublas::matrix<double> > image) {
	// we use a built-in Igor operation to handle this case
	// the approach is quite clumsy: we create a temporary wave in Igor that has the same values as the frame we are interested in
	// the igor routine then runs on that frame
	// and next we copy the data back out into a new image
	
	// this is wasteful, but the interface is the same as the other classes
	size_t x_size, y_size;
	long dimensionSizes[MAX_DIMENSIONS+1];
	long indices[MAX_DIMENSIONS];
	double value[2];
	int result;
	unsigned char threshold_result;
	boost::shared_ptr<ublas::matrix <unsigned char> > original_thresholded;
	boost::shared_ptr<ublas::matrix <unsigned char> > transposed_tresholded;
	
	waveHndl tmp_storage_wave;
	waveHndl thresholded_wave;
	
	x_size = image->size1();
	y_size = image->size2();
	
	// we make two images for the original and the transposed threshold
	original_thresholded = boost::shared_ptr<ublas::matrix <unsigned char> >(new ublas::matrix<unsigned char>(x_size, y_size));
	transposed_tresholded = boost::shared_ptr<ublas::matrix <unsigned char> >(new ublas::matrix<unsigned char>(x_size, y_size));
	
	dimensionSizes[0] = (long)x_size;
	dimensionSizes[1] = (long)y_size;
	dimensionSizes[2] = 0;
	
	/****** this function is not reentrant, only one thread can be past this stage. Make sure that this is so by locking the mutex ******/
	threadMutex.lock();
	
	// try to make the wave for the temporary storage
	result = MDMakeWave(&tmp_storage_wave, "tmp_thresh_storage", NULL, dimensionSizes, NT_FP64, 1);
	if (result != 0) {
		throw result;
	}
	
	// now copy the data
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {
			indices[0] = i;
			indices[1] = j;
			
			value[0] = (*image)(i, j);
			
			result = MDSetNumericWavePointValue(tmp_storage_wave, indices, value) ;
			if (result != 0) {
				throw result;
			}
		}
	}
	
	// now let Igor handle the calculation of the threshold
	result = XOPSilentCommand("ImageThreshold /Q /M=3 tmp_thresh_storage");
	if (result != 0) {
		throw result;
	}
	
	
	thresholded_wave = FetchWave("M_ImageThresh");
	if (thresholded_wave == NULL) {
		throw NOWAV;
	}
	
	// now copy the thresholded image back to the first output wave
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {
			indices[0] = i;
			indices[1] = j;
			
			result = MDGetNumericWavePointValue(thresholded_wave, indices, value) ;
			if (result != 0) {
				throw result;
			}
			
			threshold_result = (unsigned char)(value[0] + 0.5);
			
			(*original_thresholded)(i, j) = threshold_result;
		}
	}
	
	// now transpose the image
	result = XOPSilentCommand("MatrixTranspose tmp_thresh_storage");
	if (result != 0) {
		throw result;
	}
	
	// now calculate the threshold again
	result = XOPSilentCommand("ImageThreshold /Q /M=3 tmp_thresh_storage");
	if (result != 0) {
		throw result;
	}
	
	// we transpose the thresholded image back to the original orientation
	result = XOPSilentCommand("MatrixTranspose M_ImageThresh");
	if (result != 0) {
		throw result;
	}
	
	thresholded_wave = FetchWave("M_ImageThresh");
	if (thresholded_wave == NULL) {
		throw NOWAV;
	}
	
	// now copy the thresholded image back to the second output wave
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {
			indices[0] = i;
			indices[1] = j;
			
			result = MDGetNumericWavePointValue(thresholded_wave, indices, value) ;
			if (result != 0) {
				throw result;
			}
			
			threshold_result = (unsigned char)(value[0] + 0.5);
			
			(*transposed_tresholded)(i, j) = threshold_result;
		}
	}
	
	// now construct the combined threshold image by AND'ing the original and transpose together
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {
			if ((*original_thresholded)(i, j) < 128) {	// below the threshold, we skip it
				continue;
			} else {
				// is the transposed matrix also above the threshold?
				if ((*transposed_tresholded)(i, j) < 128) {	// below the threshold. We should not include this point
					(*original_thresholded)(i, j) = 0;
				} else {
					continue;
				}
			}
		}
	}
	
	
	
	// now delete the temporary storage wave
	result = KillWave(tmp_storage_wave);
	if (result != 0) {
		throw result;
	}
	
	/****** release the mutex ******/
	threadMutex.unlock();
	
	return original_thresholded;
}

boost::shared_ptr<ublas::matrix <unsigned char> > ThresholdImage_Igor_Fuzzy1::do_thresholding(boost::shared_ptr<ublas::matrix<double> > image) {
	// we use a built-in Igor operation to handle this case
	// the approach is quite clumsy: we create a temporary wave in Igor that has the same values as the frame we are interested in
	// the igor routine then runs on that frame
	// and next we copy the data back out into a new image
	
	// this is wasteful, but the interface is the same as the other classes
	size_t x_size, y_size;
	long dimensionSizes[MAX_DIMENSIONS+1];
	long indices[MAX_DIMENSIONS];
	double value[2];
	int result;
	unsigned char threshold_result;
	
	waveHndl tmp_storage_wave;
	
	x_size = image->size1();
	y_size = image->size2();
	
	dimensionSizes[0] = (long)x_size;
	dimensionSizes[1] = (long)y_size;
	dimensionSizes[2] = 0;
	
	/****** this function is not reentrant, only one thread can be past this stage. Make sure that this is so by locking the mutex ******/
	threadMutex.lock();
	
	// try to make the wave for the temporary storage
	result = MDMakeWave(&tmp_storage_wave, "tmp_thresh_storage", NULL, dimensionSizes, NT_FP64, 1);
	if (result != 0) {
		throw result;
	}
	
	boost::shared_ptr<ublas::matrix <unsigned char> > thresholded_image = boost::shared_ptr<ublas::matrix <unsigned char> >(new ublas::matrix<unsigned char>(x_size, y_size));
	
	// now copy the data
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {
			indices[0] = i;
			indices[1] = j;
			
			value[0] = (*image)(i, j);
			
			result = MDSetNumericWavePointValue(tmp_storage_wave, indices, value) ;
			if (result != 0) {
				throw result;
			}
		}
	}
	
	// now let Igor handle the calculation of the threshold
	result = XOPSilentCommand("ImageThreshold /Q /O /M=4 tmp_thresh_storage");
	if (result != 0) {
		throw result;
	}
	
	// now copy the data back to the output matrix
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {
			indices[0] = i;
			indices[1] = j;
			
			result = MDGetNumericWavePointValue(tmp_storage_wave, indices, value) ;
			if (result != 0) {
				throw result;
			}
			
			threshold_result = (unsigned char)(value[0] + 0.5);
			
			(*thresholded_image)(i, j) = threshold_result;
		}
	}
	
	// now delete the temporary storage wave
	result = KillWave(tmp_storage_wave);
	if (result != 0) {
		throw result;
	}
	
	/****** release the mutex ******/
	threadMutex.unlock();
	
	return thresholded_image;	
}

boost::shared_ptr<ublas::matrix <unsigned char> > ThresholdImage_Igor_Fuzzy2::do_thresholding(boost::shared_ptr<ublas::matrix<double> > image) {
	// we use a built-in Igor operation to handle this case
	// the approach is quite clumsy: we create a temporary wave in Igor that has the same values as the frame we are interested in
	// the igor routine then runs on that frame
	// and next we copy the data back out into a new image
	
	// this is wasteful, but the interface is the same as the other classes
	size_t x_size, y_size;
	long dimensionSizes[MAX_DIMENSIONS+1];
	long indices[MAX_DIMENSIONS];
	double value[2];
	int result;
	unsigned char threshold_result;
	
	waveHndl tmp_storage_wave;
	
	x_size = image->size1();
	y_size = image->size2();
	
	dimensionSizes[0] = (long)x_size;
	dimensionSizes[1] = (long)y_size;
	dimensionSizes[2] = 0;
	
	/****** this function is not reentrant, only one thread can be past this stage. Make sure that this is so by locking the mutex ******/
	threadMutex.lock();
	
	// try to make the wave for the temporary storage
	result = MDMakeWave(&tmp_storage_wave, "tmp_thresh_storage", NULL, dimensionSizes, NT_FP64, 1);
	if (result != 0) {
		throw result;
	}
	
	boost::shared_ptr<ublas::matrix <unsigned char> > thresholded_image = boost::shared_ptr<ublas::matrix <unsigned char> >(new ublas::matrix<unsigned char>(x_size, y_size));
	
	// now copy the data
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {
			indices[0] = i;
			indices[1] = j;
			
			value[0] = (*image)(i, j);
			
			result = MDSetNumericWavePointValue(tmp_storage_wave, indices, value) ;
			if (result != 0) {
				throw result;
			}
		}
	}
	
	// now let Igor handle the calculation of the threshold
	result = XOPSilentCommand("ImageThreshold /Q /O /M=5 tmp_thresh_storage");
	if (result != 0) {
		throw result;
	}
	
	// now copy the data back to the output matrix
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {
			indices[0] = i;
			indices[1] = j;
			
			result = MDGetNumericWavePointValue(tmp_storage_wave, indices, value) ;
			if (result != 0) {
				throw result;
			}
			
			threshold_result = (unsigned char)(value[0] + 0.5);
			
			(*thresholded_image)(i, j) = threshold_result;
		}
	}
	
	// now delete the temporary storage wave
	result = KillWave(tmp_storage_wave);
	if (result != 0) {
		throw result;
	}
	
	/****** release the mutex ******/
	threadMutex.unlock();
	
	return thresholded_image;
}
#endif // WITH_IGOR

boost::shared_ptr<ublas::matrix <unsigned char> > ThresholdImage_Isodata::do_thresholding(boost::shared_ptr<ublas::matrix<double> > image) {
	gsl_histogram *hist;
	boost::shared_ptr<ublas::matrix <unsigned char> > threshold_image;
	
	size_t x_size = image->size1();
	size_t y_size = image->size2();
	
	size_t number_of_bins = 256;
	size_t current_threshold_bin = 127;
	int n_iterations = 0;
	int max_iters = 50;
	double lower_mean, upper_mean;
	double sum, denominator_sum;
	double current_threshold = -1, previous_threshold;
	double lower_bin_limit, upper_bin_limit;
	size_t bin_threshold;
	double intensity_threshold;
	int converged = 0;
	
	// since this is a histogram-based approach we start by constructing the histogram
	hist = make_histogram_from_matrix(image, number_of_bins);
	
	// because this approach is based on thresholding it makes sense to only express the threshold in terms of bins, stored in "current_threshold_bin".
	// a value of 0 for "current_threshold_bin" means that all bins at index 0 and higher are considered to be 'signal', not 'background'.
	
	while ((converged == 0) && (n_iterations < max_iters)) {
		
		previous_threshold = current_threshold;
		
		// calculate the lower mean
		sum = 0;
		denominator_sum = 0;
		for (size_t i = 0; i < current_threshold_bin; i++) {
			sum += (double)i * gsl_histogram_get(hist, i);
			denominator_sum += gsl_histogram_get(hist, i);
		}
		
		lower_mean = sum / denominator_sum;
		
		// calculate the upper mean
		sum = 0;
		denominator_sum = 0;
		for (size_t i = current_threshold_bin; i < number_of_bins; i++) {
			sum += (double)i * gsl_histogram_get(hist, i);
			denominator_sum += gsl_histogram_get(hist, i);
		}
		
		upper_mean = sum / denominator_sum;
		
		current_threshold = (lower_mean + upper_mean) / 2.0;
		current_threshold_bin = floor(current_threshold);
		
		if (floor(current_threshold + 0.5) == floor(previous_threshold + 0.5)) {
			bin_threshold = floor(current_threshold + 0.5);
			converged = 1;
		}
		
	}
	
	if (converged == 0) {	// the iterations did not converge, there is no clear threshold
		// to indicate this we set everything to 'off' (0)
		threshold_image = boost::shared_ptr<ublas::matrix <unsigned char> >(new ublas::matrix<unsigned char>(x_size, y_size));
		threshold_image->clear();
		return threshold_image;
	}
	
	// now translate the threshold value to an intensity instead of being in bins
	gsl_histogram_get_range(hist, bin_threshold, &lower_bin_limit, &upper_bin_limit);
	intensity_threshold = lower_bin_limit;
	
	// get another threshold class to do the work for us
	ThresholdImage_Direct thresholder(intensity_threshold);
	
	try {
		threshold_image = thresholder.do_thresholding(image);
	}
	catch (std::bad_alloc) {
		gsl_histogram_free(hist);
		throw;
	}
	
	gsl_histogram_free(hist);
	
	return threshold_image;
}

boost::shared_ptr<ublas::matrix <unsigned char> > ThresholdImage_Triangle::do_thresholding(boost::shared_ptr<ublas::matrix<double> > image) {
	gsl_histogram *hist;
	size_t number_of_bins = 256;
	size_t maximum_bin;
	double max_val, max_bin_double;
	double end_val, end_bin_double;
	double slope, intercept, perpendicular_slope, perpendicular_intercept;
	double current_bin_value, double_i;
	double intercept_x, intercept_y;
	double distance;
	double max_distance = -1;
	size_t max_index;
	double lower_bin_limit, upper_bin_limit, intensity_threshold;
	
	boost::shared_ptr<ublas::matrix <unsigned char> > threshold_image;
	size_t x_size = image->size1();
	size_t y_size = image->size2();
	
	// since this is a histogram-based approach we start by constructing the histogram
	
	hist = make_histogram_from_matrix(image, number_of_bins);
	
	maximum_bin = gsl_histogram_max_bin(hist);
	max_bin_double = (double)maximum_bin;
	max_val = gsl_histogram_max_val(hist);
	end_val = gsl_histogram_get(hist, number_of_bins - 1);	// the bin that contains the largest intensity is the last bin in the histogram
	end_bin_double = (double)(number_of_bins - 1);
	
	// catch an unlikely case where the maximum corresponds to the last bin
	if (maximum_bin == (number_of_bins - 1)) {
		gsl_histogram_free(hist);
		threshold_image = boost::shared_ptr<ublas::matrix <unsigned char> >(new ublas::matrix<unsigned char>(x_size, y_size));
		threshold_image->clear();
		return threshold_image;
	}
	
	// calculate the line that connects the maximum and highest-intensity value
	slope = (end_val - max_val) / (end_bin_double - max_bin_double);
	intercept = max_val - slope * max_bin_double;
	
	// catch an unlikely case where the connecting line is flat (the histogram is apparently uniform)
	if (slope == 0) {
		gsl_histogram_free(hist);
		threshold_image = boost::shared_ptr<ublas::matrix <unsigned char> >(new ublas::matrix<unsigned char>(x_size, y_size));
		threshold_image->clear();
		return threshold_image;
	}
	
	
	// calculate the slope of a line perpendicular to the connecting line
	perpendicular_slope = - 1.0 / slope;
	
	for (size_t i = maximum_bin + 1; i < number_of_bins; i++) {	// determine the minimum distance in the triangle
		
		// what is the intercept for the perpendicular line if it has to go through the bin that we're currently looking at?
		current_bin_value = gsl_histogram_get(hist, i);
		double_i = (double)i;
		
		perpendicular_intercept = current_bin_value - perpendicular_slope * double_i;
		
		// where does the perpendicular line intercept the connecting line?
		// x = (b1 - b2) / (a1 - a2) and y = a1 * x + b1
		intercept_x = (intercept - perpendicular_intercept) / (perpendicular_slope - slope);
		intercept_y = slope * intercept_x + intercept;
		
		// what is the distance to the connecting line?
		distance = sqrt((intercept_x - double_i) * (intercept_x - double_i) + (intercept_y - current_bin_value) * (intercept_y - current_bin_value));
		
		if (distance > max_distance) {
			max_index = i;
			max_distance = distance;
		}
	}
	
	// translate the maximal index to a threshold value
	gsl_histogram_get_range(hist, max_index, &lower_bin_limit, &upper_bin_limit);
	intensity_threshold = lower_bin_limit;
	
	// get another threshold class to do the work for us
	ThresholdImage_Direct thresholder(intensity_threshold);
	
	try {
		threshold_image = thresholder.do_thresholding(image);
	}
	catch (std::bad_alloc) {
		gsl_histogram_free(hist);
		throw;
	}
	
	gsl_histogram_free(hist);
	
	return threshold_image;
}


boost::shared_ptr<ublas::matrix <unsigned char> > ThresholdImage_GLRT::do_thresholding(boost::shared_ptr<ublas::matrix<double> > image) {
	// the code is based on a series of matlab files sent by Didier Marguet, corresponding author of the original paper
	boost::shared_ptr<ublas::matrix <unsigned char> > threshold_image;
	size_t x_size = image->size1();
	size_t y_size = image->size2();
	
	boost::shared_ptr<ublas::matrix<double> > averages(new ublas::matrix<double> (x_size, y_size));	// contains an estimation of the average value at every position, calculation over the number of pixels in the window
	boost::shared_ptr<ublas::matrix<double> > image_squared(new ublas::matrix<double> (x_size, y_size));
	boost::shared_ptr<ublas::matrix<double> > summed_squares(new ublas::matrix<double> (x_size, y_size));
	boost::shared_ptr<ublas::matrix<double> > null_hypothesis(new ublas::matrix<double> (x_size, y_size));
	boost::shared_ptr<ublas::matrix<double> > Gaussian_window(new ublas::matrix<double> (x_size, y_size));
	boost::shared_ptr<ublas::matrix<double> > image_Gaussian_convolved(new ublas::matrix<double> (x_size, y_size));	// this is 'alpha' in the original matlab code
	boost::shared_ptr<ublas::matrix<double> > hypothesis_test(new ublas::matrix<double> (x_size, y_size));	// this is 'test' in the original matlab code
	
	
	double average = 0;
	size_t window_size = 13;	// this is the size of the window over which we calculate the hypotheses
	size_t half_window_size = window_size / 2;	// integer division takes care of the floor() aspect
	size_t window_pixels = window_size * window_size;
	double double_window_pixels = (double)(window_pixels);
	double current_value;
	double distance_x, distance_y;
	double sum;
	double sum_squared_Gaussian;
	
	
	threshold_image = boost::shared_ptr<ublas::matrix <unsigned char> >(new ublas::matrix<unsigned char>(x_size, y_size));
	threshold_image->clear();
	
	averages->clear();
	summed_squares->clear();
	hypothesis_test->clear();
	
	// calculate the square of the pixel values
	// we'll use this later
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {
			current_value = (*image)(i, j);
			current_value = current_value * current_value;
			(*image_squared)(i, j) = current_value;
		}
	}
	
	// NULL HYPOTHESIS: there is no emitter at a certain position
	
	// start by estimating the mean at every position
	// in the original code this done by a convolution of a unity matrix with the window size-> This is done using an FFT in the matlab code, but we'll do it directly
	
	for (size_t k = half_window_size; k < x_size - half_window_size; k++) {
		for (size_t l = half_window_size; l < y_size - half_window_size; l++) {
			// these two loops run over the entire image-> We exclude the edges where the window doesn't fit
			
			// now we loop over the size of the window, to determine the average at every pixel
			average = 0;
			
			for (size_t i = k - half_window_size; i <= k + half_window_size; i++) {
				for (size_t j = l - half_window_size; j <= l + half_window_size; j++) {
					average += (*image)(i, j);
				}
			}
			average /= double_window_pixels;
			
			(*averages)(k, l) = average;
		}
	}
	
	// now we do the same, but for the square of the pixel values, and we don't divide (no averaging)
	for (size_t k = half_window_size; k < x_size - half_window_size; k++) {
		for (size_t l = half_window_size; l < y_size - half_window_size; l++) {
			// these two loops run over the entire image-> We exclude the edges where the window doesn't fit
			current_value = 0;
			for (size_t i = k - half_window_size; i <= k + half_window_size; i++) {
				for (size_t j = l - half_window_size; j <= l + half_window_size; j++) {
					current_value += (*image_squared)(i, j);
				}
			}
			(*summed_squares)(k, l) = current_value;	// summed_squares is now equal to "Sim2" in the orignal matlab code
		}
	}
	
	// now calculate the null hypothesis image-> This is T_sig0_2 in the original matlab source
	for (size_t k = half_window_size; k < x_size - half_window_size; k++) {
		for (size_t l = half_window_size; l < y_size - half_window_size; l++) {
			current_value = (*summed_squares)(k, l) - double_window_pixels * (*averages)(k, l) * (*averages)(k, l);
			(*null_hypothesis)(k, l) = current_value;
		}
	}
	
	// calculate the hypothesis H1 that there is an emitter
	
	// store the values of a Gaussian in a matrix with the size of a window
	// this is so we can cache the values, and do not to calculate it every time
	sum = 0;
	for (size_t i = 0; i < window_size; i++) {
		for (size_t j = 0; j < window_size; j++) {
			// the Gaussian is assumed to be in the center of the window
			distance_x = (double)half_window_size - (double)i;
			distance_y = (double)half_window_size - (double)j;
			current_value = 1.0 / (1.77245385 * gaussianWidth) * exp(- 1.0 / (2.0 * gaussianWidth * gaussianWidth) * (distance_x * distance_x + distance_y * distance_y));
			
			(*Gaussian_window)(i, j) = current_value;
			
			sum += current_value;	// we will use this below
		}
	}
	
	// now we re-normalize this Gaussian matrix
	// at this point Gaussian_window becomes equal to 'gc' in the original matlab code
	sum /= double_window_pixels;
	sum_squared_Gaussian = 0;
	for (size_t i = 0; i < window_size; i++) {
		for (size_t j = 0; j < window_size; j++) {
			current_value = (*Gaussian_window)(i, j);
			current_value = current_value - sum;
			(*Gaussian_window)(i, j) = current_value;
			sum_squared_Gaussian += current_value * current_value;	// this is 'Sgc2' in the original code
		}
	}
	
	// now we need to again convolve this Gaussian_window ('gc') with the original image-> As before, we'll do it directly
	// these two loops run over the entire image-> We exclude the edges where the window doesn't fit
	for (size_t k = half_window_size; k < x_size - half_window_size; k++) {
		for (size_t l = half_window_size; l < y_size - half_window_size; l++) {
			
			sum = 0;
			for (size_t i = k - half_window_size; i <= k + half_window_size; i++) {
				for (size_t j = l - half_window_size; j <= l + half_window_size; j++) {
					current_value = (*Gaussian_window)(i - k + half_window_size, j - l + half_window_size);
					current_value *= (*image)(i, j);
					sum += current_value;
				}
			}
			sum /= sum_squared_Gaussian;
			(*image_Gaussian_convolved)(k, l) = sum;
		}
	}
	
	// "image_Gaussian_convolved" is now equal to 'alpha' in the original matlab code
	
	// calculate the image that will determine whether to accept or reject the null hypothesis
	for (size_t k = half_window_size; k < x_size - half_window_size; k++) {
		for (size_t l = half_window_size; l < y_size - half_window_size; l++) {
			current_value = 1 - (sum_squared_Gaussian * (*image_Gaussian_convolved)(k, l) * (*image_Gaussian_convolved)(k, l)) / (*null_hypothesis)(k , l);
			current_value = (current_value > 0) * current_value + (current_value <= 0);	// the equivalent of test = (test > 0) ->* test + (test <= 0) in the original code
			current_value = - double_window_pixels * log(current_value);
			(*hypothesis_test)(k, l) = current_value;
		}
	}
	
	// at this point 'hypothesis_test' is equal to 'carte_MV' in the original image
	// check where we have to reject the hypothesis test
	for (size_t k = half_window_size; k < x_size - half_window_size; k++) {
		for (size_t l = half_window_size; l < y_size - half_window_size; l++) {
			if ((*hypothesis_test)(k, l) > PFA) {
				(*threshold_image)(k, l) = 255;
			}
		}
	}
	
	return threshold_image;
}


boost::shared_ptr<ublas::matrix <unsigned char> > ThresholdImage_GLRT_FFT::do_thresholding(boost::shared_ptr<ublas::matrix<double> > image) {
	// the code is based on a series of matlab files sent by Didier Marguet, corresponding author of the original paper
	boost::shared_ptr<ublas::matrix <unsigned char> > threshold_image;
	size_t x_size = image->size1();
	size_t y_size = image->size2();
	
	boost::shared_ptr<ublas::matrix<double> > averages;
	boost::shared_ptr<ublas::matrix<double> > image_squared;
	boost::shared_ptr<ublas::matrix<double> > summed_squares;
	boost::shared_ptr<ublas::matrix<double> > null_hypothesis;
	boost::shared_ptr<ublas::matrix<double> > image_Gaussian_convolved;
	boost::shared_ptr<ublas::matrix<double> > hypothesis_test;
		
	size_t window_size = ceil(4 * this->gaussianWidth);
	if ((window_size % 2) == 0)	// window_size must be odd
		window_size += 1;
	
	size_t half_window_size = window_size / 2;	// integer division takes care of the floor() aspect
	size_t center_x = x_size / 2;
	size_t center_y = y_size / 2;
	size_t window_pixels = window_size * window_size;
	double double_window_pixels = (double)(window_pixels);
	double current_value;
	double distance_x, distance_y;
	double sum;
	
	threshold_image = boost::shared_ptr<ublas::matrix <unsigned char> >(new ublas::matrix<unsigned char>(x_size, y_size));
	threshold_image->clear();
	
	averages = boost::shared_ptr<ublas::matrix<double> >(new ublas::matrix<double>(x_size, y_size));
	averages->clear();
	
	image_squared = boost::shared_ptr<ublas::matrix<double> >(new ublas::matrix<double>(x_size, y_size));
	
	summed_squares = boost::shared_ptr<ublas::matrix<double> >(new ublas::matrix<double>(x_size, y_size));
	summed_squares->clear();
	
	null_hypothesis = boost::shared_ptr<ublas::matrix<double> >(new ublas::matrix<double>(x_size, y_size));
	
	image_Gaussian_convolved = boost::shared_ptr<ublas::matrix<double> >(new ublas::matrix<double>(x_size, y_size));	// this is 'alpha' in the original matlab code
	
	hypothesis_test = boost::shared_ptr<ublas::matrix<double> >(new ublas::matrix<double>(x_size, y_size));	// this is 'test' in the original matlab code
	hypothesis_test->clear();
	
	// calculate the square of the pixel values
	// we'll use this later
	ublas::noalias(*image_squared) = element_prod((*image), (*image));
	
	// NULL HYPOTHESIS: there is no emitter at a certain position
	
	// start by estimating the mean at every position
	// in the original code this done by a convolution of a unity matrix with the window size. We now do this using an FFT-based approach
	
	
	// convolve the image with a "box function", that will get us the average
	// if we have a lot of images then we only need to make this kernel once
	// but we need to make sure that only a single thread at a time is making the kernel, the others should block
	AverageKernelMutex.lock();
	
	if ((average_kernel.get() == NULL) || (averageKernelXSize != x_size) || (averageKernelYSize != y_size)) {
		averageCalculationMutex.lock();	// get a unique lock
		
		average_kernel = boost::shared_ptr<ublas::matrix<double> >(new ublas::matrix<double>(x_size, y_size));
		
		average_kernel->clear();
		for (size_t i = center_x - half_window_size; i <= center_x + half_window_size; i++) {
			for (size_t j = center_y - half_window_size; j <= center_y + half_window_size; j++) {
				(*average_kernel)(i, j) = 1;
			}
		}
		
		averageKernelXSize = x_size;
		averageKernelYSize = y_size;
		
		averageCalculationMutex.unlock();
	}
	
	averageCalculationMutex.lock_shared();
	AverageKernelMutex.unlock();
	
	averages = matrixConvolver.ConvolveMatricesWithFFT(image, average_kernel);
	
	// do the same for the squared image
	summed_squares = matrixConvolver.ConvolveMatricesWithFFT(image_squared, average_kernel);
	
	averageCalculationMutex.unlock_shared();
	
	// normalize the result, so that we get averages
	(*averages) /= double_window_pixels;
	
	// now calculate the null hypothesis image. This is T_sig0_2 in the original matlab source
	ublas::noalias(*null_hypothesis) = (*summed_squares) - (element_prod((*averages), (*averages)) * double_window_pixels);
	//(*null_hypothesis) = (*summed_squares) - (((*averages) * (*averages)).MultiplyWithScalar(double_window_pixels));
	
	// calculate the hypothesis H1 that there is an emitter
	
	// create a Gaussian kernel for the convolution
	// we only need to do this once if we are looking at a series of images
	// but make sure that only a single thread is doing this
	GaussianKernelMutex.lock();
	
	if ((Gaussian_kernel.get() == NULL) || (GaussianKernelXSize != x_size) || (GaussianKernelYSize != y_size)) {
		gaussianCalculationMutex.lock();	// get a unique lock
		Gaussian_kernel = boost::shared_ptr<ublas::matrix<double> >(new ublas::matrix<double>(x_size, y_size));
		
		sum = 0;
		Gaussian_kernel->clear();
		
		for (size_t i = center_x - half_window_size; i <= center_x + half_window_size; i++) {
			for (size_t j = center_y - half_window_size; j <= center_y + half_window_size; j++) {
				distance_x = (double)center_x - (double)i;
				distance_y = (double)center_y - (double)j;
				(*Gaussian_kernel)(i, j) = 1.0 / (1.77245385 * gaussianWidth) * exp(- 1.0 / (2.0 * gaussianWidth * gaussianWidth) * (distance_x * distance_x + distance_y * distance_y));
				sum += (*Gaussian_kernel)(i, j);
			}
		}
		
		// now we re-normalize this Gaussian matrix
		// at this point Gaussian_window becomes equal to 'gc' in the original matlab code
		sum /= double_window_pixels;
		sum_squared_Gaussian = 0;
		for (size_t i = center_x - half_window_size; i <= center_x + half_window_size; i++) {
			for (size_t j = center_y - half_window_size; j <= center_y + half_window_size; j++) {
				(*Gaussian_kernel)(i, j) = (*Gaussian_kernel)(i, j) - sum;
				sum_squared_Gaussian += (*Gaussian_kernel)(i, j) * (*Gaussian_kernel)(i, j);	// this is 'Sgc2' in the original code
			}
		}
		
		GaussianKernelXSize = x_size;
		GaussianKernelYSize = y_size;
		
		gaussianCalculationMutex.unlock();
	}
	
	// now we need to again convolve this Gaussian_window ('gc') with the original image. 
	// we now do this using the FFT
	
	gaussianCalculationMutex.lock_shared();
	GaussianKernelMutex.unlock();
	
	image_Gaussian_convolved = matrixConvolver.ConvolveMatricesWithFFT(image, Gaussian_kernel);
	gaussianCalculationMutex.unlock_shared();
	
	// now normalize this convolved image so that it becomes equal to 'alpha' in the original matlab code
	(*image_Gaussian_convolved) /= sum_squared_Gaussian;
	
	// calculate the image that will determine whether to accept or reject the null hypothesis
	for (size_t k = half_window_size; k < x_size - half_window_size; k++) {
		for (size_t l = half_window_size; l < y_size - half_window_size; l++) {
			current_value = 1 - (sum_squared_Gaussian * (*image_Gaussian_convolved)(k, l) * (*image_Gaussian_convolved)(k, l)) / (*null_hypothesis)(k , l);
			current_value = (current_value > 0) * current_value + (current_value <= 0);	// the equivalent of test = (test > 0) .* test + (test <= 0) in the original code
			current_value = - double_window_pixels * log(current_value);
			(*hypothesis_test)(k, l) = current_value;
		}
	}
	
	// at this point 'hypothesis_test' is equal to 'carte_MV' in the original image
	// check where we have to reject the hypothesis test
	for (size_t k = half_window_size; k < x_size - half_window_size; k++) {
		for (size_t l = half_window_size; l < y_size - half_window_size; l++) {
			if ((*hypothesis_test)(k, l) > PFA) {
				(*threshold_image)(k, l) = (unsigned char)255;
			}
		}
	}
	
	return threshold_image;
}


boost::shared_ptr<ublas::matrix<double> > ThresholdImage_Preprocessor_MedianFilter::do_preprocessing(boost::shared_ptr<ublas::matrix<double> > image) {
	
	size_t kernel_size = kernel_x_size * kernel_y_size;
	size_t half_kernel_size_x = kernel_x_size / 2;
	size_t half_kernel_size_y = kernel_y_size / 2;
	size_t x_size = image->size1();
	size_t y_size = image->size2();
	size_t offset;
	double value, median;
	
	gsl_vector *median_environment;
	boost::shared_ptr<ublas::matrix<double> > filtered_image;
	
	// allocate a gsl_vector with the correct size
	median_environment = gsl_vector_alloc(kernel_size);
	size_t sorted_center = kernel_size / 2;
	
	// make a copy of the image
	// this copy will be median-filtered
	// close to the edges (where the kernel doesn't fit we will not modify the image)
	filtered_image = boost::shared_ptr<ublas::matrix<double> >(new ublas::matrix<double>(x_size, y_size));
	
	for (size_t i = 0; i < x_size; ++i) {
		for (size_t j = 0; j < y_size; ++j) {
			(*filtered_image)(i, j) = (*image)(i, j);
		}
	}
	
	
	// the main loop
	// for now we handle this using a direct, slow calculation
	
	for (size_t i = half_kernel_size_x; i < x_size - half_kernel_size_x; i++) {
		for (size_t j = half_kernel_size_y; j < y_size - half_kernel_size_y; j++) {
			
			offset = 0;
			for (size_t k = i - half_kernel_size_x; k <= i + half_kernel_size_x; k++) {
				for (size_t l = j - half_kernel_size_y; l <= j + half_kernel_size_y; l++) {
					value = (*image)(k, l);
					gsl_vector_set(median_environment, offset, value);
					offset++;
				}
			}
			gsl_sort_vector(median_environment);
			median = gsl_vector_get(median_environment, sorted_center);
			(*filtered_image)(i, j) = median;
		}
	}
	
	gsl_vector_free(median_environment);
	
	return filtered_image;
}


void ThresholdImage_Preprocessor_GaussianSmoothing::generate_Gaussian_kernel(size_t x_size, size_t y_size) {
	
	size_t window_size = 31;
	size_t half_window_size = window_size / 2;
	size_t center_x = x_size / 2;
	size_t center_y = y_size / 2;
	double current_value, distance_x, distance_y;
	
	boost::shared_ptr<ublas::matrix<double> > Gaussian_window(new ublas::matrix<double>(window_size, window_size));
	
	Gaussian_kernel = boost::shared_ptr<ublas::matrix<double> >(new ublas::matrix<double>(x_size, y_size));
	
	
	// calculate the values of a Gaussian with the correct width in a smaller window
	for (size_t i = 0; i < window_size; i++) {
		for (size_t j = 0; j < window_size; j++) {
			// the Gaussian is assumed to be in the center of the window
			distance_x = (double)half_window_size - (double)i;
			distance_y = (double)half_window_size - (double)j;
			current_value = 1.0 / (6.28318531 * width * width) * exp(- 1.0 / (2.0 * width * width) * (distance_x * distance_x + distance_y * distance_y));
			// normalized Gaussian in two dimensions
			
			(*Gaussian_window)(i, j) = current_value;
		}
	}
	
	// now introduce this small kernel into a larger one that is the same size as the image
	Gaussian_kernel->clear();
	
	for (size_t i = center_x - half_window_size; i <= center_x + half_window_size; i++) {
		for (size_t j = center_y - half_window_size; j <= center_y + half_window_size; j++) {
			(*Gaussian_kernel)(i, j) = (*Gaussian_window)(i - center_x + half_window_size, j - center_y + half_window_size);
		}
	}
}




boost::shared_ptr<ublas::matrix<double> > ThresholdImage_Preprocessor_GaussianSmoothing::do_preprocessing(boost::shared_ptr<ublas::matrix<double> > image) {
	
	boost::shared_ptr<ublas::matrix<double> > filtered_image;
	size_t x_size = image->size1();
	size_t y_size = image->size2();
	
	// do we already have a Gaussian kernel stored, or is this the first run?
	generateKernelMutex.lock();
	if (Gaussian_kernel.get() == NULL) {	// we don't have a kernel, we need to generate it
		
		generate_Gaussian_kernel(x_size, y_size);
		
	} else {	// we already have a kernel stored, is it the correct size?
		// if not we will calculate a new one
		if ((x_size != Gaussian_kernel->size1()) || (y_size != Gaussian_kernel->size2())) {
			generate_Gaussian_kernel(x_size, y_size);
		}
	}
	generateKernelMutex.unlock();
	
	filtered_image = matrixConvolver.ConvolveMatricesWithFFT(image, Gaussian_kernel);
	
	return filtered_image;
}


boost::shared_ptr<ublas::matrix<double> > ThresholdImage_Preprocessor_MeanFilter::do_preprocessing(boost::shared_ptr<ublas::matrix<double> > image) {
	
	size_t kernel_size = kernel_x_size * kernel_y_size;
	double double_kernel_pixels = (double)kernel_size;
	size_t half_kernel_size_x = kernel_x_size / 2;
	size_t half_kernel_size_y = kernel_y_size / 2;
	size_t x_size = image->size1();
	size_t y_size = image->size2();
	double mean;
	
	boost::shared_ptr<ublas::matrix<double> > filtered_image;
	
	// make a copy of the image
	// this copy will be mean-filtered
	// close to the edges, where the kernel doesn't fit we will not modify the image
	filtered_image = boost::shared_ptr<ublas::matrix<double> >(new ublas::matrix<double>(x_size, y_size));
	
	for (size_t i = 0; i < x_size; ++i) {
		for (size_t j = 0; j < y_size; ++j) {
			(*filtered_image)(i, j) = (*image)(i, j);
		}
	}
	
	
	// the main loop
	// for now we handle this using a direct, slow calculation
	for (size_t i = half_kernel_size_x; i < x_size - half_kernel_size_x; i++) {
		for (size_t j = half_kernel_size_y; j < y_size - half_kernel_size_y; j++) {
			
			mean = 0;
			for (size_t k = i - half_kernel_size_x; k <= i + half_kernel_size_x; k++) {
				for (size_t l = j - half_kernel_size_y; l <= j + half_kernel_size_y; l++) {
					mean += (*image)(k, l);
				}
			}
			mean /= double_kernel_pixels;
			(*filtered_image)(i, j) = mean;
		}
	}
	
	return filtered_image;
}

boost::shared_ptr<ublas::matrix <unsigned char> > ThresholdImage_Postprocessor_RemoveIsolatedPixels::do_postprocessing(boost::shared_ptr<ublas::matrix <unsigned char> > thresholded_image, boost::shared_ptr<ublas::matrix<double> > image) {
	// we don't care about the edges, they are ignored anyway in the fitting
	size_t x_size = thresholded_image->size1();
	size_t y_size = thresholded_image->size2();
	unsigned char value;
	bool neighbour_found;
	
	boost::shared_ptr<ublas::matrix <unsigned char> > processed_thresholded_image;
	
	processed_thresholded_image = boost::shared_ptr<ublas::matrix <unsigned char> >(new ublas::matrix<unsigned char>(x_size, y_size));
	
	for (size_t i = 0; i < x_size; ++i) {
		for (size_t j = 0; j < y_size; ++j) {
			(*processed_thresholded_image)(i, j) = (*thresholded_image)(i, j);
		}
	}
	
	// we will return a copy
	for (size_t i = 1; i < x_size - 1; i++) {
		for (size_t j = 1; j < y_size - 1; j++) {
			
			value = (*thresholded_image)(i, j);
			if (value < 128) {	// this is an 'off' pixel
				continue;
			}
			
			neighbour_found = 0;
			for (size_t k = i - 1; k <= i + 1; k++) {
				for (size_t l = j - 1; l <= j + 1; l++) {
					if ((k == i) && (l == j)) {
						continue;	// this is the pixel that we are considering itself, not the environment
					}
					
					if ((*thresholded_image)(k, l) > 128) {
						neighbour_found = 1;
						break;
					}
				}
			}
			
			if (neighbour_found == 0) {
				// we didn't find an active point in the neighborhood, it was an isolated pixel
				(*processed_thresholded_image)(i, j) = 0;
			}
		}
	}
	
	return processed_thresholded_image;
}

ConvolveMatricesWithFFTClass::~ConvolveMatricesWithFFTClass() {
	if (forwardPlan != NULL) {
		fftw_destroy_plan(forwardPlan);
	}
	if (reversePlan != NULL) {
		fftw_destroy_plan(reversePlan);
	}
}

boost::shared_ptr<ublas::matrix<double> > ConvolveMatricesWithFFTClass::ConvolveMatricesWithFFT(boost::shared_ptr<ublas::matrix<double> > image1, boost::shared_ptr<ublas::matrix<double> > image2) {
	size_t x_size1, y_size1, x_size2, y_size2;
	
	x_size1 = image1->size1();
	y_size1 = image1->size2();
	x_size2 = image2->size1();
	y_size2 = image2->size2();
	
	size_t n_pixels, offset;
	size_t FFT_xSize, FFT_ySize;
	size_t n_FFT_values, nColumns;
	size_t lastRow, lastCol;
	
	double *array1;
	double *array2;
	fftw_complex *array1_FFT;
	fftw_complex *array2_FFT;
	fftw_complex complex_value;
	boost::shared_ptr<ublas::matrix<double> > convolved_image(new ublas::matrix<double>(x_size1, y_size1));
	
	// are the dimensions equal?
	if ((x_size1 != x_size2) || (y_size1 != y_size2)) {
		std::string error("Tried to convolve images with unequal dimensions");
		throw DIMENSIONS_SHOULD_BE_EQUAL(error);
	}
	
	// does the image have dimension sizes that are odd? if so remove one column and/or row so that it becomes even
	// we will correct for this when returning the image by copying the values back in
	if ((x_size1 % 2) == 1) {	// odd x size
		FFT_xSize = x_size1 - 1;
	} else {	// even size, nothing needs to happen
		FFT_xSize = x_size1;
	}
	
	if ((y_size1 % 2) == 1) {	// odd y size
		FFT_ySize = y_size1 - 1;
	} else {	// even size, nothing needs to happen
		FFT_ySize = y_size1;
	}
	
	n_pixels = FFT_xSize * FFT_ySize;
	double normalization_factor = (double)(n_pixels);
	
	// convert both of the matrices to an array of doubles
	// we allocate the memory using fftw_malloc()
	// as this is recommended by the library
	array1 = (double *)fftw_malloc(sizeof(double) * n_pixels);
	if (array1 == NULL) {
		throw std::bad_alloc();
	}
	
	array2 = (double *)fftw_malloc(sizeof(double) * n_pixels);
	if (array2 == NULL) {
		fftw_free(array1);
		throw std::bad_alloc();
	}
	
	offset = 0;
	// now copy the matrices to the arrays
	
	for (size_t i = 0; i < FFT_xSize; i++) {
		for (size_t j = 0; j < FFT_ySize; j++) {
			// IMPORTANT: the data in the array is assumed to be in ROW-MAJOR order, so we loop over y first
			array1[offset] = (*image1)(i, j);
			array2[offset] = (*image2)(i, j);
			
			offset++;
		}
	}
	
	// now allocate the arrays that will hold the transformed result
	// the dimensions of these arrays are a bit unusual, and are x_size * (y_size / 2 + 1)
	n_FFT_values = FFT_xSize * (FFT_ySize / 2 + 1);
	nColumns = FFT_ySize / 2 + 1;
	
	array1_FFT = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * n_FFT_values);
	if (array1_FFT == NULL) {
		fftw_free(array1);
		fftw_free(array2);
		throw std::bad_alloc();
	}
	array2_FFT = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * n_FFT_values);
	if (array2_FFT == NULL) {
		fftw_free(array1);
		fftw_free(array2);
		fftw_free(array1_FFT);
		throw std::bad_alloc();
	}
	
	// prepare the transform and execute it on the first array
	// if there is no forward plan yet then create it
	forwardPlanMutex.lock();
	if ((forwardPlan == NULL) || (forwardPlanXSize != FFT_xSize) || (forwardPlanYSize != FFT_ySize)) {
		forwardCalculationMutex.lock();	// require exclusive ownership
		if (forwardPlan != NULL) {
			fftw_destroy_plan(forwardPlan);
		}
		
		forwardPlanXSize = FFT_xSize;
		forwardPlanYSize = FFT_ySize;
		
		forwardPlan = fftw_plan_dft_r2c_2d((int)(FFT_xSize), (int)(FFT_ySize), array1, array1_FFT, FFTW_ESTIMATE);
		forwardCalculationMutex.unlock();
	}
	
	forwardCalculationMutex.lock_shared();
	forwardPlanMutex.unlock();
	
	fftw_execute_dft_r2c(forwardPlan, array1, array1_FFT);
	
	// do the same on the second array
	fftw_execute_dft_r2c(forwardPlan, array2, array2_FFT);
	
	forwardCalculationMutex.unlock_shared();
	
	// now do the convolution
	for (size_t i = 0; i < n_FFT_values; i++) {
		complex_value[0] = array1_FFT[i][0] * array2_FFT[i][0] - array1_FFT[i][1] * array2_FFT[i][1];
		complex_value[1] = array1_FFT[i][0] * array2_FFT[i][1] + array1_FFT[i][1] * array2_FFT[i][0];
		
		// store the result in the first array
		// we add in a comb function so the origin is at the center of the image
		array1_FFT[i][0] = (2.0 * (double)((i / nColumns + i % nColumns) % 2) - 1.0) * -1.0 * complex_value[0];
		array1_FFT[i][1] = (2.0 * (double)((i / nColumns + i % nColumns) % 2) - 1.0) * -1.0 * complex_value[1];
	}
	
	// now do the reverse transform
	// we overwrite the original array
	// if there is no reverse plan yet then create it
	reversePlanMutex.lock();
	if ((reversePlan == NULL) || (reversePlanXSize != FFT_xSize) || (reversePlanYSize != FFT_ySize)) {
		reverseCalculationMutex.lock();	// require exclusive ownership
		if (reversePlan != NULL) {
			fftw_destroy_plan(reversePlan);
		}
		
		reversePlanXSize = FFT_xSize;
		reversePlanYSize = FFT_ySize;
		
		reversePlan = fftw_plan_dft_c2r_2d((int)(FFT_xSize), (int)(FFT_ySize), array1_FFT, array1, FFTW_ESTIMATE);
		reverseCalculationMutex.unlock();
	}
	
	reverseCalculationMutex.lock_shared();
	reversePlanMutex.unlock();
	
	fftw_execute_dft_c2r(reversePlan, array1_FFT, array1);
	
	reverseCalculationMutex.unlock_shared();
	
	// and store the result (we don't overwrite the input arguments)
	offset = 0;
	for (size_t i = 0; i < FFT_xSize; i++) {
		for (size_t j = 0; j < FFT_ySize; j++) {
			// the data in the array is assumed to be in ROW-MAJOR order, so we loop over x first
			// we also normalize the result
			(*convolved_image)(i, j) = array1[offset] / normalization_factor;
			
			offset++;
		}
	}
	
	// if the number of rows was odd, make the last row (not included in the fft) a copy of that before it
	if ((x_size1 % 2) == 1) {
		lastRow = x_size1 - 1;
		for (size_t j = 0; j < y_size1; ++j) {
			(*convolved_image)(lastRow, j) = (*convolved_image)(lastRow - 1, j);
		}
	}
	
	// if the number of columns was odd, make the last column (not included in the fft) a copy of that before it
	if ((y_size1 % 2) == 1) {
		lastCol = y_size1 - 1;
		for (size_t i = 0; i < x_size1; ++i) {
			(*convolved_image)(i, lastCol) = (*convolved_image)(i, lastCol - 1);
		}
	}
	
	// if both the number of columns and the number of rows was odd, then the pixel at the top left (highest x, highest y) will be incorrect
	if (((x_size1 % 2) == 1) && ((y_size1 % 2) == 1)) {
		(*convolved_image)(x_size1 - 1, y_size1 - 1) = (*convolved_image)(x_size1 - 2, y_size1 - 2);
	}
	
	// cleanup
	fftw_free(array1);
	fftw_free(array2);
	fftw_free(array1_FFT);
	fftw_free(array2_FFT);
	
	return convolved_image;
	
}


gsl_histogram * make_histogram_from_matrix(boost::shared_ptr<ublas::matrix<double> > image, size_t number_of_bins) {
	size_t x_size, y_size;
	gsl_histogram *hist;
	double min = 1e100;
	double max = -1e100;
	double current_value;
	
	x_size = image->size1();
	y_size = image->size2();
	
	std::string error;
	error = "Unable to allocate a gsl_histogram in make_histogram_from_matrix()";
	
	hist = gsl_histogram_alloc(number_of_bins);
	if (hist == NULL) {
		throw std::bad_alloc();
	}
	
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {
			current_value = (*image)(i, j);
			if (current_value < min)
				min = current_value;
			if (current_value > max)
				max = current_value;
		}
	}
	
	// adjust the histogram bins so that they range uniformly from min to max and set the values to zero
	gsl_histogram_set_ranges_uniform(hist, min, max + 1e-10);
	
	// now populate the histogram
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {
			current_value = (*image)(i, j);
			
			gsl_histogram_increment(hist, current_value);
			
		}
	}
	
	return hist;
}
