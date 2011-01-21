/*
 *  PALM_analysis_segmentation.cpp
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 27/09/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "PALM_analysis_segmentation.h"

boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > ThresholdImage_Direct::do_thresholding(boost::shared_ptr<Eigen::MatrixXd> image) {
	
	size_t x_size, y_size;
	double current_value;
	
	x_size = image->rows();
	y_size = image->cols();
	
	boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > thresholded_image(new Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>((int)x_size, (int)y_size));
	
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
boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > ThresholdImage_Igor_Iterative::do_thresholding(boost::shared_ptr<Eigen::MatrixXd> image) {
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
	
	x_size = image->rows();
	y_size = image->cols();
	
	dimensionSizes[0] = (long)x_size;
	dimensionSizes[1] = (long)y_size;
	dimensionSizes[2] = 0;
	
	/****** this function is not reentrant, only one thread can be past this stage. Make sure that this is so by locking the mutex ******/
	boost::lock_guard<boost::mutex> locker(threadMutex);
	
	
	// try to make the wave for the temporary storage
	result = MDMakeWave(&tmp_storage_wave, "tmp_thresh_storage", NULL, dimensionSizes, NT_FP64, 1);
	if (result != 0) {
		throw result;
	}
	
	boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > thresholded_image (new Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>((int)x_size, (int)y_size));
	
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
	
	return thresholded_image;
}

boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > ThresholdImage_Igor_Bimodal::do_thresholding(boost::shared_ptr<Eigen::MatrixXd> image) {
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
	
	x_size = image->rows();
	y_size = image->cols();
	
	dimensionSizes[0] = (long)x_size;
	dimensionSizes[1] = (long)y_size;
	dimensionSizes[2] = 0;
	
	/****** this function is not reentrant, only one thread can be past this stage. Make sure that this is so by locking the mutex ******/
	boost::lock_guard<boost::mutex> locker(threadMutex);
	
	// try to make the wave for the temporary storage
	result = MDMakeWave(&tmp_storage_wave, "tmp_thresh_storage", NULL, dimensionSizes, NT_FP64, 1);
	if (result != 0) {
		throw result;
	}
	
	boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > thresholded_image (new Eigen::Matrix<int , Eigen::Dynamic, Eigen::Dynamic>((int)x_size, (int)y_size));
	
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
	
	return thresholded_image;
}

boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > ThresholdImage_Igor_Adaptive::do_thresholding(boost::shared_ptr<Eigen::MatrixXd> image) {
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
	boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > original_thresholded;
	boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > transposed_tresholded;
	
	waveHndl tmp_storage_wave;
	waveHndl thresholded_wave;
	
	x_size = image->rows();
	y_size = image->cols();
	
	// we make two images for the original and the transposed threshold
	original_thresholded = boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> >(new Eigen::Matrix<int , Eigen::Dynamic, Eigen::Dynamic>((int)x_size, (int)y_size));
	transposed_tresholded = boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> >(new Eigen::Matrix<int , Eigen::Dynamic, Eigen::Dynamic>((int)x_size, (int)y_size));
	
	dimensionSizes[0] = (long)x_size;
	dimensionSizes[1] = (long)y_size;
	dimensionSizes[2] = 0;
	
	/****** this function is not reentrant, only one thread can be past this stage. Make sure that this is so by locking the mutex ******/
	boost::lock_guard<boost::mutex> locker(threadMutex);
	
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
	
	return original_thresholded;
}

boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > ThresholdImage_Igor_Fuzzy1::do_thresholding(boost::shared_ptr<Eigen::MatrixXd> image) {
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
	
	x_size = image->rows();
	y_size = image->cols();
	
	dimensionSizes[0] = (long)x_size;
	dimensionSizes[1] = (long)y_size;
	dimensionSizes[2] = 0;
	
	/****** this function is not reentrant, only one thread can be past this stage. Make sure that this is so by locking the mutex ******/
	boost::lock_guard<boost::mutex> locker(threadMutex);
	
	// try to make the wave for the temporary storage
	result = MDMakeWave(&tmp_storage_wave, "tmp_thresh_storage", NULL, dimensionSizes, NT_FP64, 1);
	if (result != 0) {
		throw result;
	}
	
	boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > thresholded_image = boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> >(new Eigen::Matrix<int , Eigen::Dynamic, Eigen::Dynamic>((int)x_size, (int)y_size));
	
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
	
	return thresholded_image;	
}

boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > ThresholdImage_Igor_Fuzzy2::do_thresholding(boost::shared_ptr<Eigen::MatrixXd> image) {
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
	
	x_size = image->rows();
	y_size = image->cols();
	
	dimensionSizes[0] = (long)x_size;
	dimensionSizes[1] = (long)y_size;
	dimensionSizes[2] = 0;
	
	/****** this function is not reentrant, only one thread can be past this stage. Make sure that this is so by locking the mutex ******/
	boost::lock_guard<boost::mutex> locker(threadMutex);
	
	// try to make the wave for the temporary storage
	result = MDMakeWave(&tmp_storage_wave, "tmp_thresh_storage", NULL, dimensionSizes, NT_FP64, 1);
	if (result != 0) {
		throw result;
	}
	
	boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > thresholded_image (new Eigen::Matrix<int , Eigen::Dynamic, Eigen::Dynamic>((int)x_size, (int)y_size));
	
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
	
	return thresholded_image;
}
#endif // WITH_IGOR

boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > ThresholdImage_Isodata::do_thresholding(boost::shared_ptr<Eigen::MatrixXd> image) {
	gsl_histogram *hist;
	boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > threshold_image;
	
	size_t x_size = image->rows();
	size_t y_size = image->cols();
	
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
		threshold_image = boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> >(new Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>((int)x_size, (int)y_size));
		threshold_image->setConstant(0.0);
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

boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > ThresholdImage_Triangle::do_thresholding(boost::shared_ptr<Eigen::MatrixXd> image) {
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
	
	boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > threshold_image;
	size_t x_size = image->rows();
	size_t y_size = image->cols();
	
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
		threshold_image = boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> >(new Eigen::Matrix<int , Eigen::Dynamic, Eigen::Dynamic>((int)x_size, (int)y_size));
		threshold_image->setConstant(0);
		return threshold_image;
	}
	
	// calculate the line that connects the maximum and highest-intensity value
	slope = (end_val - max_val) / (end_bin_double - max_bin_double);
	intercept = max_val - slope * max_bin_double;
	
	// catch an unlikely case where the connecting line is flat (the histogram is apparently uniform)
	if (slope == 0) {
		gsl_histogram_free(hist);
		threshold_image = boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> >(new Eigen::Matrix<int , Eigen::Dynamic, Eigen::Dynamic>((int)x_size, (int)y_size));
		threshold_image->setConstant(0);
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

ThresholdImage_GLRT_FFT::ThresholdImage_GLRT_FFT(double PFA_param, double width_param) {
	this->PFA = PFA_param;
	this->gaussianWidth = width_param;
	this->kernelXSize = 0;
	this->kernelYSize = 0;
	
	// calculate the window size that will be used
	this->windowSize = ceil(4 * this->gaussianWidth);
	if ((this->windowSize % 2) == 0)	// window_size must be odd
		this->windowSize += 1;
}

void ThresholdImage_GLRT_FFT::MakeKernels(size_t xSize, size_t ySize) {
	
	this->kernelXSize = xSize;
	this->kernelYSize = ySize;
	double double_window_pixels = (double)this->windowSize * (double)this->windowSize;
	size_t half_window_size = this->windowSize / 2;	// integer division takes care of the floor() aspect
	
	size_t center_x = xSize / 2;
	size_t center_y = ySize / 2;
	double sum;
	double distance_x, distance_y;
	
	// calculate the Gaussian kernel
	boost::shared_ptr<Eigen::MatrixXd> Gaussian_kernel(new Eigen::MatrixXd((int)xSize, (int)ySize));
	
	sum = 0;
	Gaussian_kernel->setConstant(0.0);
	
	for (size_t i = center_x - half_window_size; i <= center_x + half_window_size; i++) {
		for (size_t j = center_y - half_window_size; j <= center_y + half_window_size; j++) {
			distance_x = (double)center_x - (double)i;
			distance_y = (double)center_y - (double)j;
			(*Gaussian_kernel)(i, j) = 1.0 / (1.77245385 * this->gaussianWidth) * exp(- 1.0 / (2.0 * this->gaussianWidth * this->gaussianWidth) * (distance_x * distance_x + distance_y * distance_y));
			sum += (*Gaussian_kernel)(i, j);
		}
	}
	
	// now we re-normalize this Gaussian matrix
	// at this point Gaussian_window becomes equal to 'gc' in the original matlab code
	sum /= double_window_pixels;
	this->sum_squared_Gaussian = 0;
	for (size_t i = center_x - half_window_size; i <= center_x + half_window_size; i++) {
		for (size_t j = center_y - half_window_size; j <= center_y + half_window_size; j++) {
			(*Gaussian_kernel)(i, j) = (*Gaussian_kernel)(i, j) - sum;
			this->sum_squared_Gaussian += (*Gaussian_kernel)(i, j) * (*Gaussian_kernel)(i, j);	// this is 'Sgc2' in the original code
		}
	}
	
	// now calculate the FFTs of the kernel
	this->GaussianKernelFFT = this->matrixConvolver.DoForwardFFT(Gaussian_kernel);
}

boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > ThresholdImage_GLRT_FFT::do_thresholding(boost::shared_ptr<Eigen::MatrixXd> image) {
	// the code is based on a series of matlab files sent by Didier Marguet, corresponding author of the original paper
	size_t xSize = image->rows();
	size_t ySize = image->cols();
	
	boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > threshold_image;
	boost::shared_ptr<Eigen::MatrixXd> averages;
	boost::shared_ptr<Eigen::MatrixXd> image_squared;
	boost::shared_ptr<Eigen::MatrixXd> summed_squares;
	boost::shared_ptr<Eigen::MatrixXd> null_hypothesis;
	boost::shared_ptr<Eigen::MatrixXd> image_Gaussian_convolved;
	boost::shared_ptr<Eigen::MatrixXd> hypothesis_test;
	
	double current_value;
	double double_window_pixels = (double)this->windowSize * (double)this->windowSize;
	size_t half_window_size = this->windowSize / 2;
	int imageNeedsResizing = 0;
	
	threshold_image = boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> >(new Eigen::Matrix<int , Eigen::Dynamic, Eigen::Dynamic>((int)xSize, (int)ySize));
	threshold_image->setConstant(0);
	
	// if the image has odd dimensions then the convolution will be throw an error
	// so in that case the segmentation will have to run on a slightly smaller image
	// by allocating the threshold_image first we make sure that the it will have
	// the correct (original) size
	
	if (xSize % 2 != 0) {
		xSize -= 1;
		imageNeedsResizing = 1;
	}
	if (ySize % 2 != 0) {
		ySize -= 1;
		imageNeedsResizing = 1;
	}
	
	if (imageNeedsResizing == 1) {
		boost::shared_ptr<Eigen::MatrixXd> reducedImage(new Eigen::MatrixXd((int)xSize, (int)ySize));
		for (size_t i = 0; i < xSize; ++i) {
			for (size_t j = 0; j < ySize; ++j) {
				(*reducedImage)(i,j) = (*image)(i,j);
			}
		}
		image = reducedImage;	// modify the smart_ptr
	}
	
	image_squared = boost::shared_ptr<Eigen::MatrixXd>(new Eigen::MatrixXd((int)xSize, (int)ySize));
	
	// do we have kernels of the appropriate size?
	this->kernelCalculationMutex.lock();	// make sure that the kernel cannot be modified simultaneously by another thread
	if ((this->GaussianKernelFFT.get() == NULL) || (this->kernelXSize != xSize) || (this->kernelYSize != ySize)) {
		// the kernels need to be created or updated
		// no other thread can be performing a calculation
		// while this thread modifies the kernels
		{
			boost::lock_guard<boost::shared_mutex> locker(this->segmentationCalculationMutex);
			this->MakeKernels(xSize, ySize);
		}
	}
	
	// many threads can run a calculation, but only one can modify
	this->segmentationCalculationMutex.lock_shared();
	this->kernelCalculationMutex.unlock();
	
	// calculate the square of the pixel values
	// we'll use this later
	*image_squared = (*image).cwise().square();
	
	// NULL HYPOTHESIS: there is no emitter at a certain position
	
	// start by estimating the mean at every position
	// in the original code this done by a convolution of a unity matrix with the window size. We now do this using an FFT-based approach
	
	
	// convolve the image with a "box function", that will get us the average
	averages = matrixConvolver.ConvolveMatrixWithFlatKernel(image, 2 * half_window_size + 1, 2 * half_window_size + 1);
	
	// do the same for the squared image
	summed_squares = matrixConvolver.ConvolveMatrixWithFlatKernel(image_squared, 2 * half_window_size + 1, 2 * half_window_size + 1);
	
	// normalize the result, so that we get averages
	(*averages) /= double_window_pixels;
	
	// now calculate the null hypothesis image. This is T_sig0_2 in the original matlab source
	// recycle image_squared since it's already allocated and we won't use it again
	null_hypothesis = image_squared;
	*null_hypothesis = (*summed_squares) - (*averages).cwise().square() * double_window_pixels;
	
	// calculate the hypothesis H1 that there is an emitter
	
	// now we need to again convolve this Gaussian_window ('gc') with the original image. 
	// we now do this using the FFT
	image_Gaussian_convolved = matrixConvolver.ConvolveMatrixWithGivenFFT(image, this->GaussianKernelFFT, this->kernelXSize, this->kernelYSize);
	
	this->segmentationCalculationMutex.unlock_shared();
	
	// now normalize this convolved image so that it becomes equal to 'alpha' in the original matlab code
	(*image_Gaussian_convolved) /= this->sum_squared_Gaussian;
	
	// recycle the memory that has been allocated for average
	hypothesis_test = averages;
	(*hypothesis_test) = ((*image_Gaussian_convolved).cwise().square()).cwise() / (*null_hypothesis) * this->sum_squared_Gaussian;
	
	// calculate the image that will determine whether to accept or reject the null hypothesis
	for (size_t k = half_window_size + 1; k < xSize - half_window_size; k++) {
		for (size_t l = half_window_size + 1; l < ySize - half_window_size; l++) {
			current_value = 1 - (*hypothesis_test)(k, l);
			if (current_value > 0.0) {
				if (- double_window_pixels * log(current_value) > PFA)
					(*threshold_image)(k, l) = 255;
			}
		}
	}
	
	return threshold_image;
}


boost::shared_ptr<Eigen::MatrixXd> ThresholdImage_Preprocessor_MedianFilter::do_preprocessing(boost::shared_ptr<Eigen::MatrixXd> image) {
	
	size_t kernel_size = kernel_x_size * kernel_y_size;
	size_t half_kernel_size_x = kernel_x_size / 2;
	size_t half_kernel_size_y = kernel_y_size / 2;
	size_t x_size = image->rows();
	size_t y_size = image->cols();
	size_t offset;
	double value, median;
	
	gsl_vector *median_environment;
	boost::shared_ptr<Eigen::MatrixXd> filtered_image;
	
	// allocate a gsl_vector with the correct size
	median_environment = gsl_vector_alloc(kernel_size);
	size_t sorted_center = kernel_size / 2;
	
	// make a copy of the image
	// this copy will be median-filtered
	// close to the edges (where the kernel doesn't fit we will not modify the image)
	filtered_image = boost::shared_ptr<Eigen::MatrixXd>(new Eigen::MatrixXd((int)x_size, (int)y_size));
	
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
	
	boost::shared_ptr<Eigen::MatrixXd> Gaussian_window(new Eigen::MatrixXd((int)window_size, (int)window_size));
	
	Gaussian_kernel = boost::shared_ptr<Eigen::MatrixXd>(new Eigen::MatrixXd((int)x_size, (int)y_size));
	
	
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
	Gaussian_kernel->setConstant(0.0);
	
	for (size_t i = center_x - half_window_size; i <= center_x + half_window_size; i++) {
		for (size_t j = center_y - half_window_size; j <= center_y + half_window_size; j++) {
			(*Gaussian_kernel)(i, j) = (*Gaussian_window)(i - center_x + half_window_size, j - center_y + half_window_size);
		}
	}
}




boost::shared_ptr<Eigen::MatrixXd> ThresholdImage_Preprocessor_GaussianSmoothing::do_preprocessing(boost::shared_ptr<Eigen::MatrixXd> image) {
	
	boost::shared_ptr<Eigen::MatrixXd> filtered_image;
	size_t x_size = image->rows();
	size_t y_size = image->cols();
	
	// do we already have a Gaussian kernel stored, or is this the first run?
	{
		boost::lock_guard<boost::mutex> locker(generateKernelMutex);
		if (Gaussian_kernel.get() == NULL) {	// we don't have a kernel, we need to generate it
			
			generate_Gaussian_kernel((int)x_size, (int)y_size);
			
		} else {	// we already have a kernel stored, is it the correct size?
			// if not we will calculate a new one
			if ((x_size != Gaussian_kernel->rows()) || (y_size != Gaussian_kernel->cols())) {
				generate_Gaussian_kernel((int)x_size, (int)y_size);
			}
		}
	}
	
	filtered_image = matrixConvolver.ConvolveMatricesWithFFT(image, Gaussian_kernel);
	
	return filtered_image;
}


boost::shared_ptr<Eigen::MatrixXd> ThresholdImage_Preprocessor_MeanFilter::do_preprocessing(boost::shared_ptr<Eigen::MatrixXd> image) {
	
	size_t kernel_size = kernel_x_size * kernel_y_size;
	double double_kernel_pixels = (double)kernel_size;
	size_t half_kernel_size_x = kernel_x_size / 2;
	size_t half_kernel_size_y = kernel_y_size / 2;
	size_t x_size = image->rows();
	size_t y_size = image->cols();
	double mean;
	
	boost::shared_ptr<Eigen::MatrixXd> filtered_image;
	
	// make a copy of the image
	// this copy will be mean-filtered
	// close to the edges, where the kernel doesn't fit we will not modify the image
	filtered_image = boost::shared_ptr<Eigen::MatrixXd>(new Eigen::MatrixXd((int)x_size, (int)y_size));
	
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

boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > ThresholdImage_Postprocessor_RemoveIsolatedPixels::do_postprocessing(boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > thresholded_image, boost::shared_ptr<Eigen::MatrixXd> image) {
	// we don't care about the edges, they are ignored anyway in the fitting
	size_t x_size = thresholded_image->rows();
	size_t y_size = thresholded_image->cols();
	unsigned char value;
	bool neighbour_found;
	
	boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > processed_thresholded_image;
	
	processed_thresholded_image = boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> >(new Eigen::Matrix<int , Eigen::Dynamic, Eigen::Dynamic>((int)x_size, (int)y_size));
	
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

boost::mutex ConvolveMatricesWithFFTClass::FFTWPlannerMutex;

ConvolveMatricesWithFFTClass::~ConvolveMatricesWithFFTClass() {
}

boost::shared_ptr<Eigen::MatrixXd> ConvolveMatrixWithSmallKernel(boost::shared_ptr<Eigen::MatrixXd> image, boost::shared_ptr<Eigen::MatrixXd> kernel) {
	// implement a direct convolution with sufficiently small kernel
	size_t kernelXSize = kernel->rows();
	size_t kernelYSize = kernel->cols();
	size_t imageXSize = image->rows();
	size_t imageYSize = image->rows();
	
	if ((kernelXSize % 2 != 1) || (kernelYSize % 2 != 1))
		throw std::runtime_error("Tried to use a kernel with even dimensions in ConvolveMatrixWithSmallKernel");
	
	if ((kernelXSize > 20) || (kernelYSize > 20))
		throw std::runtime_error("Tried to use a large kernel in ConvolveMatrixWithSmallKernel");
	
	boost::shared_ptr<Eigen::MatrixXd> convolvedImage(new Eigen::MatrixXd((int)imageXSize, (int)imageYSize));
	
	size_t halfKernelXSize = kernelXSize / 2;
	size_t halfKernelYSize = kernelYSize / 2;
	double sum;
	
	for (size_t i = halfKernelXSize; i < imageXSize - halfKernelXSize; ++i) {
		for (size_t j = halfKernelYSize; j < imageYSize - halfKernelYSize; ++j) {
			sum = 0.0;
			for (size_t k = 0; k < halfKernelXSize; ++k) {
				for (size_t l = 0; l < halfKernelYSize; ++l) {
					sum += (*image)(i - halfKernelXSize + k, j - halfKernelYSize + l) * (*kernel)(k, l);
				}
			}
			(*convolvedImage)(i, j) = sum;
		}
	}
	
	return convolvedImage;
}

boost::shared_ptr<Eigen::MatrixXd> ConvolveMatricesWithFFTClass::ConvolveMatricesWithFFT(boost::shared_ptr<Eigen::MatrixXd> image1, boost::shared_ptr<Eigen::MatrixXd> image2) {
	size_t x_size1, y_size1, x_size2, y_size2;
	
	x_size1 = image1->rows();
	y_size1 = image1->cols();
	x_size2 = image2->rows();
	y_size2 = image2->cols();
	
	size_t n_FFT_values, nColumns;
	
	fftw_complex complex_value;
	
	// are the dimensions equal?
	if ((x_size1 != x_size2) || (y_size1 != y_size2)) {
		std::string error("Tried to convolve images with unequal dimensions");
		throw DIMENSIONS_SHOULD_BE_EQUAL(error);
	}
	
	if ((x_size1 % 2 != 0) || (y_size1 % 2 != 0))
		throw std::runtime_error("tried to convolve images with uneven dimensions");
	
	// do the forward transforms
	boost::shared_ptr<fftw_complex> array1_FFT = DoForwardFFT(image1);
	boost::shared_ptr<fftw_complex> array2_FFT = DoForwardFFT(image2);
	
	n_FFT_values = x_size1 * y_size1;
	nColumns = y_size1 / 2 + 1;
	
	// now do the convolution
	for (size_t i = 0; i < n_FFT_values; i++) {
		complex_value[0] = array1_FFT.get()[i][0] * array2_FFT.get()[i][0] - array1_FFT.get()[i][1] * array2_FFT.get()[i][1];
		complex_value[1] = array1_FFT.get()[i][0] * array2_FFT.get()[i][1] + array1_FFT.get()[i][1] * array2_FFT.get()[i][0];
		
		// store the result in the first array
		// we add in a comb function so the origin is at the center of the image
		array1_FFT.get()[i][0] = (2.0 * (double)((i / nColumns + i % nColumns) % 2) - 1.0) * -1.0 * complex_value[0];
		array1_FFT.get()[i][1] = (2.0 * (double)((i / nColumns + i % nColumns) % 2) - 1.0) * -1.0 * complex_value[1];
	}
	
	boost::shared_ptr<Eigen::MatrixXd> convolved_image = DoReverseFFT(array1_FFT, x_size1, y_size1);
	
	return convolved_image;
	
}

boost::shared_ptr<Eigen::MatrixXd> ConvolveMatricesWithFFTClass::ConvolveMatrixWithGivenFFT(boost::shared_ptr<Eigen::MatrixXd> image, boost::shared_ptr<fftw_complex> array2_FFT, size_t FFT_xSize2, size_t FFT_ySize2) {
	size_t x_size1, y_size1;
	
	x_size1 = image->rows();
	y_size1 = image->cols();
	
	size_t n_FFT_values, nColumns;
	
	fftw_complex complex_value;
	
	if ((x_size1 != FFT_xSize2) || (y_size1 != FFT_ySize2)) {
		std::string error("Tried to convolve images with unequal dimensions with given FFT");
		throw DIMENSIONS_SHOULD_BE_EQUAL(error);
	}
	
	if ((x_size1 % 2 != 0) || (y_size1 % 2 != 0))
		throw std::runtime_error("tried to convolve images with uneven dimensions");
	
	// do the forward transforms
	boost::shared_ptr<fftw_complex> array1_FFT = DoForwardFFT(image);
	
	n_FFT_values = x_size1 * y_size1;
	nColumns = y_size1 / 2 + 1;
	
	// now do the convolution
	for (size_t i = 0; i < n_FFT_values; i++) {
		complex_value[0] = (array1_FFT.get()[i][0] * array2_FFT.get()[i][0]) - (array1_FFT.get()[i][1] * array2_FFT.get()[i][1]);
		complex_value[1] = (array1_FFT.get()[i][0] * array2_FFT.get()[i][1]) + (array1_FFT.get()[i][1] * array2_FFT.get()[i][0]);
		
		// store the result in the first array
		// we add in a comb function so the origin is at the center of the image
		array1_FFT.get()[i][0] = (2.0 * (double)((i / nColumns + i % nColumns) % 2) - 1.0) * -1.0 * complex_value[0];
		array1_FFT.get()[i][1] = (2.0 * (double)((i / nColumns + i % nColumns) % 2) - 1.0) * -1.0 * complex_value[1];
	}
	
	boost::shared_ptr<Eigen::MatrixXd> convolved_image = DoReverseFFT(array1_FFT, x_size1, y_size1);
	
	return convolved_image;
}

boost::shared_ptr<Eigen::MatrixXd> ConvolveMatricesWithFFTClass::ConvolveMatrixWithFlatKernel(boost::shared_ptr<Eigen::MatrixXd> image, size_t kernelXSize, size_t kernelYSize) {
	size_t xSize = image->rows();
	size_t ySize = image->cols();
	
	if ((kernelXSize % 2 != 1) || (kernelYSize % 2 != 1)) {
		throw std::runtime_error("A kernel with even dimensions was passed to ConvolveMatrixWithFlatKernel");
	}
	
	// calculate an accumulated image
	boost::shared_ptr<Eigen::MatrixXd> accumulatedImage(new Eigen::MatrixXd((int)xSize, (int)ySize));
	boost::shared_ptr<Eigen::MatrixXd> convolvedImage(new Eigen::MatrixXd((int)xSize, (int)ySize));
	
	// populate the upper row of the matrix
	memcpy(&(accumulatedImage->data()[0]), &(image->data()[0]), ySize * sizeof(double));
	
	// first loop: calculate the sum along the columns
	for (size_t i = 1; i < xSize; ++i) {
		for (size_t j = 0; j < ySize; ++j) {
			(*accumulatedImage)(i, j) = (*accumulatedImage)(i - 1, j) + (*image)(i, j);
		}
	}
	
	// second loop: calculate the sum of every pixel along the rows
	for (size_t i = 0; i < xSize; ++i) {
		for (size_t j = 1; j < ySize; ++j) {
			(*accumulatedImage)(i, j) = (*accumulatedImage)(i, j) + (*accumulatedImage)(i, j - 1);
		}
	}
	
	// do the actual convolution
	// the value of the convolution in the boundary region is undefined
	for (size_t i = kernelXSize / 2 + 1; i < xSize - (kernelXSize / 2); ++i) {
		for (size_t j = kernelYSize / 2 + 1; j < ySize - (kernelYSize / 2); ++j) {
			(*convolvedImage)(i, j) = (*accumulatedImage)(i + kernelXSize / 2, j + kernelXSize / 2) - (*accumulatedImage)(i - kernelXSize / 2 - 1, j + kernelXSize / 2)
			- (*accumulatedImage)(i + kernelXSize / 2, j - kernelXSize / 2 - 1) + (*accumulatedImage)(i - kernelXSize / 2 - 1, j - kernelXSize / 2 - 1);
		}
	}
	
	return convolvedImage;
}

boost::shared_ptr<fftw_complex> ConvolveMatricesWithFFTClass::DoForwardFFT(boost::shared_ptr<Eigen::MatrixXd> image) {
	
	size_t xSize = image->rows();
	size_t ySize = image->cols();
	size_t nPixels = xSize * ySize;
	fftw_plan forwardPlan;
	
	boost::shared_ptr<fftw_complex> array_FFT((fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nPixels), fftw_free);
	if (array_FFT.get() == NULL) {
		throw std::bad_alloc();
	}
	
	{
		boost::lock_guard<boost::mutex> locker(FFTWPlannerMutex);
		forwardPlan = fftw_plan_dft_r2c_2d((int)(xSize), (int)(ySize), &(image->data()[0]), array_FFT.get(), FFTW_ESTIMATE);
	}
	
	fftw_execute(forwardPlan);
	
	{
		boost::lock_guard<boost::mutex> locker(FFTWPlannerMutex);
		fftw_destroy_plan(forwardPlan);
	}
	
	return array_FFT;
}

boost::shared_ptr<Eigen::MatrixXd> ConvolveMatricesWithFFTClass::DoReverseFFT(boost::shared_ptr<fftw_complex> array_FFT, size_t xSize, size_t ySize) {
	
	boost::shared_ptr<Eigen::MatrixXd> image(new Eigen::MatrixXd((int)xSize, (int)ySize));
	
	double normalization_factor = (double)(xSize * ySize);
	fftw_plan reversePlan;
	
	{
		boost::lock_guard<boost::mutex> locker(FFTWPlannerMutex);
		reversePlan = fftw_plan_dft_c2r_2d((int)(xSize), (int)(ySize), array_FFT.get(), &(image->data()[0]), FFTW_ESTIMATE);
	}
	
	fftw_execute(reversePlan);
	
	{
		boost::lock_guard<boost::mutex> locker(FFTWPlannerMutex);
		fftw_destroy_plan(reversePlan);
	}
	
	*image /= normalization_factor;
	
	return image;
}


gsl_histogram * make_histogram_from_matrix(boost::shared_ptr<Eigen::MatrixXd> image, size_t number_of_bins) {
	size_t x_size, y_size;
	gsl_histogram *hist;
	double min = 1e100;
	double max = -1e100;
	double current_value;
	
	x_size = image->rows();
	y_size = image->cols();
	
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
