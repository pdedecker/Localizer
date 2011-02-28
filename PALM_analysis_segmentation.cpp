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
	
	for (size_t j = 0; j < y_size; j++) {
		for (size_t i = 0; i < x_size; i++) {
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
	this->useFFT = 1;
	
	// calculate the window size that will be used
	this->windowSize = ceil(4 * this->gaussianWidth);
	if ((this->windowSize % 2) == 0)	// window_size must be odd
		this->windowSize += 1;
}

void ThresholdImage_GLRT_FFT::MakeKernels(size_t xSize, size_t ySize) {
	
	this->kernelXSize = xSize;
	this->kernelYSize = ySize;
	double double_window_pixels = (double)this->windowSize * (double)this->windowSize;
	int half_window_size = this->windowSize / 2;	// integer division takes care of the floor() aspect
	
	size_t center_x = xSize / 2;
	size_t center_y = ySize / 2;
	double sum;
	double distance_x, distance_y;
	
	// decide whether to use the direct convolution or one based on the FFT
	int sizeIsPowerOfTwo = (((xSize & (xSize - 1)) == 0) && ((ySize & (ySize - 1)) == 0));
	if (sizeIsPowerOfTwo)
		this->useFFT = 1;
	else if (this->windowSize <= 7)
		this->useFFT = 0;
	else
		this->useFFT = 1;
	
	// calculate the Gaussian kernel
	// make both a small kernel suitable for direct convolution
	// and a large kernel (same size as the image) suitable for FFT-based convolution
	this->smallGaussianKernel = boost::shared_ptr<Eigen::MatrixXd> (new Eigen::MatrixXd((int)this->windowSize, (int)this->windowSize));
	boost::shared_ptr<Eigen::MatrixXd> Gaussian_kernel(new Eigen::MatrixXd((int)xSize, (int)ySize));
	
	sum = 0;
	Gaussian_kernel->setConstant(0.0);
	
	// start out by calculating the small kernel, then introduce it into the larger one
	for (int j = - 1.0 * half_window_size; j <= half_window_size; ++j) {
		for (int i = -1.0 * half_window_size; i <= half_window_size; ++i) {
			distance_x = (double)i;
			distance_y = (double)j;
			(*this->smallGaussianKernel)(i + half_window_size, j + half_window_size) = 1.0 / (1.77245385 * this->gaussianWidth) * exp(- 1.0 / (2.0 * this->gaussianWidth * this->gaussianWidth) * (distance_x * distance_x + distance_y * distance_y));
			sum += (*this->smallGaussianKernel)(i + half_window_size, j + half_window_size);
		}
	}
	
	// now we re-normalize this Gaussian matrix
	// at this point Gaussian_window becomes equal to 'gc' in the original matlab code
	sum /= double_window_pixels;
	this->sum_squared_Gaussian = 0;
	for (size_t j = 0; j < this->windowSize; ++j) {
		for (size_t i = 0; i < this->windowSize; ++i) {
			(*this->smallGaussianKernel)(i, j) = (*this->smallGaussianKernel)(i, j) - sum;
			this->sum_squared_Gaussian += (*this->smallGaussianKernel)(i, j) * (*this->smallGaussianKernel)(i, j);
		}
	}
	
	// introduce the small kernel into the big one
	for (size_t j = 0; j < this->windowSize; j++) {
		for (size_t i = 0; i < this->windowSize; i++) {
			(*Gaussian_kernel)(i + center_x - half_window_size, j + center_y - half_window_size) = (*this->smallGaussianKernel)(i, j);
		}
	}
	
	// now calculate the FFT of the big kernel
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
		boost::shared_ptr<Eigen::MatrixXd> reducedImage(GetSegmentationMatrix((int)xSize, (int)ySize), FreeSegmentationMatrix);
		for (size_t j = 0; j < ySize; ++j) {
			for (size_t i = 0; i < xSize; ++i) {
				(*reducedImage)(i,j) = (*image)(i,j);
			}
		}
		image = reducedImage;	// modify the smart_ptr
	}
	
	image_squared = boost::shared_ptr<Eigen::MatrixXd>(GetSegmentationMatrix((int)xSize, (int)ySize), FreeSegmentationMatrix);
	
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
	// dependingo on the image or kernel size choose either direct or FFT-based convolution
	if (this->useFFT == 1)
		image_Gaussian_convolved = matrixConvolver.ConvolveMatrixWithGivenFFT(image, this->GaussianKernelFFT, this->kernelXSize, this->kernelYSize);
	else
		image_Gaussian_convolved = matrixConvolver.ConvolveMatrixWithSmallKernel(image, this->smallGaussianKernel);
	
	this->segmentationCalculationMutex.unlock_shared();
	
	// now normalize this convolved image so that it becomes equal to 'alpha' in the original matlab code
	(*image_Gaussian_convolved) /= this->sum_squared_Gaussian;
	
	// recycle the memory that has been allocated for average
	hypothesis_test = averages;
	(*hypothesis_test) = ((*image_Gaussian_convolved).cwise().square()).cwise() / (*null_hypothesis) * this->sum_squared_Gaussian;
	
	// calculate the threshold value that will serve to accept or reject the presence of an emitter
	// details are in the GLRTSpeedup.tex file in the project folder
	double GLRTThreshold = 1.0 - exp(- PFA / double_window_pixels);
	
	// calculate the image that will determine whether to accept or reject the null hypothesis
	for (size_t l = half_window_size + 1; l < ySize - half_window_size; l++) {
		for (size_t k = half_window_size + 1; k < xSize - half_window_size; k++) {
			if (((*hypothesis_test)(k, l) < 1.0) && ((*hypothesis_test)(k, l) > GLRTThreshold)) {
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
	
	for (size_t j = 0; j < y_size; ++j) {
		for (size_t i = 0; i < x_size; ++i) {
			(*filtered_image)(i, j) = (*image)(i, j);
		}
	}
	
	
	// the main loop
	// for now we handle this using a direct, slow calculation
	
	for (size_t j = half_kernel_size_y; j < y_size - half_kernel_size_y; j++) {
		for (size_t i = half_kernel_size_x; i < x_size - half_kernel_size_x; i++) {
			
			offset = 0;
			for (size_t l = j - half_kernel_size_y; l <= j + half_kernel_size_y; l++) {
				for (size_t k = i - half_kernel_size_x; k <= i + half_kernel_size_x; k++) {
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
	for (size_t j = 0; j < window_size; j++) {
		for (size_t i = 0; i < window_size; i++) {
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
	
	for (size_t j = center_y - half_window_size; j <= center_y + half_window_size; j++) {
		for (size_t i = center_x - half_window_size; i <= center_x + half_window_size; i++) {
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
	
	for (size_t j = 0; j < y_size; ++j) {
		for (size_t i = 0; i < x_size; ++i) {
			(*filtered_image)(i, j) = (*image)(i, j);
		}
	}
	
	
	// the main loop
	// for now we handle this using a direct, slow calculation
	for (size_t j = half_kernel_size_y; j < y_size - half_kernel_size_y; j++) {
		for (size_t i = half_kernel_size_x; i < x_size - half_kernel_size_x; i++) {
			
			mean = 0;
			for (size_t l = j - half_kernel_size_y; l <= j + half_kernel_size_y; l++) {
				for (size_t k = i - half_kernel_size_x; k <= i + half_kernel_size_x; k++) {
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
	
	for (size_t j = 0; j < y_size; ++j) {
		for (size_t i = 0; i < x_size; ++i) {
			(*processed_thresholded_image)(i, j) = (*thresholded_image)(i, j);
		}
	}
	
	// we will return a copy
	for (size_t j = 1; j < y_size - 1; j++) {
		for (size_t i = 1; i < x_size - 1; i++) {
			
			value = (*thresholded_image)(i, j);
			if (value < 128) {	// this is an 'off' pixel
				continue;
			}
			
			neighbour_found = 0;
			for (size_t l = j - 1; l <= j + 1; l++) {
				for (size_t k = i - 1; k <= i + 1; k++) {
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

boost::shared_ptr<Eigen::MatrixXd> ConvolveMatricesWithFFTClass::ConvolveMatrixWithSmallKernel(boost::shared_ptr<Eigen::MatrixXd> image, boost::shared_ptr<Eigen::MatrixXd> kernel) {
	// implement a direct convolution with sufficiently small kernel
	size_t kernelXSize = kernel->rows();
	size_t kernelYSize = kernel->cols();
	size_t imageXSize = image->rows();
	size_t imageYSize = image->cols();
	
	if ((kernelXSize % 2 != 1) || (kernelYSize % 2 != 1))
		throw std::runtime_error("Tried to use a kernel with even dimensions in ConvolveMatrixWithSmallKernel");
	
	if ((kernelXSize > 20) || (kernelYSize > 20))
		throw std::runtime_error("Tried to use a large kernel in ConvolveMatrixWithSmallKernel");
	
	boost::shared_ptr<Eigen::MatrixXd> convolvedImage(GetSegmentationMatrix((int)imageXSize, (int)imageYSize), FreeSegmentationMatrix);
	
	size_t halfKernelXSize = kernelXSize / 2;
	size_t halfKernelYSize = kernelYSize / 2;
	double sum;
	
	for (size_t j = halfKernelYSize; j < imageYSize - halfKernelYSize; ++j) {
		for (size_t i = halfKernelXSize; i < imageXSize - halfKernelXSize; ++i) {
			sum = 0.0;
			for (size_t l = 0; l < kernelYSize; ++l) {
				for (size_t k = 0; k < kernelXSize; ++k) {
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
	nColumns = x_size1 / 2 + 1;	// the order of the dimensions will be swapped because of column-major ordering
	
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
	nColumns = x_size1 / 2 + 1;	// the order of the dimensions will be swapped because of column-major ordering
	
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
	boost::shared_ptr<Eigen::MatrixXd> accumulatedImage(GetSegmentationMatrix((int)xSize, (int)ySize), FreeSegmentationMatrix);
	boost::shared_ptr<Eigen::MatrixXd> convolvedImage(GetSegmentationMatrix((int)xSize, (int)ySize), FreeSegmentationMatrix);
	
	// populate the first column of the matrix
	memcpy(accumulatedImage->data(), image->data(), xSize * sizeof(double));
	
	// first loop: calculate the sum of every pixel along the rows
	for (size_t j = 1; j < ySize; ++j) {
		for (size_t i = 0; i < xSize; ++i) {
			(*accumulatedImage)(i, j) = (*image)(i, j) + (*accumulatedImage)(i, j - 1);
		}
	}
	
	// second loop: calculate the sum along the columns
	for (size_t j = 0; j < ySize; ++j) {
		for (size_t i = 1; i < xSize; ++i) {
			(*accumulatedImage)(i, j) = (*accumulatedImage)(i - 1, j) + (*accumulatedImage)(i, j);
		}
	}
	
	// do the actual convolution
	// the value of the convolution in the boundary region is undefined
	for (size_t j = kernelYSize / 2 + 1; j < ySize - (kernelYSize / 2); ++j) {
		for (size_t i = kernelXSize / 2 + 1; i < xSize - (kernelXSize / 2); ++i) {
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
		// pass the image dimensions in opposite order to account for column-major ordering
		forwardPlan = fftw_plan_dft_r2c_2d((int)(ySize), (int)(xSize), image->data(), array_FFT.get(), FFTW_ESTIMATE);
	}
	
	fftw_execute(forwardPlan);
	
	{
		boost::lock_guard<boost::mutex> locker(FFTWPlannerMutex);
		fftw_destroy_plan(forwardPlan);
	}
	
	return array_FFT;
}

boost::shared_ptr<Eigen::MatrixXd> ConvolveMatricesWithFFTClass::DoReverseFFT(boost::shared_ptr<fftw_complex> array_FFT, size_t xSize, size_t ySize) {
	
	boost::shared_ptr<Eigen::MatrixXd> image(GetSegmentationMatrix((int)xSize, (int)ySize), FreeSegmentationMatrix);
	
	double normalization_factor = (double)(xSize * ySize);
	fftw_plan reversePlan;
	
	{
		boost::lock_guard<boost::mutex> locker(FFTWPlannerMutex);
		// pass the image dimensions in opposite order to account for column-major ordering
		reversePlan = fftw_plan_dft_c2r_2d((int)(ySize), (int)(xSize), array_FFT.get(), image->data(), FFTW_ESTIMATE);
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
	
	for (size_t j = 0; j < y_size; j++) {
		for (size_t i = 0; i < x_size; i++) {
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
	for (size_t j = 0; j < y_size; j++) {
		for (size_t i = 0; i < x_size; i++) {
			current_value = (*image)(i, j);
			
			gsl_histogram_increment(hist, current_value);
			
		}
	}
	
	return hist;
}

/**
 * A global instance of MatrixRecycler to be used in the
 * segmentation
 */
boost::shared_ptr<MatrixRecycler> segmentationMatrixRecycler(new MatrixRecycler);

/**
 * A function that will handle allocation of memory from segmentationMatrixRecycler
 */
Eigen::MatrixXd* GetSegmentationMatrix(size_t nRows, size_t nCols) {
	return segmentationMatrixRecycler->getMatrix(nRows, nCols);
}

/**
 * A function that will handle freeing of memory from segmentationMatrixRecycler
 */
void FreeSegmentationMatrix(Eigen::MatrixXd* matrixToFree) {
	segmentationMatrixRecycler->freeMatrix(matrixToFree);
}


