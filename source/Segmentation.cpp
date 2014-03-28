/*
 Copyright 2008-2014 Peter Dedecker.
 
 This file is part of Localizer.
 
 Localizer is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 Localizer is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with Localizer.  If not, see <http://www.gnu.org/licenses/>.
 
 
 Additional permission under GNU GPL version 3 section 7
 
 If you modify this Program, or any covered work, by
 linking or combining it with libraries required for interaction
 with analysis programs such as Igor Pro or Matlab, or to acquire
 data from or control hardware related to an experimental measurement,
 the licensors of this Program grant you additional permission
 to convey the resulting work.
 */

#include "Segmentation.h"
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_blas.h>

std::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > ThresholdImage_Direct::do_thresholding(ImagePtr image) {
	
	size_t x_size, y_size;
	double current_value;
	
	x_size = image->rows();
	y_size = image->cols();
	
	std::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > thresholded_image(new Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>((int)x_size, (int)y_size));
	
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

std::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > ThresholdImage_Isodata::do_thresholding(ImagePtr image) {
	gsl_histogram *hist;
	std::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > threshold_image;
	
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
		threshold_image = std::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> >(new Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>((int)x_size, (int)y_size));
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

std::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > ThresholdImage_Triangle::do_thresholding(ImagePtr image) {
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
	
	std::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > threshold_image;
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
		threshold_image = std::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> >(new Eigen::Matrix<int , Eigen::Dynamic, Eigen::Dynamic>((int)x_size, (int)y_size));
		threshold_image->setConstant(0);
		return threshold_image;
	}
	
	// calculate the line that connects the maximum and highest-intensity value
	slope = (end_val - max_val) / (end_bin_double - max_bin_double);
	intercept = max_val - slope * max_bin_double;
	
	// catch an unlikely case where the connecting line is flat (the histogram is apparently uniform)
	if (slope == 0) {
		gsl_histogram_free(hist);
		threshold_image = std::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> >(new Eigen::Matrix<int , Eigen::Dynamic, Eigen::Dynamic>((int)x_size, (int)y_size));
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
	this->smallGaussianKernel = ImagePtr (new Image((int)this->windowSize, (int)this->windowSize));
	ImagePtr Gaussian_kernel(new Image((int)xSize, (int)ySize));
	
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

std::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > ThresholdImage_GLRT_FFT::do_thresholding(ImagePtr image) {
	// the code is based on a series of matlab files sent by Didier Marguet, corresponding author of the original paper
	size_t xSize = image->rows();
	size_t ySize = image->cols();
	
	std::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > threshold_image;
	ImagePtr averages;
	ImagePtr image_squared;
	ImagePtr summed_squares;
	ImagePtr null_hypothesis;
	ImagePtr image_Gaussian_convolved;
	ImagePtr hypothesis_test;
	
	double double_window_pixels = (double)this->windowSize * (double)this->windowSize;
	size_t half_window_size = this->windowSize / 2;
	int imageNeedsResizing = 0;
	
	threshold_image = std::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> >(new Eigen::Matrix<int , Eigen::Dynamic, Eigen::Dynamic>((int)xSize, (int)ySize));
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
		ImagePtr reducedImage(GetRecycledMatrix((int)xSize, (int)ySize), FreeRecycledMatrix);
		for (size_t j = 0; j < ySize; ++j) {
			for (size_t i = 0; i < xSize; ++i) {
				(*reducedImage)(i,j) = (*image)(i,j);
			}
		}
		image = reducedImage;	// modify the smart_ptr
	}
	
	image_squared = ImagePtr(GetRecycledMatrix((int)xSize, (int)ySize), FreeRecycledMatrix);
	
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
	*image_squared = (*image).array().square();
	
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
	*null_hypothesis = (*summed_squares) - (*averages).array().square().matrix() * double_window_pixels;
	
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
	(*hypothesis_test) = ((*image_Gaussian_convolved).array().square()) / (*null_hypothesis).array() * this->sum_squared_Gaussian;
	
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

std::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > ThresholdImage_SmoothSigma::do_thresholding(ImagePtr image) {
	if (image->rows() * image->cols() <= 1)
		throw std::runtime_error("image too small");
	
	int xSize = image->rows();
	int ySize = image->cols();
	
	// if the image has odd dimensions then the convolution will be throw an error
	// so in that case the segmentation will have to run on a slightly smaller image
	// by allocating the threshold_image first we make sure that the it will have
	// the correct (original) size
	std::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > thresholdedImage(new Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>(xSize, ySize));
	int imageNeedsResizing = 0;
	if (xSize % 2 != 0) {
		xSize -= 1;
		imageNeedsResizing = 1;
	}
	if (ySize % 2 != 0) {
		ySize -= 1;
		imageNeedsResizing = 1;
	}
	
	if (imageNeedsResizing == 1) {
		ImagePtr reducedImage(GetRecycledMatrix((int)xSize, (int)ySize), FreeRecycledMatrix);
		for (int j = 0; j < ySize; ++j) {
			for (int i = 0; i < xSize; ++i) {
				(*reducedImage)(i,j) = (*image)(i,j);
			}
		}
		image = reducedImage;	// modify the smart_ptr
	}
	
	// do the required kernels exist?
	{
		boost::lock_guard<boost::mutex> lock(_kernelCalculationMutex);
		if (_averageKernel.get() == NULL) {
			assert(_smoothingKernel.get() == NULL);
			int gaussKernelSize = ceil(_pdfStdDev * 3);
			if (gaussKernelSize % 2 == 0)
				gaussKernelSize += 1;
			
			int avgKernelSize = 3 * gaussKernelSize;
			
			int halfGaussKernelSize = gaussKernelSize / 2;
			int halfAvgKernelSize = avgKernelSize / 2;
			
			_smoothingKernel = ImagePtr(new Image(image->rows(), image->cols()));
			_averageKernel = ImagePtr(new Image(image->rows(), image->cols()));
			_smoothingKernel->setConstant(0.0);
			_averageKernel->setConstant(0.0);
			
			int imageCenterX = image->rows() / 2;
			int imageCenterY = image->cols() / 2;
			
			// gaussian kernel
			double offsetX, offsetY;
			for (int i = 0; i < gaussKernelSize; ++i) {
				for (int j = 0; j < gaussKernelSize; ++j) {
					offsetX = static_cast<double>(i - halfGaussKernelSize);
					offsetY = static_cast<double>(i - halfGaussKernelSize);
					(*_smoothingKernel)(i + imageCenterX - halfGaussKernelSize, j + imageCenterY - halfGaussKernelSize) = exp(- (offsetX * offsetX + offsetY * offsetY) / (2.0 * _pdfStdDev));
				}
			}
			
			// average kernel
			for (int i = 0; i < avgKernelSize; ++i) {
				for (int j = 0; j < avgKernelSize; ++j) {
					offsetX = static_cast<double>(i - halfAvgKernelSize);
					offsetY = static_cast<double>(i - halfAvgKernelSize);
					(*_averageKernel)(i + imageCenterX - halfAvgKernelSize, j + imageCenterY - halfAvgKernelSize) = 1 - exp(- (offsetX * offsetX + offsetY * offsetY) / (2.0 * _pdfStdDev));
				}
			}
			_avgKernelSum = _averageKernel->sum();
		}
	}
	
	// subtract the average
	ImagePtr averageConvolved = matrixConvolver.ConvolveMatricesWithFFT(image, _averageKernel);
	(*averageConvolved) /= _avgKernelSum;
	
	ImagePtr averageSubtracted(GetRecycledMatrix(image->rows(), image->cols()), FreeRecycledMatrix);
	(*averageSubtracted) = (*image) - (*averageConvolved);
	
	// and smooth it
	int nPointsInImage = image->rows() * image->cols();
	ImagePtr smoothedSubtracted = matrixConvolver.ConvolveMatricesWithFFT(averageSubtracted, _smoothingKernel);
	double average = smoothedSubtracted->sum() / static_cast<double>(nPointsInImage);
	(*smoothedSubtracted).array() -= average;
	
	// finally get the standard deviation
	double stdDev = 0;
	double averageForStdDev = (*smoothedSubtracted).sum() / static_cast<double>(nPointsInImage);
	
	double* imageData = image->data();
	for (int i = 0; i < nPointsInImage; ++i) {
		stdDev += (*imageData - averageForStdDev) * (*imageData - averageForStdDev);
		++imageData;
	}
	stdDev /= static_cast<double>(nPointsInImage - 1);
	stdDev = sqrt(stdDev);
	
	double threshold = _multiplicationFactor * stdDev;
	
	for (int j = 0; j < ySize; ++j) {
		for (int i = 0; i < xSize; ++i) {
			(*thresholdedImage)(i, j) = ((*smoothedSubtracted)(i, j) >= threshold) ? 255 : 0;
		}
	}
	
	return thresholdedImage;
}

ImagePtr ThresholdImage_Preprocessor_MedianFilter::do_preprocessing(ImagePtr image) {
	
	size_t kernel_size = kernel_x_size * kernel_y_size;
	size_t half_kernel_size_x = kernel_x_size / 2;
	size_t half_kernel_size_y = kernel_y_size / 2;
	size_t x_size = image->rows();
	size_t y_size = image->cols();
	size_t offset;
	double value, median;
	
	gsl_vector *median_environment;
	ImagePtr filtered_image;
	
	// allocate a gsl_vector with the correct size
	median_environment = gsl_vector_alloc(kernel_size);
	size_t sorted_center = kernel_size / 2;
	
	// make a copy of the image
	// this copy will be median-filtered
	// close to the edges (where the kernel doesn't fit we will not modify the image)
	filtered_image = ImagePtr(new Image((int)x_size, (int)y_size));
	
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
	
	ImagePtr Gaussian_window(new Image((int)window_size, (int)window_size));
	
	Gaussian_kernel = ImagePtr(new Image((int)x_size, (int)y_size));
	
	
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




ImagePtr ThresholdImage_Preprocessor_GaussianSmoothing::do_preprocessing(ImagePtr image) {
	
	ImagePtr filtered_image;
	int x_size = image->rows();
	int y_size = image->cols();
	
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


ImagePtr ThresholdImage_Preprocessor_MeanFilter::do_preprocessing(ImagePtr image) {
	
	size_t kernel_size = kernel_x_size * kernel_y_size;
	double double_kernel_pixels = (double)kernel_size;
	size_t half_kernel_size_x = kernel_x_size / 2;
	size_t half_kernel_size_y = kernel_y_size / 2;
	size_t x_size = image->rows();
	size_t y_size = image->cols();
	double mean;
	
	ImagePtr filtered_image;
	
	// make a copy of the image
	// this copy will be mean-filtered
	// close to the edges, where the kernel doesn't fit we will not modify the image
	filtered_image = ImagePtr(new Image((int)x_size, (int)y_size));
	
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

std::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > ThresholdImage_Postprocessor_RemoveIsolatedPixels::do_postprocessing(std::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > thresholded_image, ImagePtr image) {
	// we don't care about the edges, they are ignored anyway in the fitting
	size_t x_size = thresholded_image->rows();
	size_t y_size = thresholded_image->cols();
	unsigned char value;
	bool neighbour_found;
	
	std::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > processed_thresholded_image;
	
	processed_thresholded_image = std::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> >(new Eigen::Matrix<int , Eigen::Dynamic, Eigen::Dynamic>((int)x_size, (int)y_size));
	
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

gsl_histogram * make_histogram_from_matrix(ImagePtr image, size_t number_of_bins) {
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



