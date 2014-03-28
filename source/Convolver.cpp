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

#include "Convolver.h"

boost::mutex ConvolveMatricesWithFFTClass::FFTWPlannerMutex;

ConvolveMatricesWithFFTClass::~ConvolveMatricesWithFFTClass() {
}

ImagePtr ConvolveMatricesWithFFTClass::ConvolveMatrixWithSmallKernel(ImagePtr image, ImagePtr kernel) {
	// implement a direct convolution with sufficiently small kernel
	size_t kernelXSize = kernel->rows();
	size_t kernelYSize = kernel->cols();
	size_t imageXSize = image->rows();
	size_t imageYSize = image->cols();
	
	if ((kernelXSize % 2 != 1) || (kernelYSize % 2 != 1))
		throw std::runtime_error("Tried to use a kernel with even dimensions in ConvolveMatrixWithSmallKernel");
	
	if ((kernelXSize > 20) || (kernelYSize > 20))
		throw std::runtime_error("Tried to use a large kernel in ConvolveMatrixWithSmallKernel");
	
	ImagePtr convolvedImage(GetRecycledMatrix((int)imageXSize, (int)imageYSize), FreeRecycledMatrix);
	
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

ImagePtr ConvolveMatricesWithFFTClass::ConvolveMatricesWithFFT(ImagePtr image1, ImagePtr image2) {
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
	std::shared_ptr<fftw_complex> array1_FFT = DoForwardFFT(image1);
	std::shared_ptr<fftw_complex> array2_FFT = DoForwardFFT(image2);
	
	n_FFT_values = x_size1 * y_size1;
	nColumns = x_size1 / 2 + 1;	// the order of the dimensions will be swapped because of column-major ordering
	
	// now do the convolution
	for (size_t i = 0; i < n_FFT_values; i++) {
        //complex_value[0] = array1Ptr[2*i] * array2Ptr[2*i] - array1Ptr[2*i+1] * array2Ptr[2*i+1];
        //complex_value[1] = array1Ptr[2*i] * array2Ptr[2*i+1] + array1Ptr[2*i+1] * array2Ptr[2*i];
		complex_value[0] = array1_FFT.get()[i][0] * (array2_FFT.get()[i])[0] - array1_FFT.get()[i][1] * array2_FFT.get()[i][1];
		complex_value[1] = array1_FFT.get()[i][0] * array2_FFT.get()[i][1] + array1_FFT.get()[i][1] * array2_FFT.get()[i][0];
		
		// store the result in the first array
		// we add in a comb function so the origin is at the center of the image
        //array1Ptr[2*i] = (2.0 * (double)((i / nColumns + i % nColumns) % 2) - 1.0) * -1.0 * complex_value[0];
        //array2Ptr[2*i+1] = (2.0 * (double)((i / nColumns + i % nColumns) % 2) - 1.0) * -1.0 * complex_value[1];
		array1_FFT.get()[i][0] = (2.0 * (double)((i / nColumns + i % nColumns) % 2) - 1.0) * -1.0 * complex_value[0];
		array1_FFT.get()[i][1] = (2.0 * (double)((i / nColumns + i % nColumns) % 2) - 1.0) * -1.0 * complex_value[1];
	}
	
	ImagePtr convolved_image = DoReverseFFT(array1_FFT, x_size1, y_size1);
	
	return convolved_image;
	
}

ImagePtr ConvolveMatricesWithFFTClass::ConvolveMatrixWithGivenFFT(ImagePtr image, std::shared_ptr<fftw_complex> array2_FFT, size_t FFT_xSize2, size_t FFT_ySize2) {
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
	std::shared_ptr<fftw_complex> array1_FFT = DoForwardFFT(image);
	
	n_FFT_values = x_size1 * y_size1;
	nColumns = x_size1 / 2 + 1;	// the order of the dimensions will be swapped because of column-major ordering
	
	// now do the convolution
	for (size_t i = 0; i < n_FFT_values; i++) {
        //complex_value[0] = array1Ptr[2*i] * array2Ptr[2*i] - array1Ptr[2*i+1] * array2Ptr[2*i+1];
        //complex_value[0] = array1Ptr[2*i] * array2Ptr[2*i+1] - array1Ptr[2*i+1] * array2Ptr[2*i];
		complex_value[0] = (array1_FFT.get()[i][0] * array2_FFT.get()[i][0]) - (array1_FFT.get()[i][1] * array2_FFT.get()[i][1]);
		complex_value[1] = (array1_FFT.get()[i][0] * array2_FFT.get()[i][1]) + (array1_FFT.get()[i][1] * array2_FFT.get()[i][0]);
		
		// store the result in the first array
		// we add in a comb function so the origin is at the center of the image
        //array1Ptr[2*i] = (2.0 * (double)((i / nColumns + i % nColumns) % 2) - 1.0) * -1.0 * complex_value[0];
        //array1Ptr[2*i+1] = (2.0 * (double)((i / nColumns + i % nColumns) % 2) - 1.0) * -1.0 * complex_value[1];
		array1_FFT.get()[i][0] = (2.0 * (double)((i / nColumns + i % nColumns) % 2) - 1.0) * -1.0 * complex_value[0];
		array1_FFT.get()[i][1] = (2.0 * (double)((i / nColumns + i % nColumns) % 2) - 1.0) * -1.0 * complex_value[1];
	}
	
	ImagePtr convolved_image = DoReverseFFT(array1_FFT, x_size1, y_size1);
	
	return convolved_image;
}

ImagePtr ConvolveMatricesWithFFTClass::ConvolveMatrixWithFlatKernel(ImagePtr image, size_t kernelXSize, size_t kernelYSize) {
	size_t xSize = image->rows();
	size_t ySize = image->cols();
	
	if ((kernelXSize % 2 != 1) || (kernelYSize % 2 != 1)) {
		throw std::runtime_error("A kernel with even dimensions was passed to ConvolveMatrixWithFlatKernel");
	}
	
	// calculate an accumulated image
	ImagePtr accumulatedImage(GetRecycledMatrix((int)xSize, (int)ySize), FreeRecycledMatrix);
	ImagePtr convolvedImage(GetRecycledMatrix((int)xSize, (int)ySize), FreeRecycledMatrix);
	
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

std::shared_ptr<fftw_complex> ConvolveMatricesWithFFTClass::DoForwardFFT(ImagePtr image) {
	
	size_t xSize = image->rows();
	size_t ySize = image->cols();
	size_t nPixels = xSize * ySize;
	fftw_plan forwardPlan;
	
	std::shared_ptr<fftw_complex> array_FFT((fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nPixels), fftw_free);
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

ImagePtr ConvolveMatricesWithFFTClass::DoReverseFFT(std::shared_ptr<fftw_complex> array_FFT, size_t xSize, size_t ySize) {
	
	ImagePtr image(GetRecycledMatrix((int)xSize, (int)ySize), FreeRecycledMatrix);
	
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
