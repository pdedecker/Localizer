/*
 Copyright 2008-2011 Peter Dedecker.
 
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

#ifndef PALM_ANALYSIS_SEGMENTATION_H
#define PALM_ANALYSIS_SEGMENTATION_H

#include <algorithm>
#include "boost/smart_ptr.hpp"
#include "boost/thread/mutex.hpp"
#include "boost/thread/shared_mutex.hpp"
#include <eigen3/Eigen/Eigen>
#include "PALM_analysis_storage.h"
#include "PALM_analysis_MatrixRecycler.h"
#include <fftw3.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_histogram.h>

#ifdef WITH_IGOR
#include "XOPStandardHeaders.h"
#endif

class ThresholdImage {
public:
	ThresholdImage() {;}
	virtual ~ThresholdImage() {;}
	
	virtual boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > do_thresholding(ImagePtr image) = 0;
};

class ThresholdImage_Direct : public ThresholdImage {
	// direct thresholding, based on an absolute threshold
public:
	ThresholdImage_Direct(double thresh_parameter) {threshold = thresh_parameter;}
	~ThresholdImage_Direct() {;}
	
	boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > do_thresholding(ImagePtr image);
	
protected:
	double threshold;
};

class ThresholdImage_Isodata : public ThresholdImage {
public:
	ThresholdImage_Isodata() {;}
	~ThresholdImage_Isodata() {;}
	
	boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > do_thresholding(ImagePtr image);
};

class ThresholdImage_Triangle : public ThresholdImage {
public:
	ThresholdImage_Triangle() {;}
	
	boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > do_thresholding(ImagePtr image);
};

class ConvolveMatricesWithFFTClass {
public:
	ConvolveMatricesWithFFTClass() {;}
	~ConvolveMatricesWithFFTClass();
	
	ImagePtr ConvolveMatrixWithSmallKernel(ImagePtr image, ImagePtr kernel);
	ImagePtr ConvolveMatricesWithFFT(ImagePtr image1, ImagePtr image2);
	ImagePtr ConvolveMatrixWithGivenFFT(ImagePtr image, boost::shared_ptr<fftw_complex> array2_FFT, size_t FFT_xSize2, size_t FFT_ySize2);
	
	// http://www.leptonica.com/convolution.html
	ImagePtr ConvolveMatrixWithFlatKernel(ImagePtr image, size_t kernelXSize, size_t kernelYSize);
	
	/**
	 * Get the FFT of a single image, possibly for later use in a convolution.
	 * The x and y size of the transformed image will be returned by reference in FFT_xSize and FFT_ySize.
	 */
	boost::shared_ptr<fftw_complex> DoForwardFFT(ImagePtr image);
	
	/**
	 * Calculate the reverse FFT
	 */
	ImagePtr DoReverseFFT(boost::shared_ptr<fftw_complex> array_FFT, size_t xSize, size_t ySize);
protected:
	static boost::mutex FFTWPlannerMutex;
};

class ThresholdImage_GLRT_FFT : public ThresholdImage {
public:
	ThresholdImage_GLRT_FFT(double PFA_param, double width_param);
	~ThresholdImage_GLRT_FFT() {;}
	
	boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > do_thresholding(ImagePtr image);
	
protected:
	void MakeKernels(size_t xSize, size_t ySize);
	
	ConvolveMatricesWithFFTClass matrixConvolver;
	
	double PFA;
	double gaussianWidth;
	
	boost::mutex kernelCalculationMutex;
	boost::shared_mutex segmentationCalculationMutex;
	
	int useFFT;
	
	size_t kernelXSize, kernelYSize;
	ImagePtr smallGaussianKernel;
	boost::shared_ptr<fftw_complex> GaussianKernelFFT;
	
	double sum_squared_Gaussian;
	size_t windowSize;
};

class ThresholdImage_SmoothSigma : public ThresholdImage {
public:
	ThresholdImage_SmoothSigma(double pdfStdDev, double multiplicationFactor) : _pdfStdDev(pdfStdDev), _multiplicationFactor(multiplicationFactor) {}
	~ThresholdImage_SmoothSigma() {;}
	
	boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > do_thresholding(ImagePtr image);
	
protected:
	void _makeKernels(double psfWidth);
	
	double _pdfStdDev;
	double _multiplicationFactor;
	
	ImagePtr _averageKernel;
	ImagePtr _smoothingKernel;
	double _avgKernelSum;
	
	ConvolveMatricesWithFFTClass matrixConvolver;
	
	boost::mutex _kernelCalculationMutex;
};

class ThresholdImage_Preprocessor {
public:
	ThresholdImage_Preprocessor() {;}
	virtual ~ThresholdImage_Preprocessor() {;}
	
	virtual ImagePtr do_preprocessing(ImagePtr image) = 0;
};

class ThresholdImage_Preprocessor_DoNothing : public ThresholdImage_Preprocessor {	// WARNING: we assume that both x and y are odd
public:
	ThresholdImage_Preprocessor_DoNothing() {;}
	~ThresholdImage_Preprocessor_DoNothing() {;}
	
	ImagePtr do_preprocessing(ImagePtr image) {return image;}
};

class ThresholdImage_Preprocessor_MedianFilter : public ThresholdImage_Preprocessor {	// WARNING: we assume that both x and y are odd
public:
	ThresholdImage_Preprocessor_MedianFilter(unsigned x, size_t y) {kernel_x_size = x; kernel_y_size = y;}
	~ThresholdImage_Preprocessor_MedianFilter() {;}
	
	ImagePtr do_preprocessing(ImagePtr image);
	
protected:
	size_t kernel_x_size, kernel_y_size;
};


class ThresholdImage_Preprocessor_GaussianSmoothing : public ThresholdImage_Preprocessor {
public:
	ThresholdImage_Preprocessor_GaussianSmoothing(double Gaussian_width) {width = Gaussian_width;}	// we calculate the convolution kernel once,
	// when the first image is supplied
	~ThresholdImage_Preprocessor_GaussianSmoothing() {;}
	
	ImagePtr do_preprocessing(ImagePtr image);
	
protected:
	void generate_Gaussian_kernel(size_t x_size, size_t y_size);
	ConvolveMatricesWithFFTClass matrixConvolver;
	
	boost::mutex generateKernelMutex;
	
	ImagePtr Gaussian_kernel;	// automatically initialized to a NULL pointer
	double width;
	size_t kernel_x_size;
	size_t kernel_y_size;
};


class ThresholdImage_Preprocessor_MeanFilter : public ThresholdImage_Preprocessor {	// WARNING: we assume that both x and y are odd
public:
	ThresholdImage_Preprocessor_MeanFilter(unsigned x, size_t y) {kernel_x_size = x; kernel_y_size = y;}
	~ThresholdImage_Preprocessor_MeanFilter() {;}
	
	ImagePtr do_preprocessing(ImagePtr image);
	
protected:
	size_t kernel_x_size, kernel_y_size;
};


class ThresholdImage_Postprocessor {
public:
	ThresholdImage_Postprocessor() {;}
	virtual ~ThresholdImage_Postprocessor() {;}
	
	virtual boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > do_postprocessing(boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > thresholded_image, ImagePtr image) = 0;
};

class ThresholdImage_Postprocessor_DoNothing : public ThresholdImage_Postprocessor {	// remove pixels that are considered to be 'on'
	// but do not have any active neighbours
public:
	ThresholdImage_Postprocessor_DoNothing() {;}
	~ThresholdImage_Postprocessor_DoNothing() {;}
	
	boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > do_postprocessing(boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > thresholded_image, ImagePtr image) {return thresholded_image;}
};

class ThresholdImage_Postprocessor_RemoveIsolatedPixels : public ThresholdImage_Postprocessor {	// remove pixels that are considered to be 'on'
	// but do not have any active neighbours
public:
	ThresholdImage_Postprocessor_RemoveIsolatedPixels() {;}
	~ThresholdImage_Postprocessor_RemoveIsolatedPixels() {;}
	
	boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > do_postprocessing(boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > thresholded_image, ImagePtr image);
};

gsl_histogram * make_histogram_from_matrix(ImagePtr image, size_t number_of_bins);

#endif

