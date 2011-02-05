/*
 *  PALM_analysis_segmentation.h
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 27/09/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef PALM_ANALYSIS_SEGMENTATION_H
#define PALM_ANALYSIS_SEGMENTATION_H

#include <algorithm>
#include "boost/smart_ptr.hpp"
#include "boost/thread/mutex.hpp"
#include "boost/thread/shared_mutex.hpp"
#include <Eigen/Eigen>
#include "PALM_analysis_storage.h"
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
	
	virtual boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > do_thresholding(boost::shared_ptr<Eigen::MatrixXd> image) = 0;
};

class ThresholdImage_Direct : public ThresholdImage {
	// direct thresholding, based on an absolute threshold
public:
	ThresholdImage_Direct(double thresh_parameter) {threshold = thresh_parameter;}
	~ThresholdImage_Direct() {;}
	
	boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > do_thresholding(boost::shared_ptr<Eigen::MatrixXd> image);
	
protected:
	double threshold;
};

#ifdef WITH_IGOR
class ThresholdImage_Igor_Iterative : public ThresholdImage {
	// making use of a built-in Igor operation
public:
	ThresholdImage_Igor_Iterative() {;}
	~ThresholdImage_Igor_Iterative() {;}
	
	boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > do_thresholding(boost::shared_ptr<Eigen::MatrixXd> image);
	
protected:
	boost::mutex threadMutex;
};

class ThresholdImage_Igor_Bimodal : public ThresholdImage {
	// making use of a built-in Igor operation
public:
	ThresholdImage_Igor_Bimodal() {;}
	~ThresholdImage_Igor_Bimodal() {;}
	
	boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > do_thresholding(boost::shared_ptr<Eigen::MatrixXd> image);
	
protected:
	boost::mutex threadMutex;
};

class ThresholdImage_Igor_Adaptive : public ThresholdImage {
	// making use of a built-in Igor operation
public:
	ThresholdImage_Igor_Adaptive() {;}
	~ThresholdImage_Igor_Adaptive() {;}
	
	boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > do_thresholding(boost::shared_ptr<Eigen::MatrixXd> image);
	
protected:
	boost::mutex threadMutex;
};

class ThresholdImage_Igor_Fuzzy1 : public ThresholdImage {
	// making use of a built-in Igor operation
public:
	ThresholdImage_Igor_Fuzzy1() {;}
	~ThresholdImage_Igor_Fuzzy1() {;}
	
	boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > do_thresholding(boost::shared_ptr<Eigen::MatrixXd> image);
	
protected:
	boost::mutex threadMutex;
};

class ThresholdImage_Igor_Fuzzy2 : public ThresholdImage {
	// making use of a built-in Igor operation
public:
	ThresholdImage_Igor_Fuzzy2() {;}
	~ThresholdImage_Igor_Fuzzy2() {;}
	
	boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > do_thresholding(boost::shared_ptr<Eigen::MatrixXd> image);
	
protected:
	boost::mutex threadMutex;
};
#endif // WITH_IGOR

class ThresholdImage_Isodata : public ThresholdImage {
public:
	ThresholdImage_Isodata() {;}
	~ThresholdImage_Isodata() {;}
	
	boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > do_thresholding(boost::shared_ptr<Eigen::MatrixXd> image);
};

class ThresholdImage_Triangle : public ThresholdImage {
public:
	ThresholdImage_Triangle() {;}
	
	boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > do_thresholding(boost::shared_ptr<Eigen::MatrixXd> image);
};

class ConvolveMatricesWithFFTClass {
public:
	ConvolveMatricesWithFFTClass() {;}
	~ConvolveMatricesWithFFTClass();
	
	boost::shared_ptr<Eigen::MatrixXd> ConvolveMatrixWithSmallKernel(boost::shared_ptr<Eigen::MatrixXd> image, boost::shared_ptr<Eigen::MatrixXd> kernel);
	boost::shared_ptr<Eigen::MatrixXd> ConvolveMatricesWithFFT(boost::shared_ptr<Eigen::MatrixXd> image1, boost::shared_ptr<Eigen::MatrixXd> image2);
	boost::shared_ptr<Eigen::MatrixXd> ConvolveMatrixWithGivenFFT(boost::shared_ptr<Eigen::MatrixXd> image, boost::shared_ptr<fftw_complex> array2_FFT, size_t FFT_xSize2, size_t FFT_ySize2);
	
	// http://www.leptonica.com/convolution.html
	boost::shared_ptr<Eigen::MatrixXd> ConvolveMatrixWithFlatKernel(boost::shared_ptr<Eigen::MatrixXd> image, size_t kernelXSize, size_t kernelYSize);
	
	/**
	 * Get the FFT of a single image, possibly for later use in a convolution.
	 * The x and y size of the transformed image will be returned by reference in FFT_xSize and FFT_ySize.
	 */
	boost::shared_ptr<fftw_complex> DoForwardFFT(boost::shared_ptr<Eigen::MatrixXd> image);
	
	/**
	 * Calculate the reverse FFT
	 */
	boost::shared_ptr<Eigen::MatrixXd> DoReverseFFT(boost::shared_ptr<fftw_complex> array_FFT, size_t xSize, size_t ySize);
protected:
	static boost::mutex FFTWPlannerMutex;
};

class ThresholdImage_GLRT_FFT : public ThresholdImage {
public:
	ThresholdImage_GLRT_FFT(double PFA_param, double width_param);
	~ThresholdImage_GLRT_FFT() {;}
	
	boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > do_thresholding();
	boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > do_thresholding(boost::shared_ptr<Eigen::MatrixXd> image);
	
protected:
	void MakeKernels(size_t xSize, size_t ySize);
	
	ConvolveMatricesWithFFTClass matrixConvolver;
	
	double PFA;
	double gaussianWidth;
	
	boost::mutex kernelCalculationMutex;
	boost::shared_mutex segmentationCalculationMutex;
	
	int useFFT;
	
	size_t kernelXSize, kernelYSize;
	boost::shared_ptr<Eigen::MatrixXd> smallGaussianKernel;
	boost::shared_ptr<fftw_complex> GaussianKernelFFT;
	
	double sum_squared_Gaussian;
	size_t windowSize;
};

class ThresholdImage_Preprocessor {
public:
	ThresholdImage_Preprocessor() {;}
	virtual ~ThresholdImage_Preprocessor() {;}
	
	virtual boost::shared_ptr<Eigen::MatrixXd> do_preprocessing(boost::shared_ptr<Eigen::MatrixXd> image) = 0;
};

class ThresholdImage_Preprocessor_DoNothing : public ThresholdImage_Preprocessor {	// WARNING: we assume that both x and y are odd
public:
	ThresholdImage_Preprocessor_DoNothing() {;}
	~ThresholdImage_Preprocessor_DoNothing() {;}
	
	boost::shared_ptr<Eigen::MatrixXd> do_preprocessing(boost::shared_ptr<Eigen::MatrixXd> image) {return image;}
};

class ThresholdImage_Preprocessor_MedianFilter : public ThresholdImage_Preprocessor {	// WARNING: we assume that both x and y are odd
public:
	ThresholdImage_Preprocessor_MedianFilter(unsigned x, size_t y) {kernel_x_size = x; kernel_y_size = y;}
	~ThresholdImage_Preprocessor_MedianFilter() {;}
	
	boost::shared_ptr<Eigen::MatrixXd> do_preprocessing(boost::shared_ptr<Eigen::MatrixXd> image);
	
protected:
	size_t kernel_x_size, kernel_y_size;
};


class ThresholdImage_Preprocessor_GaussianSmoothing : public ThresholdImage_Preprocessor {
public:
	ThresholdImage_Preprocessor_GaussianSmoothing(double Gaussian_width) {width = Gaussian_width;}	// we calculate the convolution kernel once,
	// when the first image is supplied
	~ThresholdImage_Preprocessor_GaussianSmoothing() {;}
	
	boost::shared_ptr<Eigen::MatrixXd> do_preprocessing(boost::shared_ptr<Eigen::MatrixXd> image);
	
protected:
	void generate_Gaussian_kernel(size_t x_size, size_t y_size);
	ConvolveMatricesWithFFTClass matrixConvolver;
	
	boost::mutex generateKernelMutex;
	
	boost::shared_ptr<Eigen::MatrixXd> Gaussian_kernel;	// automatically initialized to a NULL pointer
	double width;
	size_t kernel_x_size;
	size_t kernel_y_size;
};


class ThresholdImage_Preprocessor_MeanFilter : public ThresholdImage_Preprocessor {	// WARNING: we assume that both x and y are odd
public:
	ThresholdImage_Preprocessor_MeanFilter(unsigned x, size_t y) {kernel_x_size = x; kernel_y_size = y;}
	~ThresholdImage_Preprocessor_MeanFilter() {;}
	
	boost::shared_ptr<Eigen::MatrixXd> do_preprocessing(boost::shared_ptr<Eigen::MatrixXd> image);
	
protected:
	size_t kernel_x_size, kernel_y_size;
};


class ThresholdImage_Postprocessor {
public:
	ThresholdImage_Postprocessor() {;}
	virtual ~ThresholdImage_Postprocessor() {;}
	
	virtual boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > do_postprocessing(boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > thresholded_image, boost::shared_ptr<Eigen::MatrixXd> image) = 0;
};

class ThresholdImage_Postprocessor_DoNothing : public ThresholdImage_Postprocessor {	// remove pixels that are considered to be 'on'
	// but do not have any active neighbours
public:
	ThresholdImage_Postprocessor_DoNothing() {;}
	~ThresholdImage_Postprocessor_DoNothing() {;}
	
	boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > do_postprocessing(boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > thresholded_image, boost::shared_ptr<Eigen::MatrixXd> image) {return thresholded_image;}
};

class ThresholdImage_Postprocessor_RemoveIsolatedPixels : public ThresholdImage_Postprocessor {	// remove pixels that are considered to be 'on'
	// but do not have any active neighbours
public:
	ThresholdImage_Postprocessor_RemoveIsolatedPixels() {;}
	~ThresholdImage_Postprocessor_RemoveIsolatedPixels() {;}
	
	boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > do_postprocessing(boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > thresholded_image, boost::shared_ptr<Eigen::MatrixXd> image);
};

gsl_histogram * make_histogram_from_matrix(boost::shared_ptr<Eigen::MatrixXd> image, size_t number_of_bins);

#endif

