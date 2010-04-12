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

#include "boost/smart_ptr.hpp"
#include "boost/thread/mutex.hpp"
#include "boost/thread/shared_mutex.hpp"
#include <boost/numeric/ublas/matrix.hpp>
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

namespace ublas = boost::numeric::ublas;

class ThresholdImage {
public:
	ThresholdImage() {;}
	virtual ~ThresholdImage() {;}
	
	virtual boost::shared_ptr<ublas::matrix <unsigned char> > do_thresholding(boost::shared_ptr<ublas::matrix<double> > image) = 0;
};

class ThresholdImage_Direct : public ThresholdImage {
	// direct thresholding, based on an absolute threshold
public:
	ThresholdImage_Direct(double thresh_parameter) {threshold = thresh_parameter;}
	~ThresholdImage_Direct() {;}
	
	boost::shared_ptr<ublas::matrix <unsigned char> > do_thresholding(boost::shared_ptr<ublas::matrix<double> > image);
	
protected:
	double threshold;
};

#ifdef WITH_IGOR
class ThresholdImage_Igor_Iterative : public ThresholdImage {
	// making use of a built-in Igor operation
public:
	ThresholdImage_Igor_Iterative() {;}
	~ThresholdImage_Igor_Iterative() {;}
	
	boost::shared_ptr<ublas::matrix <unsigned char> > do_thresholding(boost::shared_ptr<ublas::matrix<double> > image);
	
protected:
	boost::mutex threadMutex;
};

class ThresholdImage_Igor_Bimodal : public ThresholdImage {
	// making use of a built-in Igor operation
public:
	ThresholdImage_Igor_Bimodal() {;}
	~ThresholdImage_Igor_Bimodal() {;}
	
	boost::shared_ptr<ublas::matrix <unsigned char> > do_thresholding(boost::shared_ptr<ublas::matrix<double> > image);
	
protected:
	boost::mutex threadMutex;
};

class ThresholdImage_Igor_Adaptive : public ThresholdImage {
	// making use of a built-in Igor operation
public:
	ThresholdImage_Igor_Adaptive() {;}
	~ThresholdImage_Igor_Adaptive() {;}
	
	boost::shared_ptr<ublas::matrix <unsigned char> > do_thresholding(boost::shared_ptr<ublas::matrix<double> > image);
	
protected:
	boost::mutex threadMutex;
};

class ThresholdImage_Igor_Fuzzy1 : public ThresholdImage {
	// making use of a built-in Igor operation
public:
	ThresholdImage_Igor_Fuzzy1() {;}
	~ThresholdImage_Igor_Fuzzy1() {;}
	
	boost::shared_ptr<ublas::matrix <unsigned char> > do_thresholding(boost::shared_ptr<ublas::matrix<double> > image);
	
protected:
	boost::mutex threadMutex;
};

class ThresholdImage_Igor_Fuzzy2 : public ThresholdImage {
	// making use of a built-in Igor operation
public:
	ThresholdImage_Igor_Fuzzy2() {;}
	~ThresholdImage_Igor_Fuzzy2() {;}
	
	boost::shared_ptr<ublas::matrix <unsigned char> > do_thresholding(boost::shared_ptr<ublas::matrix<double> > image);
	
protected:
	boost::mutex threadMutex;
};
#endif // WITH_IGOR

class ThresholdImage_Isodata : public ThresholdImage {
public:
	ThresholdImage_Isodata() {;}
	~ThresholdImage_Isodata() {;}
	
	boost::shared_ptr<ublas::matrix <unsigned char> > do_thresholding(boost::shared_ptr<ublas::matrix<double> > image);
};

class ThresholdImage_Triangle : public ThresholdImage {
public:
	ThresholdImage_Triangle() {;}
	
	boost::shared_ptr<ublas::matrix <unsigned char> > do_thresholding(boost::shared_ptr<ublas::matrix<double> > image);
};

class ThresholdImage_GLRT : public ThresholdImage {
public:
	ThresholdImage_GLRT(double PFA_param, double width_param) {PFA = PFA_param; gaussianWidth = width_param;} 
	
	boost::shared_ptr<ublas::matrix <unsigned char> > do_thresholding(boost::shared_ptr<ublas::matrix<double> > image);
	
protected:
	double PFA;
	double gaussianWidth;
};

class ConvolveMatricesWithFFTClass {
public:
	ConvolveMatricesWithFFTClass() {forwardPlan = NULL; reversePlan = NULL; forwardPlanXSize = 0; forwardPlanYSize = 0; reversePlanXSize = 0; reversePlanYSize = 0;}
	~ConvolveMatricesWithFFTClass();
	
	boost::shared_ptr<ublas::matrix<double> > ConvolveMatricesWithFFT(boost::shared_ptr<ublas::matrix<double> > image1, boost::shared_ptr<ublas::matrix<double> > image2);
	
protected:
	fftw_plan forwardPlan;
	fftw_plan reversePlan;
	
	size_t forwardPlanXSize;
	size_t forwardPlanYSize;
	size_t reversePlanXSize;
	size_t reversePlanYSize;
	
	boost::mutex forwardPlanMutex;
	boost::mutex reversePlanMutex;
	
	boost::shared_mutex forwardCalculationMutex;
	boost::shared_mutex reverseCalculationMutex;
	
	/** the reason for the many mutexes is for thread safety: because the creation of fftw plans is not thread safe, only a single thread can create a plan
	 at any given time. This is assured by exclusively locking the planMutexes. 
	 But it's possible that one thread wants to create a plan, while another thread is already running a calculation with the previous one (for example images with different
	 sizes are being passed). In that case we need to make sure that the calculations with the previous plan are finished before we make a new one. We could do this with a 
	 standard mutex, but at the same time we need to allow different calculations with the same plan to proceed in parallel. So we use a shared mutex: every thread that enters
	 a calculation does a shared_lock() on the mutex, which is unlocked when the calculation is finished. But to create a plan we require exclusive ownership. This means that
	 when we get the lock all the calculation threads have finished, and it is safe to create a different plan
	 ***/
};

class ThresholdImage_GLRT_FFT : public ThresholdImage {
public:
	ThresholdImage_GLRT_FFT(double PFA_param, double width_param) {PFA = PFA_param; gaussianWidth = width_param; averageKernelXSize = 0; averageKernelYSize = 0; GaussianKernelXSize = 0; GaussianKernelYSize = 0;} 
	~ThresholdImage_GLRT_FFT() {;}
	
	boost::shared_ptr<ublas::matrix <unsigned char> > do_thresholding();
	boost::shared_ptr<ublas::matrix <unsigned char> > do_thresholding(boost::shared_ptr<ublas::matrix<double> > image);
	
protected:
	ConvolveMatricesWithFFTClass matrixConvolver;
	
	double PFA;
	double gaussianWidth;
	boost::shared_ptr<ublas::matrix<double> > Gaussian_kernel;
	boost::shared_ptr<ublas::matrix<double> > average_kernel;
	size_t averageKernelXSize, averageKernelYSize;
	size_t GaussianKernelXSize, GaussianKernelYSize;
	double sum_squared_Gaussian;
	
	boost::mutex GaussianKernelMutex;
	boost::mutex AverageKernelMutex;
	
	boost::shared_mutex gaussianCalculationMutex;
	boost::shared_mutex averageCalculationMutex;
	
	// for an explanation of the many mutexes see ConvolveMatricesWithFFTClass
};

class ThresholdImage_Preprocessor {
public:
	ThresholdImage_Preprocessor() {;}
	virtual ~ThresholdImage_Preprocessor() {;}
	
	virtual boost::shared_ptr<ublas::matrix<double> > do_preprocessing(boost::shared_ptr<ublas::matrix<double> > image) = 0;
};

class ThresholdImage_Preprocessor_MedianFilter : public ThresholdImage_Preprocessor {	// WARNING: we assume that both x and y are odd
public:
	ThresholdImage_Preprocessor_MedianFilter(unsigned x, size_t y) {kernel_x_size = x; kernel_y_size = y;}
	~ThresholdImage_Preprocessor_MedianFilter() {;}
	
	boost::shared_ptr<ublas::matrix<double> > do_preprocessing(boost::shared_ptr<ublas::matrix<double> > image);
	
protected:
	size_t kernel_x_size, kernel_y_size;
};


class ThresholdImage_Preprocessor_GaussianSmoothing : public ThresholdImage_Preprocessor {
public:
	ThresholdImage_Preprocessor_GaussianSmoothing(double Gaussian_width) {width = Gaussian_width;}	// we calculate the convolution kernel once,
	// when the first image is supplied
	~ThresholdImage_Preprocessor_GaussianSmoothing() {;}
	
	boost::shared_ptr<ublas::matrix<double> > do_preprocessing(boost::shared_ptr<ublas::matrix<double> > image);
	
protected:
	void generate_Gaussian_kernel(size_t x_size, size_t y_size);
	ConvolveMatricesWithFFTClass matrixConvolver;
	
	boost::mutex generateKernelMutex;
	
	boost::shared_ptr<ublas::matrix<double> > Gaussian_kernel;	// automatically initialized to a NULL pointer
	double width;
	size_t kernel_x_size;
	size_t kernel_y_size;
};


class ThresholdImage_Preprocessor_MeanFilter : public ThresholdImage_Preprocessor {	// WARNING: we assume that both x and y are odd
public:
	ThresholdImage_Preprocessor_MeanFilter(unsigned x, size_t y) {kernel_x_size = x; kernel_y_size = y;}
	~ThresholdImage_Preprocessor_MeanFilter() {;}
	
	boost::shared_ptr<ublas::matrix<double> > do_preprocessing(boost::shared_ptr<ublas::matrix<double> > image);
	
protected:
	size_t kernel_x_size, kernel_y_size;
};


class ThresholdImage_Postprocessor {
public:
	ThresholdImage_Postprocessor() {;}
	virtual ~ThresholdImage_Postprocessor() {;}
	
	virtual boost::shared_ptr<ublas::matrix <unsigned char> > do_postprocessing(boost::shared_ptr<ublas::matrix <unsigned char> > thresholded_image, boost::shared_ptr<ublas::matrix<double> > image) = 0;
};

class ThresholdImage_Postprocessor_RemoveIsolatedPixels : public ThresholdImage_Postprocessor {	// remove pixels that are considered to be 'on'
	// but do not have any active neighbours
public:
	ThresholdImage_Postprocessor_RemoveIsolatedPixels() {;}
	~ThresholdImage_Postprocessor_RemoveIsolatedPixels() {;}
	
	boost::shared_ptr<ublas::matrix <unsigned char> > do_postprocessing(boost::shared_ptr<ublas::matrix <unsigned char> > thresholded_image, boost::shared_ptr<ublas::matrix<double> > image);
};

gsl_histogram * make_histogram_from_matrix(boost::shared_ptr<ublas::matrix<double> > image, size_t number_of_bins);

#endif

