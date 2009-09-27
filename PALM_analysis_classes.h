/*
 *  PALM_analysis_classes.h
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 25/01/2008.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef PALM_ANALYSIS_CLASSES
#define PALM_ANALYSIS_CLASSES

#include <vector>
#include <queue>
#include <list>
#include <string>
#include <stdexcept>
#include <algorithm>
#include "boost/smart_ptr.hpp"
#include "boost/thread/mutex.hpp"
#include "boost/thread/shared_mutex.hpp"
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_blas.h>
#include "PALM_analysis_defines.h"
#include "PALM_analysis_errors.h"
#include "PALM_analysis_storage.h"
#include "PALM_analysis_FileIO.h"
#include <fftw3.h>
#include "tiffio.h"
#include "stdint.h"
#include "XOPStandardHeaders.h"

#define GSL_RANGE_CHECK_OFF	// this is not required since PALMMatrix<double> does range checks



using namespace std;


/**** format description: the format in which the data is returned from the fitting functions is as follows: ******/
// each position has 11 entries, based on the Gaussian fitting
// these entries are:
// AMPLITUDE	WIDTH		X		Y		OFFSET		AMP_ERROR	WIDTH_ERROR		X_ERROR		Y_ERROR		OFFSET_ERROR		# OF ITERATIONS


/**
 *@brief A class that handles the actual PALM analysis
 *
 *
 * The 'PALMAnalysisController' class must be initialized with pointers to the thresholder, particle finder, fitter, and output writer objects.
 * It will coordinate the PALM analysis between those objects and take care of all the details
 */
class PALMAnalysisController {
public:
	PALMAnalysisController(boost::shared_ptr<ImageLoader> imageLoader_rhs, boost::shared_ptr<ThresholdImage> thresholder_rhs,
						   boost::shared_ptr<ThresholdImage_Preprocessor> thresholdImagePreprocessor_rhs,
						   boost::shared_ptr<ThresholdImage_Postprocessor> thresholdImagePostprocessor_rhs,
						   boost::shared_ptr<ParticleFinder> particleFinder_rhs, boost::shared_ptr<FitPositions> fitPositions_rhs,
						   boost::shared_ptr<PALMResultsWriter> resultsWriter_rhs);
	~PALMAnalysisController() {;}
	
	void DoPALMAnalysis();
	
protected:
	size_t nImages;
	
	std::queue<size_t> framesToBeProcessed;
	std::list<PALMResults> fittedPositions;
	
	boost::shared_ptr<ImageLoader> imageLoader;
	boost::shared_ptr<ThresholdImage> thresholder;
	boost::shared_ptr<ThresholdImage_Preprocessor> thresholdImagePreprocessor;
	boost::shared_ptr<ThresholdImage_Postprocessor> thresholdImagePostprocessor;
	boost::shared_ptr<ParticleFinder> particleFinder;
	boost::shared_ptr<FitPositions> fitPositions;
	boost::shared_ptr<PALMResultsWriter> resultsWriter;
	
	boost::mutex acquireFrameForProcessingMutex;
	boost::mutex addPALMResultsMutex;
	
	/**
	 * The function that will be handle the image processing in a different thread
	 */
	void ThreadPoolWorker();
	
	/**
	 * Compare two PALMResults to see which one originated from an earlier frame
	 */
	int ComparePALMResults(boost::shared_ptr<PALMResults> result1, boost::shared_ptr<PALMResults> result2) {
		assert(result1->getFrameIndex() != result2->getFrameIndex();
		if (result1->getFrameIndex() < result2->getFrameIndex()) {
			return 1;
		} else {
			return 0;
		}
	}
};

/**
 *@brief Stores the result from fitting a particular frame, remembering both the frame index as well as the fitted positions
 */
class PALMResults {
public:
	PALMResults(size_t frameIndex_rhs, boost::shared_ptr<PALMMatrix<double> > fittedPositions_rhs) {
		frameIndex = frameIndex_rhs;
		fittedPositions = fittedPositions_rhs;
	}
	~PALMResults() {;}
	
	size_t getFrameIndex() const {return frameIndex;}
	size_t getNumberOfFittedPositions const {return fittedPositions->ySize();}
	boost::shared_ptr<PALMMatrix<double> > getFittedPositions() const {return fittedPositions;}
	
protected:
	size_t frameIndex;
	boost::shared_ptr<PALMMatrix<double> > fittedPositions;
}

/**
 * @brief An abstract base class from which all other CCD processor classes must derive
 * 
 * The 'CCDImagesProcessor' family of classes contains algorithms that accept a CCD image stack as input (encapsulated in a 'ImageLoader' object)
 * and produce another image stack as output.
 *
 * This class is abstract, meaning that it can not be instantiated. However, it provides the virtual function 'convert_images()', 
 * which is where derived classes should do their work. This function does not accept arguments, instead the parameters a processor needs
 * should be passed in the constructor.
 */
class CCDImagesProcessor {
public:
	CCDImagesProcessor() {;}
	
	virtual ~CCDImagesProcessor() {;}
	
	virtual int convert_images() = 0;
protected:
	size_t total_number_of_images;
	size_t x_size;
	size_t y_size;
	
	ImageLoader *image_loader;
	OutputWriter *output_writer;
};

/**
 * @brief Calculate the average image of all frames and subtract it from every frame in the CCD image stack
 */
class CCDImagesProcessorAverageSubtraction : public CCDImagesProcessor {	// subtracts averages or partial averages from the image trace
public:
	CCDImagesProcessorAverageSubtraction(ImageLoader *i_loader, OutputWriter *o_writer, size_t nFramesAveraging);
	~CCDImagesProcessorAverageSubtraction() {output_writer->flush_and_close();}

	int convert_images();
	
	size_t get_n_frames_averaging() const {return n_frames_averaging;}
	
protected:
	size_t n_frames_averaging;	// how many frames do we average over?
	
	void subtract_average_of_entire_trace();
};

/**
 * @brief Calculates the average image of all frames and subtract it from every frame in the CCD image stack.
 */
class CCDImagesProcessorDifferenceImage : public CCDImagesProcessor {
public:
	CCDImagesProcessorDifferenceImage(ImageLoader *i_loader, OutputWriter *o_writer);
	~CCDImagesProcessorDifferenceImage() {output_writer->flush_and_close();}
	
	int convert_images();
	
// there are no protected members
};

/**
 * @brief Converts the input CCD image stack into a standard format, currently an uncompressed TIFF.
 */
class CCDImagesProcessorConvertToSimpleFileFormat : public CCDImagesProcessor {
public:
	CCDImagesProcessorConvertToSimpleFileFormat(ImageLoader *i_loader, OutputWriter *o_writer);
	~CCDImagesProcessorConvertToSimpleFileFormat() {output_writer->flush_and_close();}
	
	int convert_images();
	
	// there are no protected members
};

/**
 * @brief Generates a new CCD image stack that is a cropped version of the input image stack.
 * The boundary box for the cropping is specified by the parameters in the constructor.
 */
class CCDImagesProcessorCrop : public CCDImagesProcessor {
public:
	CCDImagesProcessorCrop(ImageLoader *i_loader, OutputWriter *o_writer, size_t startX, size_t endX, size_t startY, size_t endY);
	~CCDImagesProcessorCrop() {output_writer->flush_and_close();}
	
	int convert_images();
	
protected:
	size_t startX, endX, startY, endY;	// the points between which we should crop
	size_t croppedXSize, croppedYSize;
};

/**
 * @brief An abstract base class from which all other 'particle finder' classes must derive
 * 
 * This family of classes contains algorithms that accept a single binary segment images as input and reduce this to a list of (x, y) positions
 * where a localization fit will be attempted.
 *
 * This class is abstract, meaning that it can not be instantiated. However, it provides the virtual function 'findPositions()', 
 * which is where derived classes should do their work. Any necessary parameters should be passed in the constructors of derived classes.
 */
class ParticleFinder {
public:
	ParticleFinder() {minDistanceFromEdge = 0;}
	virtual ~ParticleFinder() {;}
	
	virtual boost::shared_ptr<PALMMatrix<double> > findPositions(boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<PALMMatrix <unsigned char> > threshold_image) = 0;
protected:
	double minDistanceFromEdge;
};

/**
 * @brief Assumes that all active pixels within a given radius belong to the same particle
 */
class ParticleFinder_radius : public ParticleFinder {
public:
	ParticleFinder_radius(double dist, double rhs) {minDistanceFromEdge = dist; radius = rhs;}
	~ParticleFinder_radius() {;}
	
	boost::shared_ptr<PALMMatrix<double> > findPositions(boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<PALMMatrix <unsigned char> > threshold_image);
protected:
	double radius;
};

/**
 * @brief Assumes that all neighbouring active pixels (four-way adjacency) belong to the same particle
 */
class ParticleFinder_adjacent4 : public ParticleFinder {
public:
	ParticleFinder_adjacent4(double rhs) {minDistanceFromEdge = rhs;}
	~ParticleFinder_adjacent4() {;}
	
	boost::shared_ptr<PALMMatrix<double> > findPositions(boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<PALMMatrix <unsigned char> > threshold_image);
	
protected:
	void growParticle(position centerPosition, list<position> &positionsInCurrentParticle, boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<PALMMatrix <unsigned char> > threshold_image, boost::shared_ptr<PALMMatrix<long> > mapped_image);

};

/**
 * @brief Assumes that all neighbouring active pixels (eight-way adjacency) belong to the same particle
 */
class ParticleFinder_adjacent8 : public ParticleFinder {
public:
	ParticleFinder_adjacent8(double rhs) {minDistanceFromEdge = rhs;}
	~ParticleFinder_adjacent8() {;}
	
	boost::shared_ptr<PALMMatrix<double> > findPositions(boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<PALMMatrix <unsigned char> > threshold_image);
	
protected:
	void growParticle(position centerPosition, list<position> &positionsInCurrentParticle, boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<PALMMatrix <unsigned char> > threshold_image, boost::shared_ptr<PALMMatrix<long> > mapped_image);
	
};

/**
 * @brief An abstract base class from which all other 'FitPositions' classes must derive
 * 
 * This family of classes contains algorithms that accept a single CCD image and a list of positions as input and return a list of subpixel-localized positions.
 *
 * This class is abstract, meaning that it can not be instantiated. However, it provides the virtual function 'fit_positions()', 
 * which is where derived classes should do their work. Any necessary parameters should be passed in the constructors of derived classes.
 */
class FitPositions {
public:
	FitPositions() {;}
	virtual ~FitPositions() {;}
	
	boost::shared_ptr<PALMMatrix<double> > fit_positions(const boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<PALMMatrix<double> > positions);
	virtual boost::shared_ptr<PALMMatrix<double> > fit_positions(const boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<PALMMatrix<double> > positions, size_t startPos, size_t endPos) = 0;
	// the second function fits the positions between startPos and endPos (indices in the array passed in positions)
	// it's mainly provided to help with multithreading
	
protected:
	/**
	 *@brief Accept raw fitted positions and remove those that were not successful (signaled by a negative amplitude)
	 *
	 * Due to using a not-resizable storage (a matrix) for the fit results, we cannot easily delete unsuccessful fits
	 * this function makes a copy of the fitted positions and returns a possibly smaller matrix that does not contain
	 * the unsuccessful fits
	 */
	boost::shared_ptr<PALMMatrix<double> > RemoveUnsuccessfulFits(boost::shared_ptr<PALMMatrix<double> > fittedPositions);
};

/**
 * @brief Localizes the particles using nonlinear least-squares Levenberg-Marquardt fitting.
 */
class FitPositionsGaussian : public FitPositions {
public:
	// FitPositionsGaussian() {sigma = 0.1;}
	// FitPositionsGaussian(double rhs) {sigma = rhs;}
	FitPositionsGaussian(size_t cutoff_radius_rhs, double r_initial_rhs, double sigma_rhs) {cutoff_radius = cutoff_radius_rhs; r_initial = r_initial_rhs; sigma = sigma_rhs;}
	~FitPositionsGaussian() {;}
	
	boost::shared_ptr<PALMMatrix<double> > fit_positions(const boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<PALMMatrix<double> > positions, size_t startPos, size_t endPos);
	
protected:
	double sigma;
	size_t cutoff_radius;
	double r_initial;
};

class FitPositionsGaussian_FixedWidth : public FitPositions {
public:
	FitPositionsGaussian_FixedWidth(size_t cutoff_radius_rhs, double r_initial_rhs, double sigma_rhs) {cutoff_radius = cutoff_radius_rhs; r_initial = r_initial_rhs; sigma = sigma_rhs;}
	~FitPositionsGaussian_FixedWidth() {;}
	
	boost::shared_ptr<PALMMatrix<double> > fit_positions(const boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<PALMMatrix<double> > positions, size_t startPos, size_t endPos);
	
protected:
	double sigma;
	size_t cutoff_radius;
	double r_initial;
};

class FitPositionsMultiplication : public FitPositions {
	// fits the positions by doing an interative multiplication of the data with a Gaussian at the current best-guess position
	// if this converges then we assume that we have found the actual position
	
	// a description is given in Thompson, Biophys J 82:2775 2002
public:
	FitPositionsMultiplication(size_t cutoff_radius_rhs, double r_initial_rhs, double convergence_rhs) {cutoff_radius = cutoff_radius_rhs; r_initial = r_initial_rhs; convergence_threshold = convergence_rhs;}
	~FitPositionsMultiplication() {;}
	
	// r_initial should be the standard deviation of the Gaussian
	boost::shared_ptr<PALMMatrix<double> > fit_positions(const boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<PALMMatrix<double> > positions, size_t startPos, size_t endPos);
	
protected:
	int multiply_with_gaussian(boost::shared_ptr<PALMMatrix<double> > original_image, boost::shared_ptr<PALMMatrix<double> > masked_image, double x, double y,
							   double std_dev, double background, double amplitude);
	// masked_image should be provided with the same dimensions as original_image. It will be overwritten with the contents of the multiplication
	int determine_x_y_position(boost::shared_ptr<PALMMatrix<double> > masked_image, double &x, double &y);
	
	double convergence_threshold;
	size_t cutoff_radius;
	double r_initial;
};

class FitPositionsCentroid : public FitPositions {
	// fits the positions by calculating a centroid for the pixel values
	
public:
	FitPositionsCentroid(size_t cutoff_radius_rhs) {cutoff_radius = cutoff_radius_rhs;}
	~FitPositionsCentroid() {;}
	
	// r_initial should be the standard deviation of the Gaussian
	boost::shared_ptr<PALMMatrix<double> > fit_positions(const boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<PALMMatrix<double> > positions, size_t startPos, size_t endPos);
	
protected:
	size_t cutoff_radius;
};


class ThresholdImage {
public:
	ThresholdImage() {;}
	virtual ~ThresholdImage() {;}
	
	virtual boost::shared_ptr<PALMMatrix <unsigned char> > do_thresholding(boost::shared_ptr<PALMMatrix<double> > image) = 0;
};

class ThresholdImage_Direct : public ThresholdImage {
	// direct thresholding, based on an absolute threshold
public:
	ThresholdImage_Direct(double thresh_parameter) {threshold = thresh_parameter;}
	~ThresholdImage_Direct() {;}
	
	boost::shared_ptr<PALMMatrix <unsigned char> > do_thresholding(boost::shared_ptr<PALMMatrix<double> > image);
	
protected:
	double threshold;
};

class ThresholdImage_Igor_Iterative : public ThresholdImage {
	// making use of a built-in Igor operation
public:
	ThresholdImage_Igor_Iterative() {;}
	~ThresholdImage_Igor_Iterative() {;}
	
	boost::shared_ptr<PALMMatrix <unsigned char> > do_thresholding(boost::shared_ptr<PALMMatrix<double> > image);
	
protected:
	boost::mutex threadMutex;
};

class ThresholdImage_Igor_Bimodal : public ThresholdImage {
	// making use of a built-in Igor operation
public:
	ThresholdImage_Igor_Bimodal() {;}
	~ThresholdImage_Igor_Bimodal() {;}
	
	boost::shared_ptr<PALMMatrix <unsigned char> > do_thresholding(boost::shared_ptr<PALMMatrix<double> > image);
	
protected:
	boost::mutex threadMutex;
};

class ThresholdImage_Igor_Adaptive : public ThresholdImage {
	// making use of a built-in Igor operation
public:
	ThresholdImage_Igor_Adaptive() {;}
	~ThresholdImage_Igor_Adaptive() {;}
	
	boost::shared_ptr<PALMMatrix <unsigned char> > do_thresholding(boost::shared_ptr<PALMMatrix<double> > image);
	
protected:
	boost::mutex threadMutex;
};

class ThresholdImage_Igor_Fuzzy1 : public ThresholdImage {
	// making use of a built-in Igor operation
public:
	ThresholdImage_Igor_Fuzzy1() {;}
	~ThresholdImage_Igor_Fuzzy1() {;}
	
	boost::shared_ptr<PALMMatrix <unsigned char> > do_thresholding(boost::shared_ptr<PALMMatrix<double> > image);
	
protected:
	boost::mutex threadMutex;
};

class ThresholdImage_Igor_Fuzzy2 : public ThresholdImage {
	// making use of a built-in Igor operation
public:
	ThresholdImage_Igor_Fuzzy2() {;}
	~ThresholdImage_Igor_Fuzzy2() {;}
	
	boost::shared_ptr<PALMMatrix <unsigned char> > do_thresholding(boost::shared_ptr<PALMMatrix<double> > image);
	
protected:
	boost::mutex threadMutex;
};

class ThresholdImage_Isodata : public ThresholdImage {
public:
	ThresholdImage_Isodata() {;}
	~ThresholdImage_Isodata() {;}
	
	boost::shared_ptr<PALMMatrix <unsigned char> > do_thresholding(boost::shared_ptr<PALMMatrix<double> > image);
};

class ThresholdImage_Triangle : public ThresholdImage {
public:
	ThresholdImage_Triangle() {;}
	
	boost::shared_ptr<PALMMatrix <unsigned char> > do_thresholding(boost::shared_ptr<PALMMatrix<double> > image);
};

class ThresholdImage_GLRT : public ThresholdImage {
public:
	ThresholdImage_GLRT(double PFA_param, double width_param) {PFA = PFA_param; gaussianWidth = width_param;} 
	
	boost::shared_ptr<PALMMatrix <unsigned char> > do_thresholding(boost::shared_ptr<PALMMatrix<double> > image);
	
protected:
	double PFA;
	double gaussianWidth;
};

class ConvolveMatricesWithFFTClass {
public:
	ConvolveMatricesWithFFTClass() {forwardPlan = NULL; reversePlan = NULL; forwardPlanXSize = 0; forwardPlanYSize = 0; reversePlanXSize = 0; reversePlanYSize = 0;}
	~ConvolveMatricesWithFFTClass();
	
	boost::shared_ptr<PALMMatrix<double> > ConvolveMatricesWithFFT(boost::shared_ptr<PALMMatrix<double> > image1, boost::shared_ptr<PALMMatrix<double> > image2);
	
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
	
	boost::shared_ptr<PALMMatrix <unsigned char> > do_thresholding();
	boost::shared_ptr<PALMMatrix <unsigned char> > do_thresholding(boost::shared_ptr<PALMMatrix<double> > image);
	
protected:
	ConvolveMatricesWithFFTClass matrixConvolver;
	
	double PFA;
	double gaussianWidth;
	boost::shared_ptr<PALMMatrix<double> > Gaussian_kernel;
	boost::shared_ptr<PALMMatrix<double> > average_kernel;
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
	
	virtual boost::shared_ptr<PALMMatrix<double> > do_preprocessing(boost::shared_ptr<PALMMatrix<double> > image) = 0;
};

class ThresholdImage_Preprocessor_MedianFilter : public ThresholdImage_Preprocessor {	// WARNING: we assume that both x and y are odd
public:
	ThresholdImage_Preprocessor_MedianFilter(unsigned x, size_t y) {kernel_x_size = x; kernel_y_size = y;}
	~ThresholdImage_Preprocessor_MedianFilter() {;}
	
	boost::shared_ptr<PALMMatrix<double> > do_preprocessing(boost::shared_ptr<PALMMatrix<double> > image);
	
protected:
	size_t kernel_x_size, kernel_y_size;
};


class ThresholdImage_Preprocessor_GaussianSmoothing : public ThresholdImage_Preprocessor {
public:
	ThresholdImage_Preprocessor_GaussianSmoothing(double Gaussian_width) {width = Gaussian_width;}	// we calculate the convolution kernel once,
																									// when the first image is supplied
	~ThresholdImage_Preprocessor_GaussianSmoothing() {;}
	
	boost::shared_ptr<PALMMatrix<double> > do_preprocessing(boost::shared_ptr<PALMMatrix<double> > image);
	
protected:
	void generate_Gaussian_kernel(size_t x_size, size_t y_size);
	ConvolveMatricesWithFFTClass matrixConvolver;
	
	boost::mutex generateKernelMutex;
	
	boost::shared_ptr<PALMMatrix<double> > Gaussian_kernel;	// automatically initialized to a NULL pointer
	double width;
	size_t kernel_x_size;
	size_t kernel_y_size;
};


class ThresholdImage_Preprocessor_MeanFilter : public ThresholdImage_Preprocessor {	// WARNING: we assume that both x and y are odd
public:
	ThresholdImage_Preprocessor_MeanFilter(unsigned x, size_t y) {kernel_x_size = x; kernel_y_size = y;}
	~ThresholdImage_Preprocessor_MeanFilter() {;}
	
	boost::shared_ptr<PALMMatrix<double> > do_preprocessing(boost::shared_ptr<PALMMatrix<double> > image);
	
protected:
	size_t kernel_x_size, kernel_y_size;
};


class ThresholdImage_Postprocessor {
public:
	ThresholdImage_Postprocessor() {;}
	virtual ~ThresholdImage_Postprocessor() {;}
	
	virtual boost::shared_ptr<PALMMatrix <unsigned char> > do_postprocessing(boost::shared_ptr<PALMMatrix <unsigned char> > thresholded_image, boost::shared_ptr<PALMMatrix<double> > image) = 0;
};

class ThresholdImage_Postprocessor_RemoveIsolatedPixels : public ThresholdImage_Postprocessor {	// remove pixels that are considered to be 'on'
																								// but do not have any active neighbours
public:
	ThresholdImage_Postprocessor_RemoveIsolatedPixels() {;}
	~ThresholdImage_Postprocessor_RemoveIsolatedPixels() {;}
	
	boost::shared_ptr<PALMMatrix <unsigned char> > do_postprocessing(boost::shared_ptr<PALMMatrix <unsigned char> > thresholded_image, boost::shared_ptr<PALMMatrix<double> > image);
};

class ThresholdImage_Postprocessor_RemovePixelsBelowMean : public ThresholdImage_Postprocessor {	// remove pixels that are considered to be 'on'
																									// but that are below the mean in the image
public:
	ThresholdImage_Postprocessor_RemovePixelsBelowMean() {;}
	~ThresholdImage_Postprocessor_RemovePixelsBelowMean() {;}
	
	boost::shared_ptr<PALMMatrix <unsigned char> > do_postprocessing(boost::shared_ptr<PALMMatrix <unsigned char> > thresholded_image, boost::shared_ptr<PALMMatrix<double> > image);
};


// the routines below are used in the least-squares fitting of a Gaussian to the spots

int Gauss_2D_fit_function(const gsl_vector *params, void *measured_intensities, gsl_vector *deviations);
int Gauss_2D_fit_function_FixedWidth(const gsl_vector *params, void *measured_intensities, gsl_vector *deviations);

int Gauss_2D_fit_function_Jacobian(const gsl_vector *params, void *measured_intensities, gsl_matrix *jacobian);
int Gauss_2D_fit_function_FixedWidth(const gsl_vector *params, void *measured_intensities, gsl_matrix *jacobian);

int Gauss_2D_fit_function_and_Jacobian(const gsl_vector *params, void *measured_intensities_struct, gsl_vector *model_values, gsl_matrix *jacobian);
int Gauss_2D_fit_function_and_Jacobian_FixedWidth(const gsl_vector *params, void *measured_intensities_struct, gsl_vector *model_values, gsl_matrix *jacobian);

// the routines below are simply adapted versions of the least-squares routines above, but have been 'tweaked' to approximate Poissonian instead of Gaussian error distributions

// int Gauss_2D_Poissonian_fit_function(const gsl_vector *params, void *measured_intensities, gsl_vector *model_values);

// int Gauss_2D_Poissonian_fit_function_Jacobian(const gsl_vector *params, void *measured_intensities, gsl_matrix *jacobian);

// int Gauss_2D_Poissonian_fit_function_and_Jacobian(const gsl_vector *params, void *measured_intensities_struct, gsl_vector *model_values, gsl_matrix *jacobian);

class measured_data_Gauss_fits {
public:	
	measured_data_Gauss_fits() {;}
	~measured_data_Gauss_fits() {;}
	
	double xOffset;
	double yOffset;
	double sigma;
	double width;
	boost::shared_ptr<PALMMatrix<double> > imageSubset;
};



#endif
