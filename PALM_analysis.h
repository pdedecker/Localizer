#ifndef PALM_ANALYSIS
#define PALM_ANALYSIS

#include <iostream>
#include <fstream>
#include <vector>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_histogram.h>
#include "boost/thread.hpp"
#include "boost/bind.hpp"
#include "XOPStandardHeaders.h"
#include "PALM_analysis_errors.h"
#include "PALM_analysis_classes.h"
#include "PALM_analysis_defines.h"
#include "PALM_analysis_FileIO.h"

#define GSL_RANGE_CHECK_OFF	// this is not required since PALMMatrix<double> does range checks

using namespace std;

int load_partial_ccd_image(ImageLoader *image_loader, size_t n_start, size_t n_end);

int parse_ccd_headers(ImageLoader *image_loader);


boost::shared_ptr<PALMMatrix <unsigned char> > do_processing_and_thresholding(boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<ThresholdImage_Preprocessor>preprocessor, 
																		 boost::shared_ptr<ThresholdImage> thresholder, boost::shared_ptr<ThresholdImage_Postprocessor> postprocessor);

class END_SHOULD_BE_LARGER_THAN_START{};
class GET_NTH_IMAGE_RETURNED_NULL{};


class PALMBitmapImageDeviationCalculator {
public:
	PALMBitmapImageDeviationCalculator() {;}
	virtual ~PALMBitmapImageDeviationCalculator() {;}
	
	virtual double getDeviation(boost::shared_ptr<PALMMatrix<double> > positions, size_t index) = 0;
};

class PALMBitmapImageDeviationCalculator_FitUncertainty : public PALMBitmapImageDeviationCalculator {
public:
	PALMBitmapImageDeviationCalculator_FitUncertainty(double scaleFactor_rhs, double upperLimit_rhs) {scaleFactor = scaleFactor_rhs; upperLimit = upperLimit_rhs;}
	~PALMBitmapImageDeviationCalculator_FitUncertainty() {;}
	
	double getDeviation(boost::shared_ptr<PALMMatrix<double> > positions, size_t index) {return ((positions->get(index, 7) < upperLimit) ? (positions->get(index, 7) * scaleFactor) : -1);}
	
private:
	double scaleFactor;
	double upperLimit;
};

class PALMBitmapImageDeviationCalculator_Constant : public PALMBitmapImageDeviationCalculator {
public:
	PALMBitmapImageDeviationCalculator_Constant(double deviation_rhs) {deviation = deviation_rhs;}
	PALMBitmapImageDeviationCalculator_Constant() {;}
	
	double getDeviation(boost::shared_ptr<PALMMatrix<double> > positions, size_t index) {return deviation;}
	
private:
	double deviation;
};

class PALMBitmapImageDeviationCalculator_AmplitudeSquareRoot : public PALMBitmapImageDeviationCalculator {
public:
	PALMBitmapImageDeviationCalculator_AmplitudeSquareRoot(double PSFWidth_rhs, double scaleFactor_rhs) {PSFWidth = PSFWidth_rhs; scaleFactor = scaleFactor_rhs;}
	PALMBitmapImageDeviationCalculator_AmplitudeSquareRoot() {;}
	
	double getDeviation(boost::shared_ptr<PALMMatrix<double> > positions, size_t index) {return PSFWidth / (scaleFactor * sqrt(positions->get(index, 1)));}
	
private:
	double PSFWidth;
	double scaleFactor;
};

class calculate_PALM_bitmap_image_ThreadStartParameters {
public:
	calculate_PALM_bitmap_image_ThreadStartParameters() {;}
	~calculate_PALM_bitmap_image_ThreadStartParameters() {;}
	
	boost::shared_ptr<PALMMatrix<double> > positions;
	boost::shared_ptr<PALMVolume <unsigned short> > image;
	boost::shared_ptr<PALMMatrix<double> > totalIntensities;
	
	boost::shared_ptr<PALMMatrix<double> > colors;
	boost::shared_ptr<PALMBitmapImageDeviationCalculator> deviationCalculator;
	
	int normalizeColors;
	size_t nFrames;
	size_t startIndex;
	size_t endIndex;
	size_t imageWidth;
	size_t imageHeight;
	size_t xSize;
	size_t ySize;
	double maxAmplitude;	// maximum amplitude of a fitted Gaussian over the entire fitted positions
	double scaleFactor;
};

boost::shared_ptr<PALMVolume <unsigned short> > calculate_PALM_bitmap_image(boost::shared_ptr<PALMMatrix<double> > positions, boost::shared_ptr<PALMMatrix<double> > colors, boost::shared_ptr<PALMBitmapImageDeviationCalculator> deviationCalculator,
																size_t xSize, size_t ySize, size_t imageWidth, size_t imageHeight, int normalizeColors);

boost::shared_ptr<PALMVolume <unsigned short> > calculate_PALM_bitmap_image_parallel(boost::shared_ptr<PALMMatrix<double> > positions, boost::shared_ptr<PALMMatrix<double> > colors, boost::shared_ptr<PALMBitmapImageDeviationCalculator> deviationCalculator,
																		 size_t xSize, size_t ySize, size_t imageWidth, size_t imageHeight, int normalizeColors);

void calculate_PALM_bitmap_image_ThreadStart(boost::shared_ptr<calculate_PALM_bitmap_image_ThreadStartParameters> startParameters);


int construct_summed_intensity_trace(ImageLoader *image_loader, string output_wave_name, long startX, long startY, long endX, long endY);

int construct_average_image(ImageLoader *image_loader, string output_wave_name, long startX, long startY, long endX, long endY);

void calculateStandardDeviationImage(ImageLoader *image_loader, string output_wave_name, long startX, long startY, long endX, long endY);

gsl_histogram * make_histogram_from_matrix(boost::shared_ptr<PALMMatrix<double> > image, size_t number_of_bins);

#endif
