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
#include "PALM_analysis_classes.h"
#include "PALM_analysis_defines.h"

#define GSL_RANGE_CHECK_OFF	// this is not required since encap_gsl_matrix does range checks

using namespace std;

int load_partial_ccd_image(ImageLoader *image_loader, unsigned long n_start, unsigned long n_end);

int parse_ccd_headers(ImageLoader *image_loader);

int do_analyze_images_operation(boost::shared_ptr<ImageLoader> image_loader, const string output_wave_name, boost::shared_ptr<FitPositions> positions_fitter, 
								boost::shared_ptr<ParticleFinder> particle_finder, boost::shared_ptr<ThresholdImage_Preprocessor> preprocessor, 
								boost::shared_ptr<ThresholdImage> thresholder, boost::shared_ptr<ThresholdImage_Postprocessor> postprocessor);

int do_analyze_images_operation_parallel(boost::shared_ptr<ImageLoader> image_loader, const string output_wave_name, boost::shared_ptr<FitPositions> positions_fitter, 
										 boost::shared_ptr<ParticleFinder> particle_finder, boost::shared_ptr<ThresholdImage_Preprocessor> preprocessor, 
										 boost::shared_ptr<ThresholdImage> thresholder, boost::shared_ptr<ThresholdImage_Postprocessor> postprocessor);

class threadStartDataOld {
public:
	threadStartDataOld (boost::shared_ptr<encap_gsl_matrix> im, boost::shared_ptr<FitPositions> posFitter, boost::shared_ptr<encap_gsl_matrix> pos, boost::shared_ptr<encap_gsl_matrix> localPositions,
					unsigned long startPos, unsigned long endPos) 
					{image = im; positionsFitter = posFitter, positions = pos; localPositionsMatrix = localPositions, startPosition = startPos; endPosition = endPos;}
	~threadStartDataOld() {;}
	
	boost::shared_ptr<encap_gsl_matrix> getPositions() {return positions;}
	boost::shared_ptr<encap_gsl_matrix> getImage() {return image;}
	boost::shared_ptr<FitPositions> getPositionsFitter() {return positionsFitter;}
	boost::shared_ptr<encap_gsl_matrix> getLocalPositionsMatrix() {return localPositionsMatrix;}
	unsigned long getStartPosition() {return startPosition;}
	unsigned long getEndPosition() {return endPosition;}
	
protected:
	unsigned long startPosition;
	unsigned long endPosition;
	boost::shared_ptr<encap_gsl_matrix> positions;
	boost::shared_ptr<encap_gsl_matrix> image;
	boost::shared_ptr<FitPositions> positionsFitter;
	boost::shared_ptr<encap_gsl_matrix> localPositionsMatrix;
};

class threadStartData {
public:
	threadStartData(boost::shared_ptr<ThresholdImage> thresholder_rhs, boost::shared_ptr<ThresholdImage_Preprocessor> preprocessor_rhs,
					 boost::shared_ptr<ThresholdImage_Postprocessor> postprocessor_rhs,
					 boost::shared_ptr<ParticleFinder> particleFinder_rhs, boost::shared_ptr<FitPositions> positionsFitter_rhs) {
		thresholder = thresholder_rhs; preprocessor = preprocessor_rhs; postprocessor = postprocessor_rhs; particleFinder = particleFinder_rhs, positionsFitter = positionsFitter_rhs;
	}
	
	~threadStartData() {;}
	
	boost::shared_ptr<encap_gsl_matrix> image;
	boost::shared_ptr<ThresholdImage> thresholder;
	boost::shared_ptr<ThresholdImage_Preprocessor> preprocessor;
	boost::shared_ptr<ThresholdImage_Postprocessor> postprocessor;
	boost::shared_ptr<ParticleFinder> particleFinder;
	boost::shared_ptr<FitPositions> positionsFitter;
	boost::shared_ptr<encap_gsl_matrix> fittedPositions;	// the positions will be returned here
};
	

void fitPositionsThreadStart(boost::shared_ptr<threadStartData> data);


int do_analyze_images_operation_no_positions_finding(boost::shared_ptr<ImageLoader> image_loader, const string output_wave_name, waveHndl fitting_positions, 
													 boost::shared_ptr<FitPositions> positions_fitter);
// this function is the same as the previous one, except that it does not try to localize the positions before fitting, but assumes that the positions
// to start fitting at are provided as a 2D wave in fitting_positions

boost::shared_ptr<encap_gsl_matrix_uchar> do_processing_and_thresholding(boost::shared_ptr<encap_gsl_matrix> image, boost::shared_ptr<ThresholdImage_Preprocessor>preprocessor, 
																		 boost::shared_ptr<ThresholdImage> thresholder, boost::shared_ptr<ThresholdImage_Postprocessor> postprocessor);

class END_SHOULD_BE_LARGER_THAN_START{};
class GET_NTH_IMAGE_RETURNED_NULL{};


class PALMBitmapImageDeviationCalculator {
public:
	PALMBitmapImageDeviationCalculator();
	~PALMBitmapImageDeviationCalculator();
	
	virtual double getDeviation(boost::shared_ptr<encap_gsl_matrix> positions, size_t index) = 0;
};

class PALMBitmapImageDeviationCalculator_FitUncertainty : public PALMBitmapImageDeviationCalculator {
public:
	PALMBitmapImageDeviationCalculator_FitUncertainty(double scaleFactor_rhs, double upperLimit_rhs) {scaleFactor = scaleFactor_rhs; upperLimit = upperLimit_rhs;}
	~PALMBitmapImageDeviationCalculator_FitUncertainty() {;}
	
	double getDeviation(boost::shared_ptr<encap_gsl_matrix> positions, size_t index) {return ((positions->get(index, 7) < upperLimit) ? (positions->get(index, 7) * scaleFactor) : -1);}
	
private:
	double scaleFactor;
	double upperLimit;
};

class PALMBitmapImageDeviationCalculator_Constant : public PALMBitmapImageDeviationCalculator {
public:
	PALMBitmapImageDeviationCalculator_Constant(double deviation_rhs) {deviation = deviation_rhs;}
	PALMBitmapImageDeviationCalculator_Constant() {;}
	
	double getDeviation(boost::shared_ptr<encap_gsl_matrix> positions, size_t index) {return deviation;}
	
private:
	double deviation;
};

class PALMBitmapImageDeviationCalculator_AmplitudeSquareRoot : public PALMBitmapImageDeviationCalculator {
public:
	PALMBitmapImageDeviationCalculator_AmplitudeSquareRoot(double scaleFactor_rhs) {scaleFactor = scaleFactor_rhs;}
	PALMBitmapImageDeviationCalculator_AmplitudeSquareRoot() {;}
	
	double getDeviation(boost::shared_ptr<encap_gsl_matrix> positions, size_t index) {return scaleFactor * sqrt(positions->get(index, 7));}
	
private:
	double scaleFactor;
};


int construct_summed_intensity_trace(ImageLoader *image_loader, string output_wave_name, long startX, long startY, long endX, long endY);

int construct_average_image(ImageLoader *image_loader, string output_wave_name, long startX, long startY, long endX, long endY);

void calculateStandardDeviationImage(ImageLoader *image_loader, string output_wave_name, long startX, long startY, long endX, long endY);

gsl_histogram * make_histogram_from_matrix(boost::shared_ptr<encap_gsl_matrix> image, unsigned long number_of_bins);

#endif