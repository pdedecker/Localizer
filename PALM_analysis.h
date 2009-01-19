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

int do_analyze_images_operation_parallel2(boost::shared_ptr<ImageLoader> image_loader, const string output_wave_name, boost::shared_ptr<FitPositions> positions_fitter, 
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
	threadStartData(boost::shared_ptr<encap_gsl_matrix> image_rhs, boost::shared_ptr<FitPositions> positionsFitter_rhs, boost::shared_ptr<encap_gsl_matrix> positions_rhs, 
					boost::shared_ptr<encap_gsl_matrix> outputPositions_rhs, unsigned long start_rhs, unsigned long end_rhs) {image = image_rhs; positions = positions_rhs; outputPositions = outputPositions_rhs; positionsFitter = positionsFitter_rhs; start = start_rhs; end = end_rhs;}
	~threadStartData() {;}
	
	boost::shared_ptr<encap_gsl_matrix> image;
	boost::shared_ptr<encap_gsl_matrix> positions;
	boost:: shared_ptr<encap_gsl_matrix> outputPositions;
	boost::shared_ptr<FitPositions> positionsFitter;
	
	unsigned long start;
	unsigned long end;
};

void fitPositionsThreadStart(threadStartData data);

class threadStartData2 {
public:
	threadStartData2(boost::shared_ptr<ThresholdImage> thresholder_rhs, boost::shared_ptr<ThresholdImage_Preprocessor> preprocessor_rhs,
					 boost::shared_ptr<ThresholdImage_Postprocessor> postprocessor_rhs,
					 boost::shared_ptr<ParticleFinder> particleFinder_rhs, boost::shared_ptr<FitPositions> positionsFitter_rhs) {
		thresholder = thresholder_rhs; preprocessor = preprocessor_rhs; postprocessor = postprocessor_rhs; particleFinder = particleFinder_rhs, positionsFitter = positionsFitter_rhs;
	}
	
	~threadStartData2() {;}
	
	boost::shared_ptr<encap_gsl_matrix> image;
	boost::shared_ptr<ThresholdImage> thresholder;
	boost::shared_ptr<ThresholdImage_Preprocessor> preprocessor;
	boost::shared_ptr<ThresholdImage_Postprocessor> postprocessor;
	boost::shared_ptr<ParticleFinder> particleFinder;
	boost::shared_ptr<FitPositions> positionsFitter;
	boost::shared_ptr<encap_gsl_matrix> fittedPositions;	// the positions will be returned here
};
	

void fitPositionsThreadStart2(boost::shared_ptr<threadStartData2> data);


int do_analyze_images_operation_no_positions_finding(boost::shared_ptr<ImageLoader> image_loader, const string output_wave_name, waveHndl fitting_positions, 
													 boost::shared_ptr<FitPositions> positions_fitter);
// this function is the same as the previous one, except that it does not try to localize the positions before fitting, but assumes that the positions
// to start fitting at are provided as a 2D wave in fitting_positions

boost::shared_ptr<encap_gsl_matrix_uchar> do_processing_and_thresholding(boost::shared_ptr<encap_gsl_matrix> image, boost::shared_ptr<ThresholdImage_Preprocessor>preprocessor, 
																		 boost::shared_ptr<ThresholdImage> thresholder, boost::shared_ptr<ThresholdImage_Postprocessor> postprocessor);

class END_SHOULD_BE_LARGER_THAN_START{};
class GET_NTH_IMAGE_RETURNED_NULL{};

struct thread_start_parameters {
	gsl_matrix *image;
	unsigned long n_image;
	int analysis_method;	// pass 1 for 2D Gaussian fitting
	double treshold;
	unsigned long radius;
	unsigned long min_distance_from_edge;
	unsigned long cutoff_radius;
	double r_initial;
	double sigma;
	double background;
	struct thread_return_parameters *return_params;	// we provide this so we don't get problems with the return going out-of-scope when the thread finishes
};

struct thread_return_parameters {
	gsl_matrix *positions;
	unsigned long n_image;
};

/* gsl_matrix * find_positions(const gsl_matrix * image, double treshold, unsigned long radius, unsigned long min_distance_from_edge);
	// the return values from this function is a N * 3 matrix with intensity and x and y coordinates in the columns

// gsl_matrix * fit_Gauss_2D_to_positions(const gsl_matrix *image, gsl_matrix *positions, unsigned long cutoff_radius, double r_initial, double background, double sigma);
	// WARNING: cutoff should never be larger than radius in the previous function!

gsl_matrix * find_positions_no_sort(const gsl_matrix * image, gsl_matrix_uchar *thresholded_image, unsigned long radius, unsigned long min_distance_from_edge); */

int construct_summed_intensity_trace(ImageLoader *image_loader, string output_wave_name, long startX, long startY, long endX, long endY);

int construct_average_image(ImageLoader *image_loader, string output_wave_name, long startX, long startY, long endX, long endY);

void calculateStandardDeviationImage(ImageLoader *image_loader, string output_wave_name, long startX, long startY, long endX, long endY);

gsl_histogram * make_histogram_from_matrix(boost::shared_ptr<encap_gsl_matrix> image, unsigned long number_of_bins);

void print_localized_positions(gsl_matrix *positions);
	
void print_fitted_positions(gsl_matrix *positions);

void print_image_to_file(gsl_matrix *image, string file_path);

void print_fit_state(unsigned long iteration, gsl_multifit_fdfsolver *s);

void print_initial_parameters(gsl_vector *params);

void print_image(gsl_matrix * image);

#endif