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

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <queue>
#include <list>
#include <string>
#include <stdexcept>
#include "boost/smart_ptr.hpp"
#include "boost/thread/mutex.hpp"
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_blas.h>
#include "PALM_analysis_defines.h"
#include <fftw3.h>
#include "tiffio.h"
#include "stdint.h"
#include "XOPStandardHeaders.h"

#define GSL_RANGE_CHECK_OFF	// this is not require since encap_gsl_matrix does range checks



using namespace std;


/**** format description: the format in which the data is returned from the fitting functions is as follows: ******/
// each position has 11 entries, based on the Gaussian fitting
// these entries are:
// AMPLITUDE	WIDTH		X		Y		OFFSET		AMP_ERROR	WIDTH_ERROR		X_ERROR		Y_ERROR		OFFSET_ERROR		# OF ITERATIONS


class encap_gsl_matrix {
public:
	encap_gsl_matrix() {matrix = NULL;}
	encap_gsl_matrix(size_t x, size_t y);
	~encap_gsl_matrix();
	
	gsl_matrix* get_ptr() const {return matrix;}
	
	void set(size_t x, size_t y, double value);
	void set_all(double value) {gsl_matrix_set_all(matrix, value);}
	void operator=(gsl_matrix* rhs);
	
	double get(size_t x, size_t y);
	size_t get_x_size() const {return x_size;}
	size_t get_y_size() const {return y_size;}
	
protected:
	string error;
	gsl_matrix *matrix;
	size_t x_size, y_size;
};

class encap_gsl_matrix_uchar {
public:
	encap_gsl_matrix_uchar() {matrix = NULL;}
	encap_gsl_matrix_uchar(size_t x, size_t y);
	~encap_gsl_matrix_uchar();
	
	gsl_matrix_uchar* get_ptr() const {return matrix;}
	
	void set(size_t x, size_t y, unsigned char value);
	void set_all(unsigned char value) {gsl_matrix_uchar_set_all(matrix, value);}
	void operator=(gsl_matrix_uchar* rhs);
	
	unsigned char get(size_t x, size_t y);
	size_t get_x_size() const {return x_size;}
	size_t get_y_size() const {return y_size;}
	
protected:
	string error;
	gsl_matrix_uchar *matrix;
	size_t x_size, y_size;
};

class encap_gsl_matrix_long {
public:
	encap_gsl_matrix_long() {matrix = NULL;}
	encap_gsl_matrix_long(size_t x, size_t y);
	~encap_gsl_matrix_long();
	
	gsl_matrix_long* get_ptr() const {return matrix;}
	
	void set(size_t x, size_t y, long value);
	void set_all(long value) {gsl_matrix_long_set_all(matrix, value);}
	void operator=(gsl_matrix_long* rhs);
	
	long get(size_t x, size_t y);
	size_t get_x_size() const {return x_size;}
	size_t get_y_size() const {return y_size;}
	
protected:
	string error;
	gsl_matrix_long *matrix;
	size_t x_size, y_size;
};

class position {
public:
	position() {x = 0; y = 0; intensity = 0;}
	position(double xLoc, double yLoc) {x = xLoc; y = yLoc;}
	position(double xLoc, double yLoc, double intensity_rhs) {x = xLoc; y = yLoc; intensity = intensity_rhs;}
	~position() {;}
	
	void set_x(double xLoc) {x = xLoc;}
	void set_y(double yLoc) {y = yLoc;}
	void set_intensity(double intensity_rhs) {intensity = intensity_rhs;}
	
	double get_x() {return x;}
	double get_y() {return y;}
	double get_intensity() {return intensity;}
	
protected:
	double x;
	double y;
	double intensity;
};


class ImageLoader {
public:
	ImageLoader();
	ImageLoader(const string rhs);
	virtual ~ImageLoader();
	
	// int open_file(const string rhs);
	// int close_current_file();
	
	unsigned long get_total_number_of_images() const {return total_number_of_images;}
	unsigned long get_x_size() const {return x_size;}
	unsigned long get_y_size() const {return y_size;}
	virtual boost::shared_ptr<encap_gsl_matrix> get_nth_image(const unsigned long n) = 0;	// images are numbered from 0 to N - 1
	
protected:
	virtual int parse_header_information() = 0;
	
	unsigned long image_cache_size;
	
	string path;
	ifstream file;
	unsigned long header_length;
	unsigned long total_number_of_images;
	unsigned long x_size;
	unsigned long y_size;
	int storage_type;
	vector <boost::shared_ptr<encap_gsl_matrix> > image_cache;	// array of pointer to gsl_matrices containing cached images
	vector <unsigned long> images_in_cache;	// keeps track of the indices of the images in the cache
	// if there is no image at a particular location then this is set to -1
};

class ImageLoaderAndor : public ImageLoader {
public:
	ImageLoaderAndor(string rhs);
	ImageLoaderAndor(string rhs, unsigned long image_cache_size_rhs);
	~ImageLoaderAndor();
	
	boost::shared_ptr<encap_gsl_matrix> get_nth_image(const unsigned long n);
	
protected:
	int parse_header_information();
};

class ImageLoaderHamamatsu : public ImageLoader {
public:
	ImageLoaderHamamatsu(string rhs);
	ImageLoaderHamamatsu(string rhs, unsigned long image_cache_size);
	~ImageLoaderHamamatsu();
	
	boost::shared_ptr<encap_gsl_matrix> get_nth_image(const unsigned long n);
	
protected:
	int parse_header_information();
};

class ImageLoaderSPE : public ImageLoader {
public:
	ImageLoaderSPE(string rhs);
	ImageLoaderSPE(string rhs, unsigned long image_cache_size);
	~ImageLoaderSPE();
	
	boost::shared_ptr<encap_gsl_matrix> get_nth_image(const unsigned long n);
	
protected:
	int parse_header_information();
};

class SimpleImageLoader : public ImageLoader {	// loads data from a binary file from a square array consisting of unsigned longs in row-major order
public:
	SimpleImageLoader(string rhs);
	SimpleImageLoader(string rhs, unsigned long image_cache_size);
	~SimpleImageLoader();
	
	boost::shared_ptr<encap_gsl_matrix> get_nth_image(const unsigned long n);
	
protected:
	int parse_header_information();
};

class ImageLoaderTIFF_Igor : public ImageLoader {	// loads data from TIFF files using extensive callbacks to IGOR because the TIFF format is a mess
public:
	ImageLoaderTIFF_Igor(string rhs);
	ImageLoaderTIFF_Igor(string rhs, unsigned long image_cache_size);
	~ImageLoaderTIFF_Igor();
	
	boost::shared_ptr<encap_gsl_matrix> get_nth_image(const unsigned long n);
	
protected:
	int parse_header_information();
	
	string native_Igor_FilePath;	// in Windows format or the native HFS Macintosh format
	waveHndl tiff_cache_wave;
	
};

class ImageLoaderTIFF : public ImageLoader {	// loads data from TIFF files using the libtiff library
public:
	ImageLoaderTIFF(string rhs);
	ImageLoaderTIFF(string rhs, unsigned long image_cache_size);
	~ImageLoaderTIFF();
	
	boost::shared_ptr<encap_gsl_matrix> get_nth_image(const unsigned long n);
	
protected:
	int parse_header_information();
	
	TIFF* tiff_file;
	unsigned int bitsPerPixel;
	unsigned int sampleFormat;	// unsigned integer of floating point?
								// 1 for uint, 3 for floating point (same as tiff specification)
	vector<unsigned long> directoryIndices;	// some images in the TIFF format are not actual data, but low resolution previews of other images and such
											// this means that the image at directory i is not necessarily the i-th experimental frame
											// we correct for this using this array
	
};

class ImageLoaderIgor : public ImageLoader {
public:
	ImageLoaderIgor(string waveName);
	~ImageLoaderIgor() {;}
	
	boost::shared_ptr<encap_gsl_matrix> get_nth_image(const unsigned long n);
	
protected:
	int parse_header_information() {return 0;}
	
	waveHndl igor_data_wave;
};



class OutputWriter {
public:
	OutputWriter();
	OutputWriter(const string &rhs, int overwrite);
	virtual ~OutputWriter() {;}
	
	string get_file_path() const {return file_path;}
	unsigned long get_n_images_written() const {return n_images_written;}
	
	virtual void write_image(boost::shared_ptr<encap_gsl_matrix> new_image) = 0;
	
	virtual int flush_and_close() = 0;
	
protected:
	
	string file_path;
	ofstream file;
	
	unsigned long n_images_written;
	unsigned long x_size;
	unsigned long y_size;
	
	queue <boost::shared_ptr<encap_gsl_matrix> > image_buffer;
};


class SimpleOutputWriter : public OutputWriter {
public:
	SimpleOutputWriter(const string &rhs, int overwrite);
	~SimpleOutputWriter();
	
	void write_image(boost::shared_ptr<encap_gsl_matrix> new_image);
	
	int flush_and_close();
	
protected:
	void flush_cache();
};

class TIFFOutputWriter : public OutputWriter {
public:
	TIFFOutputWriter(const string &rhs, int overwrite);
	~TIFFOutputWriter();
	
	void write_image(boost::shared_ptr<encap_gsl_matrix> new_image);
	
	int flush_and_close();
	
protected:
	void flush_cache();
	int compress;	// if 0 then don't compress the data, otherwise compress
	
	TIFF *tiff_file;
};

class IgorOutputWriter {	// this class will take care of outputting the fitted positions to an Igor wave
public:
	IgorOutputWriter(const string rhs);
	~IgorOutputWriter();
	
	int append_new_positions(boost::shared_ptr<encap_gsl_matrix> positions);
	int write_positions_to_wave();	// to be called when the positions have been uploaded
	// writes the positions to the wave
	// the wave is only created at this point, and WILL overwrite previous waves of the same name
	
private:
	vector<boost::shared_ptr<encap_gsl_matrix> > positions_array;
	unsigned long total_number_of_positions;
	string wave_name;
};



class CCDImagesProcessor {
public:
	CCDImagesProcessor() {;}
	CCDImagesProcessor(ImageLoader *i_loader, OutputWriter *o_writer) {;}
	
	virtual ~CCDImagesProcessor() {;}
	
	virtual void set_n_parameter(double n) {;}
	
	virtual int convert_images() = 0;
protected:
	unsigned long total_number_of_images;
	unsigned long x_size;
	unsigned long y_size;
	
	ImageLoader *image_loader;
	OutputWriter *output_writer;
};

class CCDImagesProcessorAverageSubtraction : public CCDImagesProcessor {	// subtracts averages or partial averages from the image trace
public:
	CCDImagesProcessorAverageSubtraction(ImageLoader *i_loader, OutputWriter *o_writer);
	~CCDImagesProcessorAverageSubtraction() {output_writer->flush_and_close();}

	int convert_images();
	
	void set_n_parameter(double n);		// in this class this function sets the number of frames that we average over
	unsigned long get_n_frames_averaging() const {return n_frames_averaging;}
	
protected:
	unsigned long n_frames_averaging;	// how many frames do we average over?
	
	void subtract_average_of_entire_trace();
	void subtract_partial_average();
};

class CCDImagesProcessorDifferenceImage : public CCDImagesProcessor {	// creates a difference image, where we subtract the next frame from the current one
																		// similar to the method used in Betzig et al, Science 2006
																		// the result is a trace of length N-1
public:
	CCDImagesProcessorDifferenceImage(ImageLoader *i_loader, OutputWriter *o_writer);
	~CCDImagesProcessorDifferenceImage() {output_writer->flush_and_close();}
	
	int convert_images();
	
	void set_n_parameter(double n) {;}		// ignored in this class
	
// there are no protected members
};

class CCDImagesProcessorConvertToSimpleFileFormat : public CCDImagesProcessor {	// creates a difference image, where we subtract the next frame from the current one
	// similar to the method used in Betzig et al, Science 2006
	// the result is a trace of length N-1
public:
	CCDImagesProcessorConvertToSimpleFileFormat(ImageLoader *i_loader, OutputWriter *o_writer);
	~CCDImagesProcessorConvertToSimpleFileFormat() {output_writer->flush_and_close();}
	
	int convert_images();
	
	void set_n_parameter(double n) {;}		// ignored in this class
	
	// there are no protected members
};


class ParticleFinder {
public:
	ParticleFinder() {minDistanceFromEdge = 0;}
	virtual ~ParticleFinder() {;}
	
	virtual boost::shared_ptr<encap_gsl_matrix> findPositions(boost::shared_ptr<encap_gsl_matrix> image, boost::shared_ptr<encap_gsl_matrix_uchar> threshold_image) = 0;
protected:
	double minDistanceFromEdge;
};

class ParticleFinder_radius : public ParticleFinder {
public:
	ParticleFinder_radius(double dist, double rhs) {minDistanceFromEdge = dist; radius = rhs;}
	~ParticleFinder_radius() {;}
	
	boost::shared_ptr<encap_gsl_matrix> findPositions(boost::shared_ptr<encap_gsl_matrix> image, boost::shared_ptr<encap_gsl_matrix_uchar> threshold_image);
protected:
	double radius;
};

class ParticleFinder_adjacent4 : public ParticleFinder {
public:
	ParticleFinder_adjacent4(double rhs) {minDistanceFromEdge = rhs;}
	~ParticleFinder_adjacent4() {;}
	
	boost::shared_ptr<encap_gsl_matrix> findPositions(boost::shared_ptr<encap_gsl_matrix> image, boost::shared_ptr<encap_gsl_matrix_uchar> threshold_image);
	
protected:
	void growParticle(position centerPosition, list<position> &positionsInCurrentParticle, boost::shared_ptr<encap_gsl_matrix> image, boost::shared_ptr<encap_gsl_matrix_uchar> threshold_image, boost::shared_ptr<encap_gsl_matrix_long> mapped_image);

};

class ParticleFinder_adjacent8 : public ParticleFinder {
public:
	ParticleFinder_adjacent8(double rhs) {minDistanceFromEdge = rhs;}
	~ParticleFinder_adjacent8() {;}
	
	boost::shared_ptr<encap_gsl_matrix> findPositions(boost::shared_ptr<encap_gsl_matrix> image, boost::shared_ptr<encap_gsl_matrix_uchar> threshold_image);
	
protected:
	void growParticle(position centerPosition, list<position> &positionsInCurrentParticle, boost::shared_ptr<encap_gsl_matrix> image, boost::shared_ptr<encap_gsl_matrix_uchar> threshold_image, boost::shared_ptr<encap_gsl_matrix_long> mapped_image);
	
};

class FitPositions {
public:
	FitPositions() {;}
	virtual ~FitPositions() {;}
	
	virtual boost::shared_ptr<encap_gsl_matrix> fit_positions(const boost::shared_ptr<encap_gsl_matrix> image, boost::shared_ptr<encap_gsl_matrix> positions) = 0;
	virtual boost::shared_ptr<encap_gsl_matrix> fit_positions(const boost::shared_ptr<encap_gsl_matrix> image, boost::shared_ptr<encap_gsl_matrix> positions, unsigned long startPos, unsigned long endPos) = 0;
	// the second function fits the positions between startPos and endPos (indices in the array passed in positions)
	// it's mainly provided to help with multithreading
};

class FitPositionsGaussian : public FitPositions {
public:
	// FitPositionsGaussian() {sigma = 0.1;}
	// FitPositionsGaussian(double rhs) {sigma = rhs;}
	FitPositionsGaussian(unsigned long cutoff_radius_rhs, double r_initial_rhs, double background_rhs, double sigma_rhs) {cutoff_radius = cutoff_radius_rhs; r_initial = r_initial_rhs; background = background_rhs; sigma = sigma_rhs;}
	~FitPositionsGaussian() {;}
	
	boost::shared_ptr<encap_gsl_matrix> fit_positions(const boost::shared_ptr<encap_gsl_matrix> image, boost::shared_ptr<encap_gsl_matrix> positions);
	boost::shared_ptr<encap_gsl_matrix> fit_positions(const boost::shared_ptr<encap_gsl_matrix> image, boost::shared_ptr<encap_gsl_matrix> positions, unsigned long startPos, unsigned long endPos);
	
	void set_sigma(double rhs) {sigma = rhs;}
	double get_sigma() const {return sigma;}
	
protected:
	double sigma;
	unsigned long cutoff_radius;
	double background;
	double r_initial;
};

class FitPositionsMultiplication : public FitPositions {
	// fits the positions by doing an interative multiplication of the data with a Gaussian at the current best-guess position
	// if this converges then we assume that we have found the actual position
	
	// a description is given in Thompson, Biophys J 82:2775 2002
public:
	FitPositionsMultiplication(unsigned long cutoff_radius_rhs, double r_initial_rhs, double background_rhs, double convergence_rhs) {cutoff_radius = cutoff_radius_rhs; r_initial = r_initial_rhs; background = background_rhs; convergence_threshold = convergence_rhs;}
	~FitPositionsMultiplication() {;}
	
	boost::shared_ptr<encap_gsl_matrix> fit_positions(const boost::shared_ptr<encap_gsl_matrix> image, boost::shared_ptr<encap_gsl_matrix> positions);
	// r_initial should be the standard deviation of the Gaussian
	boost::shared_ptr<encap_gsl_matrix> fit_positions(const boost::shared_ptr<encap_gsl_matrix> image, boost::shared_ptr<encap_gsl_matrix> positions, unsigned long startPos, unsigned long endPos);
	
protected:
	int multiply_with_gaussian(boost::shared_ptr<encap_gsl_matrix> original_image, boost::shared_ptr<encap_gsl_matrix> masked_image, double x, double y,
							   double std_dev, double background, double amplitude);
	// masked_image should be provided with the same dimensions as original_image. It will be overwritten with the contents of the multiplication
	int determine_x_y_position(boost::shared_ptr<encap_gsl_matrix> masked_image, double &x, double &y);
	
	double convergence_threshold;
	unsigned long cutoff_radius;
	double background;
	double r_initial;
};

class FitPositionsCentroid : public FitPositions {
	// fits the positions by calculating a centroid for the pixel values
	
public:
	FitPositionsCentroid(unsigned long cutoff_radius_rhs) {cutoff_radius = cutoff_radius_rhs;}
	~FitPositionsCentroid() {;}
	
	boost::shared_ptr<encap_gsl_matrix> fit_positions(const boost::shared_ptr<encap_gsl_matrix> image, boost::shared_ptr<encap_gsl_matrix> positions);
	// r_initial should be the standard deviation of the Gaussian
	boost::shared_ptr<encap_gsl_matrix> fit_positions(const boost::shared_ptr<encap_gsl_matrix> image, boost::shared_ptr<encap_gsl_matrix> positions, unsigned long startPos, unsigned long endPos);
	
protected:
	unsigned long cutoff_radius;
};


class ThresholdImage {
public:
	ThresholdImage() {;}
	virtual ~ThresholdImage() {;}
	
	virtual boost::shared_ptr<encap_gsl_matrix_uchar> do_thresholding(boost::shared_ptr<encap_gsl_matrix> image) = 0;
};

class ThresholdImage_Direct : public ThresholdImage {
	// direct thresholding, based on an absolute threshold
public:
	ThresholdImage_Direct(double thresh_parameter) {threshold = thresh_parameter;}
	~ThresholdImage_Direct() {;}
	
	boost::shared_ptr<encap_gsl_matrix_uchar> do_thresholding(boost::shared_ptr<encap_gsl_matrix> image);
	
protected:
	double threshold;
};

class ThresholdImage_Igor_Iterative : public ThresholdImage {
	// making use of a built-in Igor operation
public:
	ThresholdImage_Igor_Iterative() {;}
	~ThresholdImage_Igor_Iterative() {;}
	
	boost::shared_ptr<encap_gsl_matrix_uchar> do_thresholding(boost::shared_ptr<encap_gsl_matrix> image);
	
protected:
	boost::mutex threadMutex;
};

class ThresholdImage_Igor_Bimodal : public ThresholdImage {
	// making use of a built-in Igor operation
public:
	ThresholdImage_Igor_Bimodal() {;}
	~ThresholdImage_Igor_Bimodal() {;}
	
	boost::shared_ptr<encap_gsl_matrix_uchar> do_thresholding(boost::shared_ptr<encap_gsl_matrix> image);
	
protected:
	boost::mutex threadMutex;
};

class ThresholdImage_Igor_Adaptive : public ThresholdImage {
	// making use of a built-in Igor operation
public:
	ThresholdImage_Igor_Adaptive() {;}
	~ThresholdImage_Igor_Adaptive() {;}
	
	boost::shared_ptr<encap_gsl_matrix_uchar> do_thresholding(boost::shared_ptr<encap_gsl_matrix> image);
	
protected:
	boost::mutex threadMutex;
};

class ThresholdImage_Igor_Fuzzy1 : public ThresholdImage {
	// making use of a built-in Igor operation
public:
	ThresholdImage_Igor_Fuzzy1() {;}
	~ThresholdImage_Igor_Fuzzy1() {;}
	
	boost::shared_ptr<encap_gsl_matrix_uchar> do_thresholding(boost::shared_ptr<encap_gsl_matrix> image);
	
protected:
	boost::mutex threadMutex;
};

class ThresholdImage_Igor_Fuzzy2 : public ThresholdImage {
	// making use of a built-in Igor operation
public:
	ThresholdImage_Igor_Fuzzy2() {;}
	~ThresholdImage_Igor_Fuzzy2() {;}
	
	boost::shared_ptr<encap_gsl_matrix_uchar> do_thresholding(boost::shared_ptr<encap_gsl_matrix> image);
	
protected:
	boost::mutex threadMutex;
};

class ThresholdImage_Isodata : public ThresholdImage {
public:
	ThresholdImage_Isodata() {;}
	~ThresholdImage_Isodata() {;}
	
	boost::shared_ptr<encap_gsl_matrix_uchar> do_thresholding(boost::shared_ptr<encap_gsl_matrix> image);
};

class ThresholdImage_Triangle : public ThresholdImage {
public:
	ThresholdImage_Triangle() {;}
	
	boost::shared_ptr<encap_gsl_matrix_uchar> do_thresholding(boost::shared_ptr<encap_gsl_matrix> image);
};

class ThresholdImage_MTT : public ThresholdImage {
public:
	ThresholdImage_MTT(double PFA_param, double width_param) {PFA = PFA_param; gaussianWidth = width_param;} 
	
	boost::shared_ptr<encap_gsl_matrix_uchar> do_thresholding(boost::shared_ptr<encap_gsl_matrix> image);
	
protected:
	double PFA;
	double gaussianWidth;
};

class ConvolveMatricesWithFFTClass {
public:
	ConvolveMatricesWithFFTClass() {forwardPlanExists = 0; reversePlanExists = 0; forwardPlan = NULL; reversePlan = NULL; xSize = 0; ySize = 0;}
	~ConvolveMatricesWithFFTClass();
	
	boost::shared_ptr<encap_gsl_matrix> ConvolveMatricesWithFFT(boost::shared_ptr<encap_gsl_matrix> image1, boost::shared_ptr<encap_gsl_matrix> image2);
	
protected:
	fftw_plan forwardPlan;
	fftw_plan reversePlan;
	
	int forwardPlanExists;
	int reversePlanExists;
	
	unsigned long xSize;
	unsigned long ySize;
	
	boost::mutex planMutex;
};

class ThresholdImage_MTT_FFT : public ThresholdImage {
public:
	ThresholdImage_MTT_FFT(double PFA_param, double width_param) {PFA = PFA_param; gaussianWidth = width_param; kernel_x_size = 0; kernel_y_size = 0;} 
	~ThresholdImage_MTT_FFT() {;}
	
	boost::shared_ptr<encap_gsl_matrix_uchar> do_thresholding();
	boost::shared_ptr<encap_gsl_matrix_uchar> do_thresholding(boost::shared_ptr<encap_gsl_matrix> image);
	
protected:
	ConvolveMatricesWithFFTClass matrixConvolver;
	
	double PFA;
	double gaussianWidth;
	boost::shared_ptr<encap_gsl_matrix> Gaussian_kernel;
	boost::shared_ptr<encap_gsl_matrix> average_kernel;
	unsigned long kernel_x_size, kernel_y_size;
	double sum_squared_Gaussian;
	boost::mutex GaussianKernelMutex;
	boost::mutex AverageKernelMutex;
};

class ThresholdImage_Preprocessor {
public:
	ThresholdImage_Preprocessor() {;}
	virtual ~ThresholdImage_Preprocessor() {;}
	
	virtual boost::shared_ptr<encap_gsl_matrix> do_preprocessing(boost::shared_ptr<encap_gsl_matrix> image) = 0;
};

class ThresholdImage_Preprocessor_MedianFilter : public ThresholdImage_Preprocessor {	// WARNING: we assume that both x and y are odd
public:
	ThresholdImage_Preprocessor_MedianFilter(unsigned x, unsigned long y) {kernel_x_size = x; kernel_y_size = y;}
	~ThresholdImage_Preprocessor_MedianFilter() {;}
	
	boost::shared_ptr<encap_gsl_matrix> do_preprocessing(boost::shared_ptr<encap_gsl_matrix> image);
	
protected:
	unsigned long kernel_x_size, kernel_y_size;
};


class ThresholdImage_Preprocessor_GaussianSmoothing : public ThresholdImage_Preprocessor {
public:
	ThresholdImage_Preprocessor_GaussianSmoothing(double Gaussian_width) {width = Gaussian_width;}	// we calculate the convolution kernel once,
																									// when the first image is supplied
	~ThresholdImage_Preprocessor_GaussianSmoothing() {;}
	
	boost::shared_ptr<encap_gsl_matrix> do_preprocessing(boost::shared_ptr<encap_gsl_matrix> image);
	
protected:
	void generate_Gaussian_kernel(unsigned long x_size, unsigned long y_size);
	ConvolveMatricesWithFFTClass matrixConvolver;
	
	boost::mutex generateKernelMutex;
	
	boost::shared_ptr<encap_gsl_matrix> Gaussian_kernel;	// automatically initialized to a NULL pointer
	double width;
	unsigned long kernel_x_size;
	unsigned long kernel_y_size;
};


class ThresholdImage_Preprocessor_MeanFilter : public ThresholdImage_Preprocessor {	// WARNING: we assume that both x and y are odd
public:
	ThresholdImage_Preprocessor_MeanFilter(unsigned x, unsigned long y) {kernel_x_size = x; kernel_y_size = y;}
	~ThresholdImage_Preprocessor_MeanFilter() {;}
	
	boost::shared_ptr<encap_gsl_matrix> do_preprocessing(boost::shared_ptr<encap_gsl_matrix> image);
	
protected:
	unsigned long kernel_x_size, kernel_y_size;
};


class ThresholdImage_Postprocessor {
public:
	ThresholdImage_Postprocessor() {;}
	virtual ~ThresholdImage_Postprocessor() {;}
	
	virtual boost::shared_ptr<encap_gsl_matrix_uchar> do_postprocessing(boost::shared_ptr<encap_gsl_matrix_uchar> thresholded_image, boost::shared_ptr<encap_gsl_matrix> image) = 0;
};

class ThresholdImage_Postprocessor_RemoveIsolatedPixels : public ThresholdImage_Postprocessor {	// remove pixels that are considered to be 'on'
																								// but do not have any active neighbours
public:
	ThresholdImage_Postprocessor_RemoveIsolatedPixels() {;}
	~ThresholdImage_Postprocessor_RemoveIsolatedPixels() {;}
	
	boost::shared_ptr<encap_gsl_matrix_uchar> do_postprocessing(boost::shared_ptr<encap_gsl_matrix_uchar> thresholded_image, boost::shared_ptr<encap_gsl_matrix> image);
};

class ThresholdImage_Postprocessor_RemovePixelsBelowMean : public ThresholdImage_Postprocessor {	// remove pixels that are considered to be 'on'
																									// but that are below the mean in the image
public:
	ThresholdImage_Postprocessor_RemovePixelsBelowMean() {;}
	~ThresholdImage_Postprocessor_RemovePixelsBelowMean() {;}
	
	boost::shared_ptr<encap_gsl_matrix_uchar> do_postprocessing(boost::shared_ptr<encap_gsl_matrix_uchar> thresholded_image, boost::shared_ptr<encap_gsl_matrix> image);
};


boost::shared_ptr<encap_gsl_matrix> convolve_matrices_using_fft(boost::shared_ptr<encap_gsl_matrix> image1, boost::shared_ptr<encap_gsl_matrix> image2);


// the routines below are used in the least-squares fitting of a Gaussian to the spots

int Gauss_2D_fit_function(const gsl_vector *params, void *measured_intensities, gsl_vector *model_values);

int Gauss_2D_fit_function_Jacobian(const gsl_vector *params, void *measured_intensities, gsl_matrix *jacobian);

int Gauss_2D_fit_function_and_Jacobian(const gsl_vector *params, void *measured_intensities_struct, gsl_vector *model_values, gsl_matrix *jacobian);

// the routines below are simply adapted versions of the least-squares routines above, but have been 'tweaked' to approximate Poissonian instead of Gaussian error distributions

int Gauss_2D_Poissonian_fit_function(const gsl_vector *params, void *measured_intensities, gsl_vector *model_values);

int Gauss_2D_Poissonian_fit_function_Jacobian(const gsl_vector *params, void *measured_intensities, gsl_matrix *jacobian);

int Gauss_2D_Poissonian_fit_function_and_Jacobian(const gsl_vector *params, void *measured_intensities_struct, gsl_vector *model_values, gsl_matrix *jacobian);

double * convert_gsl_matrix_to_array(gsl_matrix *matrix, unsigned long &number_of_elements);

inline int return_xy_from_array(unsigned long pos, unsigned long x_size, unsigned long y_size, double x_offset, double y_offset, double &x, double &y);
// because we have to convert the image data to a linear array for the fitting, we need some way to covert the indices
// for the array to the x,y coordinates in the image
// WARNING: we assume that the data is stored in a row-major order

gsl_vector * convert_gsl_matrix_to_vector(gsl_matrix *matrix);
// ROW MAJOR ordering


/*class OutputWriter {
public:
	OutputWriter(const unsigned long n_images);
	OutputWriter(const string path, const unsigned long n_images);
	~OutputWriter();
	
	int create_file(const string rhs);	// if there is already a file open then -1 will be returned, and nothing will happen
	int write_positions_in_single_image(const gsl_matrix *positions);	// this writes the images to the file incrementally, so the first image should be first,
	// then the second, and so on until the last
	int close_file;		// should be called after the last image has been submitted
	// if the number of images that is submitted does not equal the predicted total number of images then an error will be raised
	
protected:
	unsigned long total_number_of_images;
	unsigned long number_of_images_written;
	unsigned long number_of_positions_written;
	unsigned long *number_of_positions_array;		// this array keeps track of how many positions are found in the m-th image
	// the size of this array is equal to the total number of images
	unsigned long total_number_of_positions;
	string file_path;
	ofstream file;
};*/

class measured_data_Gauss_fits {
public:	
	measured_data_Gauss_fits();
	~measured_data_Gauss_fits();
	unsigned long get_number_of_intensities() {return number_of_intensities;}
	unsigned long get_x_size() {return x_size;}
	unsigned long get_y_size() {return y_size;}
	double get_x_offset() {return x_offset;}
	double get_y_offset() {return y_offset;}
	double * get_intensities() {return intensities;}
	double * get_sigma() {return sigma;}
	void * get_x_y_position(unsigned long pos, double &x, double &y);
	
	void set_number_of_intensities(unsigned long n) {number_of_intensities = n;}
	void set_x_size(unsigned long size) {x_size = size;}
	void set_y_size(unsigned long size) {y_size = size;}
	void set_x_offset(double offset) {x_offset = offset;}
	void set_y_offset(double offset) {y_offset = offset;}
	void set_intensities(double *intens) {intensities = intens;}
	void set_sigma(double *s) {sigma = s;}
	void set_x_y_positions(gsl_matrix *pos) {x_y_positions = pos;}
	
protected:
	unsigned long number_of_intensities;
	unsigned long x_size;
	unsigned long y_size;
	double x_offset;
	double y_offset;
	double *intensities;	// NOTE: we do not delete this data when the object is destroyed
	double *sigma;	// NOTE: we do not delete this data when the object is destroyed
	gsl_matrix *x_y_positions;
};

class OUT_OF_MEMORY {
public:
	OUT_OF_MEMORY(string rhs) {error = rhs;}
	~OUT_OF_MEMORY() {;}
	
	string get_error() {return error;}
	
protected:
	string error;
};
class CANNOT_OPEN_FILE {};
class DIMENSIONS_SHOULD_BE_EQUAL {};
class KERNEL_SIZE_SHOULD_BE_ODD {};
class IMAGE_INDEX_BEYOND_N_IMAGES {};
class GET_NTH_IMAGE_FILE_NOT_OPEN {};
class CANNOT_DETERMINE_SPE_STORAGE_TYPE {};
class CANNOT_OPEN_OUTPUT_FILE {};
class SIZE_OF_CHAR_IS_NOT_ONE_BYTE {};
class SIZE_OF_FLOAT_IS_NOT_FOUR_BYTES {};
class OUTPUT_FILE_ALREADY_EXISTS {};
class NUMBER_OF_AVERAGING_FRAMES_SHOULD_BE_ODD {};
class NUMBER_OF_AVERAGING_LARGER_THAN_N_FRAMES_IN_FILE {};

class ERROR_READING_FILE_DATA {
public:
	ERROR_READING_FILE_DATA(string rhs) {error = rhs;}
	~ERROR_READING_FILE_DATA() {;}
	
	string get_error() {return error;}
	
protected:
	string error;
};

class ERROR_WRITING_FILE_DATA {
public:
	ERROR_WRITING_FILE_DATA(string rhs) {error = rhs;}
	~ERROR_WRITING_FILE_DATA() {;}
	
	string get_error() {return error;}
	
protected:
	string error;
};



#endif