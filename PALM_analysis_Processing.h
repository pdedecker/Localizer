/*
 *  PALM_analysis_classes.h
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 25/01/2008.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef PALM_ANALYSIS_PROCESSING_H
#define PALM_ANALYSIS_PROCESSING_H

#define GSL_RANGE_CHECK_OFF	// this is not required since ublas::matrix<double> does range checks

#include <vector>
#include <queue>
#include <list>
#include <string>
#include <algorithm>
#include <boost/smart_ptr.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "PALM_analysis_defines.h"
#include "PALM_analysis_errors.h"
#include "PALM_analysis_storage.h"
#include "PALM_analysis_FileIO.h"

namespace ublas = boost::numeric::ublas;

class ImageLoader;
class ImageOutputWriter;


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
	CCDImagesProcessor(boost::shared_ptr<ImageLoader> imageLoader_rhs, boost::shared_ptr<ImageOutputWriter> outputWriter_rhs) :
		image_loader(imageLoader_rhs), output_writer(outputWriter_rhs), total_number_of_images(0), x_size(0), y_size(0) {}
	
	virtual ~CCDImagesProcessor() {;}
	
	virtual void convert_images() = 0;
protected:
	size_t total_number_of_images;
	size_t x_size;
	size_t y_size;
	
	boost::shared_ptr<ImageLoader> image_loader;
	boost::shared_ptr<ImageOutputWriter> output_writer;
};

/**
 * @brief Calculate the average image of all frames and subtract it from every frame in the CCD image stack
 */
class CCDImagesProcessorAverageSubtraction : public CCDImagesProcessor {	// subtracts averages or partial averages from the image trace
public:
	CCDImagesProcessorAverageSubtraction(boost::shared_ptr<ImageLoader> imageLoader_rhs, boost::shared_ptr<ImageOutputWriter> outputWriter_rhs, size_t nFramesAveraging_rhs) :
		CCDImagesProcessor(imageLoader_rhs, outputWriter_rhs), n_frames_averaging(nFramesAveraging_rhs) {}
	
	~CCDImagesProcessorAverageSubtraction() {;}

	void convert_images();
	
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
	CCDImagesProcessorDifferenceImage(boost::shared_ptr<ImageLoader> imageLoader_rhs, boost::shared_ptr<ImageOutputWriter> outputWriter_rhs) :
		CCDImagesProcessor(imageLoader_rhs, outputWriter_rhs) {}
	
	~CCDImagesProcessorDifferenceImage() {;}
	
	void convert_images();
	
// there are no protected members
};

/**
 * @brief Converts the input CCD image stack into a standard format, currently an uncompressed TIFF.
 */
class CCDImagesProcessorConvertToSimpleFileFormat : public CCDImagesProcessor {
public:
	CCDImagesProcessorConvertToSimpleFileFormat(boost::shared_ptr<ImageLoader> imageLoader_rhs, boost::shared_ptr<ImageOutputWriter> outputWriter_rhs) :
		CCDImagesProcessor(imageLoader_rhs, outputWriter_rhs) {}
	~CCDImagesProcessorConvertToSimpleFileFormat() {;}
	
	
	void convert_images();
	
	// there are no protected members
};

/**
 * @brief Generates a new CCD image stack that is a cropped version of the input image stack.
 * The boundary box for the cropping is specified by the parameters in the constructor.
 */
class CCDImagesProcessorCrop : public CCDImagesProcessor {
public:
	CCDImagesProcessorCrop(boost::shared_ptr<ImageLoader> imageLoader_rhs, boost::shared_ptr<ImageOutputWriter> outputWriter_rhs,
						   size_t startX_rhs, size_t endX_rhs, size_t startY_rhs, size_t endY_rhs) :
		CCDImagesProcessor(imageLoader_rhs, outputWriter_rhs), startX(startX_rhs), endX(endX_rhs), startY(startY_rhs), endY(endY_rhs) {
			if ((startX >= endX) || (startY >= endY)) {
				throw std::runtime_error("Bad limit values specified for cropping");
			}
			
			if ((endX >= x_size) || (endY >= y_size)) {
				throw std::runtime_error("Bad limit values specified for cropping");
			}
			
			croppedXSize = endX - startX + 1;
			croppedYSize = endY - startY + 1;
		}
	
	~CCDImagesProcessorCrop() {;}
	
	void convert_images();
	
protected:
	size_t startX, endX, startY, endY;	// the points between which we should crop
	size_t croppedXSize, croppedYSize;
};

/**
 * @brief Tries to convert a recorded CCD movie to represent the intensity in number of photons.
 */
class CCDImagesProcessorConvertToPhotons : public CCDImagesProcessor {
public:
	CCDImagesProcessorConvertToPhotons(boost::shared_ptr<ImageLoader> imageLoader_rhs, boost::shared_ptr<ImageOutputWriter> outputWriter_rhs,
									   double multiplicationFactor_rhs, double offset_rhs) :
		CCDImagesProcessor(imageLoader_rhs, outputWriter_rhs), multiplicationFactor(multiplicationFactor_rhs), offset(offset_rhs) {}
	~CCDImagesProcessorConvertToPhotons() {;}
	
	void convert_images();
	
protected:
	double multiplicationFactor;
	double offset;
};


#endif
