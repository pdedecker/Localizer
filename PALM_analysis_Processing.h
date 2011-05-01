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

#define GSL_RANGE_CHECK_OFF	// this is not required since Eigen::MatrixXddoes range checks

#include <vector>
#include <queue>
#include <deque>
#include <list>
#include <string>
#include <algorithm>
#include <boost/smart_ptr.hpp>
#include <Eigen/Eigen>
#include "PALM_analysis_defines.h"
#include "PALM_analysis_errors.h"
#include "PALM_analysis_storage.h"
#include "PALM_analysis_FileIO.h"

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
	CCDImagesProcessor(boost::shared_ptr<ProgressReporter> progressReporter_rhs) : progressReporter(progressReporter_rhs) {}
	virtual ~CCDImagesProcessor() {;}
	
	virtual void convert_images(boost::shared_ptr<ImageLoader> image_loader, boost::shared_ptr<ImageOutputWriter> output_writer) = 0;
	
protected:
	boost::shared_ptr<ProgressReporter> progressReporter;
};

/**
 * @brief Calculate the average image of all frames and subtract it from every frame in the CCD image stack
 */
class CCDImagesProcessorAverageSubtraction : public CCDImagesProcessor {	// subtracts averages or partial averages from the image trace
public:
	CCDImagesProcessorAverageSubtraction(boost::shared_ptr<ProgressReporter> progressReporter_rhs, size_t nFramesAveraging_rhs) : CCDImagesProcessor(progressReporter_rhs), n_frames_averaging(nFramesAveraging_rhs) {}
	~CCDImagesProcessorAverageSubtraction() {;}

	void convert_images(boost::shared_ptr<ImageLoader> image_loader, boost::shared_ptr<ImageOutputWriter> output_writer);
	
	size_t get_n_frames_averaging() const {return n_frames_averaging;}
	
protected:
	size_t n_frames_averaging;	// how many frames do we average over?
	
	void subtractAverageOfEntireMovie(boost::shared_ptr<ImageLoader> image_loader, boost::shared_ptr<ImageOutputWriter> output_writer);
	void subtractRollingAverage(boost::shared_ptr<ImageLoader> image_loader, boost::shared_ptr<ImageOutputWriter> output_writer, size_t nFramesInAverage);
};

/**
 * @brief Calculates the average image of all frames and subtract it from every frame in the CCD image stack.
 */
class CCDImagesProcessorDifferenceImage : public CCDImagesProcessor {
public:
	CCDImagesProcessorDifferenceImage(boost::shared_ptr<ProgressReporter> progressReporter_rhs) : CCDImagesProcessor(progressReporter_rhs) {}
	~CCDImagesProcessorDifferenceImage() {;}
	
	void convert_images(boost::shared_ptr<ImageLoader> image_loader, boost::shared_ptr<ImageOutputWriter> output_writer);
	
// there are no protected members
};

/**
 * @brief Converts the input CCD image stack into a standard format, currently an uncompressed TIFF.
 */
class CCDImagesProcessorConvertToSimpleFileFormat : public CCDImagesProcessor {
public:
	CCDImagesProcessorConvertToSimpleFileFormat(boost::shared_ptr<ProgressReporter> progressReporter_rhs) : CCDImagesProcessor(progressReporter_rhs) {}
	~CCDImagesProcessorConvertToSimpleFileFormat() {;}
	
	void convert_images(boost::shared_ptr<ImageLoader> image_loader, boost::shared_ptr<ImageOutputWriter> output_writer);
	
	// there are no protected members
};

/**
 * @brief Generates a new CCD image stack that is a cropped version of the input image stack.
 * The boundary box for the cropping is specified by the parameters in the constructor.
 */
class CCDImagesProcessorCrop : public CCDImagesProcessor {
public:
	CCDImagesProcessorCrop(boost::shared_ptr<ProgressReporter> progressReporter_rhs, size_t startX_rhs, size_t endX_rhs, size_t startY_rhs, size_t endY_rhs) : CCDImagesProcessor(progressReporter_rhs), startX(startX_rhs), endX(endX_rhs), startY(startY_rhs), endY(endY_rhs) {}
	~CCDImagesProcessorCrop() {;}
	
	void convert_images(boost::shared_ptr<ImageLoader> image_loader, boost::shared_ptr<ImageOutputWriter> output_writer);
	
protected:
	size_t startX, endX, startY, endY;	// the points between which we should crop
	size_t croppedXSize, croppedYSize;
};

/**
 * @brief Tries to convert a recorded CCD movie to represent the intensity in number of photons.
 */
class CCDImagesProcessorConvertToPhotons : public CCDImagesProcessor {
public:
	CCDImagesProcessorConvertToPhotons(boost::shared_ptr<ProgressReporter> progressReporter_rhs, double multiplicationFactor_rhs, double offset_rhs) : CCDImagesProcessor(progressReporter_rhs), multiplicationFactor(multiplicationFactor_rhs), offset(offset_rhs) {}
	~CCDImagesProcessorConvertToPhotons() {;}
	
	void convert_images(boost::shared_ptr<ImageLoader> image_loader, boost::shared_ptr<ImageOutputWriter> output_writer);
	
protected:
	double multiplicationFactor;
	double offset;
};


#endif
