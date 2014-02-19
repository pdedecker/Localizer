/*
 Copyright 2008-2014 Peter Dedecker.
 
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

#include "PALM_analysis_defines.h"
#include "PALM_analysis_errors.h"
#include "PALM_analysis_storage.h"

class ImageLoader;
class ImageOutputWriter;
class ProgressReporter;


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
	CCDImagesProcessor(std::shared_ptr<ProgressReporter> progressReporter_rhs) : progressReporter(progressReporter_rhs) {}
	virtual ~CCDImagesProcessor() {;}
	
	virtual void convert_images(std::shared_ptr<ImageLoader> image_loader, std::shared_ptr<ImageOutputWriter> output_writer) = 0;
	
protected:
	std::shared_ptr<ProgressReporter> progressReporter;
};

/**
 * @brief Calculate the average image of all frames and subtract it from every frame in the CCD image stack
 */
class CCDImagesProcessorAverageSubtraction : public CCDImagesProcessor {	// subtracts averages or partial averages from the image trace
public:
	CCDImagesProcessorAverageSubtraction(std::shared_ptr<ProgressReporter> progressReporter_rhs, size_t nFramesAveraging_rhs) : CCDImagesProcessor(progressReporter_rhs), n_frames_averaging(nFramesAveraging_rhs) {}
	~CCDImagesProcessorAverageSubtraction() {;}

	void convert_images(std::shared_ptr<ImageLoader> image_loader, std::shared_ptr<ImageOutputWriter> output_writer);
	
	size_t get_n_frames_averaging() const {return n_frames_averaging;}
	
protected:
	size_t n_frames_averaging;	// how many frames do we average over?
	
	void subtractAverageOfEntireMovie(std::shared_ptr<ImageLoader> image_loader, std::shared_ptr<ImageOutputWriter> output_writer);
	void subtractRollingAverage(std::shared_ptr<ImageLoader> image_loader, std::shared_ptr<ImageOutputWriter> output_writer, int nFramesInAverage);
};

/**
 * @brief Calculates the average image of all frames and subtract it from every frame in the CCD image stack.
 */
class CCDImagesProcessorDifferenceImage : public CCDImagesProcessor {
public:
	CCDImagesProcessorDifferenceImage(std::shared_ptr<ProgressReporter> progressReporter_rhs) : CCDImagesProcessor(progressReporter_rhs) {}
	~CCDImagesProcessorDifferenceImage() {;}
	
	void convert_images(std::shared_ptr<ImageLoader> image_loader, std::shared_ptr<ImageOutputWriter> output_writer);
	
// there are no protected members
};

/**
 * @brief Converts the input CCD image stack into a standard format, currently an uncompressed TIFF.
 */
class CCDImagesProcessorConvertToSimpleFileFormat : public CCDImagesProcessor {
public:
	CCDImagesProcessorConvertToSimpleFileFormat(std::shared_ptr<ProgressReporter> progressReporter_rhs) : CCDImagesProcessor(progressReporter_rhs) {}
	~CCDImagesProcessorConvertToSimpleFileFormat() {;}
	
	void convert_images(std::shared_ptr<ImageLoader> image_loader, std::shared_ptr<ImageOutputWriter> output_writer);
	
	// there are no protected members
};

/**
 * @brief Generates a new CCD image stack that is a cropped version of the input image stack.
 * The boundary box for the cropping is specified by the parameters in the constructor.
 */
class CCDImagesProcessorCrop : public CCDImagesProcessor {
public:
	CCDImagesProcessorCrop(std::shared_ptr<ProgressReporter> progressReporter_rhs, size_t startX_rhs, size_t endX_rhs, size_t startY_rhs, size_t endY_rhs) : CCDImagesProcessor(progressReporter_rhs), startX(startX_rhs), endX(endX_rhs), startY(startY_rhs), endY(endY_rhs) {}
	~CCDImagesProcessorCrop() {;}
	
	void convert_images(std::shared_ptr<ImageLoader> image_loader, std::shared_ptr<ImageOutputWriter> output_writer);
	
protected:
	int startX, endX, startY, endY;	// the points between which we should crop
	int croppedXSize, croppedYSize;
};

/**
 * @brief Tries to convert a recorded CCD movie to represent the intensity in number of photons.
 */
class CCDImagesProcessorConvertToPhotons : public CCDImagesProcessor {
public:
	CCDImagesProcessorConvertToPhotons(std::shared_ptr<ProgressReporter> progressReporter_rhs, double multiplicationFactor_rhs, double offset_rhs) : CCDImagesProcessor(progressReporter_rhs), multiplicationFactor(multiplicationFactor_rhs), offset(offset_rhs) {}
	~CCDImagesProcessorConvertToPhotons() {;}
	
	void convert_images(std::shared_ptr<ImageLoader> image_loader, std::shared_ptr<ImageOutputWriter> output_writer);
	
protected:
	double multiplicationFactor;
	double offset;
};


#endif
