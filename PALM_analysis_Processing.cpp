/*
 *  PALM_analysis_classes.cpp
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 25/01/2008.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */


#include "PALM_analysis_Processing.h"

void CCDImagesProcessorAverageSubtraction::convert_images() {
	// this function exists so that we could select between a partial or global average
	// for now only global averaging is supported
	this->n_frames_averaging = this->total_number_of_images;
	subtract_average_of_entire_trace();
}

void CCDImagesProcessorAverageSubtraction::subtract_average_of_entire_trace() {
	size_t n;
	boost::shared_ptr<ublas::matrix<double> > average_image;
	boost::shared_ptr<ublas::matrix<double> > loaded_image;
	boost::shared_ptr<ublas::matrix<double> > subtracted_image;
	
	// we pass through the images two times:
	// the first pass calculates the average,
	// the second pass subtracts it from the image
	// fortunately out intermediate format uses doubles to store the data!
	
	average_image = boost::shared_ptr<ublas::matrix<double> >(new ublas::matrix<double>(x_size, y_size));
	subtracted_image = boost::shared_ptr<ublas::matrix<double> >(new ublas::matrix<double>(x_size, y_size));
	
	std::fill(average_image->data().begin(), average_image->data().end(), double(0.0));
	
	for (n = 0; n < total_number_of_images; n++) {
		loaded_image = image_loader->get_nth_image(n);
		
		(*average_image) += (*loaded_image);
	}
	
	// now divide each point so that we get the average
	*average_image /= (double)this->n_frames_averaging;
	
	// now subtract the average for each frame
	for (n = 0; n < total_number_of_images; n++) {
		loaded_image = image_loader->get_nth_image(n);
		
		ublas::noalias(*subtracted_image) = (*loaded_image) - (*average_image);
		
		output_writer->write_image(subtracted_image);
	}
}

void CCDImagesProcessorDifferenceImage::convert_images() {
	
	boost::shared_ptr<ublas::matrix<double> > current_image;
	boost::shared_ptr<ublas::matrix<double> > next_image;
	
	// we start by loading the first image in next_image
	// this is required so the loop that follows can start properly
	
	next_image = image_loader->get_nth_image(0);
	
	
	for (size_t n = 0; n < (total_number_of_images - 1); n++) {
		
		// the previous image for this run of the loop is the image that was previously in current_image
		// so we have to shift it down
		current_image = next_image;
		
		next_image = image_loader->get_nth_image(n + 1);
		
		// now do the actual subtraction
		*current_image = *current_image - *next_image;
		
		// current_image now contains the subtracted image, we should write it to disk
		output_writer->write_image(current_image);
	}
}

void CCDImagesProcessorConvertToSimpleFileFormat::convert_images() {
	
	boost::shared_ptr<ublas::matrix<double> > current_image;
	
	for (size_t n = 0; n < total_number_of_images; ++n) {
		current_image = image_loader->get_nth_image(n);
		output_writer->write_image(current_image);
	}
}

void CCDImagesProcessorCrop::convert_images() {
	boost::shared_ptr<ublas::matrix <double> > croppedImage;
	boost::shared_ptr<ublas::matrix <double> > loadedImage;
	
	for (size_t n = 0; n < total_number_of_images; ++n) {
		loadedImage = image_loader->get_nth_image(n);
		croppedImage = boost::shared_ptr<ublas::matrix<double> > (new ublas::matrix<double> (croppedXSize, croppedYSize));
		
		for (size_t x = startX; x <= endX; ++x) {
			for (size_t y = startY; y <= endY; ++y) {
				(*croppedImage)(x - startX, y - startY) = (*loadedImage)(x, y);
			}
		}
		
		output_writer->write_image(croppedImage);
	}
}

void CCDImagesProcessorConvertToPhotons::convert_images() {
	assert(this->multiplicationFactor > 0.0);
	
	this->total_number_of_images = image_loader->get_total_number_of_images();
	this->x_size = image_loader->getXSize();
	this->y_size = image_loader->getYSize();
	
	boost::shared_ptr<ublas::matrix <double> > loadedImage;
	boost::shared_ptr<ublas::matrix <double> > convertedImage;
	
	for (size_t n = 0; n < this->total_number_of_images; ++n) {
		loadedImage = image_loader->get_nth_image(n);
		for (size_t i = 0; i < this->x_size; ++i) {
			for (size_t j = 0; j < this->y_size; ++j) {
				(*convertedImage)(i, j) = ((*loadedImage)(i, j) - this->offset >= 0) ? (((*loadedImage)(i, j) - this->offset) / this->multiplicationFactor) : 0.0;
			}
		}
		output_writer->write_image(convertedImage);
	}
}

