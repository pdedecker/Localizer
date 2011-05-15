/*
 *  PALM_analysis_classes.cpp
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 25/01/2008.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */


#include "PALM_analysis_Processing.h"

void CCDImagesProcessorAverageSubtraction::convert_images(boost::shared_ptr<ImageLoader> image_loader, boost::shared_ptr<ImageOutputWriter> output_writer) {
	// depending on the value of this->n_frames_averaging we either need to subtract the entire movie or just a small part
	this->progressReporter->CalculationStarted();
	
	if (this->n_frames_averaging == 0) {
		this->subtractAverageOfEntireMovie(image_loader, output_writer);
	} else {
		this->subtractRollingAverage(image_loader, output_writer, this->n_frames_averaging);
	}
	
	this->progressReporter->CalculationDone();
}

void CCDImagesProcessorAverageSubtraction::subtractAverageOfEntireMovie(boost::shared_ptr<ImageLoader> image_loader, boost::shared_ptr<ImageOutputWriter> output_writer) {	
	// this function exists so that we could select between a partial or global average
	// for now only global averaging is supported
	
	size_t x_size = image_loader->getXSize();
	size_t y_size = image_loader->getYSize();
	size_t total_number_of_images = image_loader->getNImages();
	boost::shared_ptr<Eigen::MatrixXd> average_image;
	boost::shared_ptr<Eigen::MatrixXd> loaded_image;
	boost::shared_ptr<Eigen::MatrixXd> subtracted_image;
	
	int abortStatus;
	
	// we pass through the images two times:
	// the first pass calculates the average,
	// the second pass subtracts it from the image
	
	if (this->n_frames_averaging == 0)
		this->n_frames_averaging = total_number_of_images;
	
	average_image = boost::shared_ptr<Eigen::MatrixXd>(new Eigen::MatrixXd((int)x_size, (int)y_size));
	subtracted_image = boost::shared_ptr<Eigen::MatrixXd>(new Eigen::MatrixXd((int)x_size, (int)y_size));
	
	average_image->setConstant(0.0);
	
	for (size_t n = 0; n < total_number_of_images; n++) {
		loaded_image = image_loader->readImage(n);
		
		(*average_image) += (*loaded_image);
		
		// test for user abort and provide a progress update every 10 frames
		if (n % 10 == 0) {
			abortStatus = this->progressReporter->UpdateCalculationProgress((double)n / (double)total_number_of_images * 50.0, 100.0);
			if (abortStatus != 0) {
				// user abort
				this->progressReporter->CalculationAborted();
				throw USER_ABORTED("aborted by user request");
			}
		}
	}
	
	// now divide each point so that we get the average
	*average_image /= (double)this->n_frames_averaging;
	
	// now subtract the average for each frame
	for (size_t n = 0; n < total_number_of_images; n++) {
		loaded_image = image_loader->readImage(n);
		
		(*subtracted_image) = (*loaded_image) - (*average_image);
		
		output_writer->write_image(subtracted_image);
		
		// test for user abort and provide a progress update every 10 frames
		if (n % 10 == 0) {
			abortStatus = this->progressReporter->UpdateCalculationProgress((double)n / (double)total_number_of_images * 50.0 + 50.0, 100.0);
			if (abortStatus != 0) {
				// user abort
				this->progressReporter->CalculationAborted();
				throw USER_ABORTED("aborted by user request");
			}
		}
	}
	
}

void CCDImagesProcessorAverageSubtraction::subtractRollingAverage(boost::shared_ptr<ImageLoader> image_loader, boost::shared_ptr<ImageOutputWriter> output_writer, size_t nFramesInAverage) {
	
	boost::shared_ptr<Eigen::MatrixXd> currentImage;
	size_t xSize = image_loader->getXSize();
	size_t ySize = image_loader->getYSize();
	size_t nFramesInMovie = image_loader->getNImages();
	int nFramesSurroundingFrame = nFramesInAverage / 2;
	int abortStatus;
	
	if (nFramesInAverage > nFramesInMovie) {
		throw std::runtime_error("The number of frames requested in the rolling average is larger than the total number of frames in the movie");
	}
	if (nFramesInAverage % 2 == 0) {
		throw std::runtime_error("Subtracting a rolling average requires an odd number of frames in the average");
	}
	
	boost::shared_ptr<Eigen::MatrixXd> averageImage (new Eigen::MatrixXd((int)xSize, (int)ySize));
	boost::shared_ptr<Eigen::MatrixXd> summedImages (new Eigen::MatrixXd((int)xSize, (int)ySize));
	boost::shared_ptr<Eigen::MatrixXd> activeImage (new Eigen::MatrixXd((int)xSize, (int)ySize));
	
	std::deque<boost::shared_ptr<Eigen::MatrixXd> > frameBuffer;
	
	summedImages->setConstant(0.0);
	
	// start by loading a full set of images
	for (size_t i = 0; i < nFramesInAverage; ++i) {
		currentImage = image_loader->readImage(i);
		frameBuffer.push_front(currentImage);
		(*summedImages) += (*currentImage);
	}
	
	// loop over all frames in the movie
	for (int n = 0; n < nFramesInMovie; ++n) {
		
		// test for user abort and provide a progress update every 10 frames
		if (n % 10 == 0) {
			abortStatus = this->progressReporter->UpdateCalculationProgress((double)n, (double)nFramesInMovie);
			if (abortStatus != 0) {
				// user abort
				this->progressReporter->CalculationAborted();
				throw USER_ABORTED("aborted by user request");
			}
		}
		
		// the queue should be as long as nFramesInAverage at all times
		assert(frameBuffer.size() == nFramesInAverage);
		
		if (n - nFramesSurroundingFrame <= 0) {
			// we're too close to the beginning of the movie to get a full rolling average
			// subtract the average made from the first set of frames instead
			*activeImage = *frameBuffer.at(n);
			
			// do not take the current frame into account when calculating the average
			*summedImages -= *activeImage;
			(*averageImage) = *summedImages / (double)(nFramesInAverage - 1);
			*summedImages += *activeImage;
			
			*activeImage -= *averageImage;
			output_writer->write_image(activeImage);
			
		} else if (n + nFramesSurroundingFrame >= nFramesInMovie) {
			// we're too close to the end of the movie to get a full rolling average
			// subtract the average made from the last set of frames instead
			*activeImage = *frameBuffer.at(nFramesInAverage - (nFramesInMovie - n));
			
			// subtract the contribution of the current frame
			*summedImages -= *activeImage;
			(*averageImage) = *summedImages / (double)(nFramesInAverage - 1);
			*summedImages += *activeImage;
			
			*activeImage -= *averageImage;
			output_writer->write_image(activeImage);
			
		} else {
			// if we're here then we're in the middle of the movie
			// first remove the contribution of the frame that will go out of scope
			// from the sum of frames
			*summedImages -= *(frameBuffer.back());
			frameBuffer.pop_back();
			
			// now add a new frame to the buffer
			frameBuffer.push_front(image_loader->readImage(n + nFramesSurroundingFrame));
			// add the contribution of the new frame
			*summedImages += *(frameBuffer.front());
			
			*activeImage = *frameBuffer.at(nFramesInAverage / 2);
			*summedImages -= *activeImage;
			(*averageImage) = *summedImages / (double)(nFramesInAverage - 1);
			*summedImages += *activeImage;
			
			*activeImage -= *averageImage;
			
			output_writer->write_image(activeImage);
		}
	}
}

void CCDImagesProcessorDifferenceImage::convert_images(boost::shared_ptr<ImageLoader> image_loader, boost::shared_ptr<ImageOutputWriter> output_writer) {
	boost::shared_ptr<Eigen::MatrixXd> current_image;
	boost::shared_ptr<Eigen::MatrixXd> next_image;
	size_t total_number_of_images = image_loader->getNImages();
	
	int abortStatus;
	
	if (total_number_of_images <= 1) {
		throw std::runtime_error("Impossible to calculate a difference image on a sequence that has less than two frames");
	}
	
	// we start by loading the first image in next_image
	// this is required so the loop that follows can start properly
	
	next_image = image_loader->readImage(0);
	
	this->progressReporter->CalculationStarted();
	
	for (size_t n = 0; n < (total_number_of_images - 1); n++) {
		// test for user abort and provide a progress update every 10 frames
		if (n % 10 == 0) {
			abortStatus = this->progressReporter->UpdateCalculationProgress((double)(n + 1), (double)total_number_of_images);
			if (abortStatus != 0) {
				// user abort
				this->progressReporter->CalculationAborted();
				throw USER_ABORTED("aborted by user request");
			}
		}
		
		// the previous image for this run of the loop is the image that was previously in current_image
		// so we have to shift it down
		current_image = next_image;
		
		next_image = image_loader->readImage(n + 1);
		
		// now do the actual subtraction
		*current_image = *current_image - *next_image;
		
		// current_image now contains the subtracted image, we should write it to disk
		output_writer->write_image(current_image);
	}
	this->progressReporter->CalculationDone();
}

void CCDImagesProcessorConvertToSimpleFileFormat::convert_images(boost::shared_ptr<ImageLoader> image_loader, boost::shared_ptr<ImageOutputWriter> output_writer) {
	boost::shared_ptr<Eigen::MatrixXd> current_image;
	size_t total_number_of_images = image_loader->getNImages();
	int abortStatus;
	
	this->progressReporter->CalculationStarted();
	
	for (size_t n = 0; n < total_number_of_images; ++n) {
		// test for user abort and provide a progress update every 10 frames
		if (n % 10 == 0) {
			abortStatus = this->progressReporter->UpdateCalculationProgress((double)n, (double)total_number_of_images);
			if (abortStatus != 0) {
				// user abort
				this->progressReporter->CalculationAborted();
				throw USER_ABORTED("aborted by user request");
			}
		}
		
		current_image = image_loader->readImage(n);
		output_writer->write_image(current_image);
	}
	this->progressReporter->CalculationDone();
}

void CCDImagesProcessorCrop::convert_images(boost::shared_ptr<ImageLoader> image_loader, boost::shared_ptr<ImageOutputWriter> output_writer) {
	size_t total_number_of_images = image_loader->getNImages();
	int abortStatus;
	
	if ((startX >= endX) || (startY >= endY)) {
		throw std::runtime_error("Bad limit values specified for cropping");
	}
	
	if ((endX >= image_loader->getXSize()) || (endY >= image_loader->getYSize())) {
		throw std::runtime_error("Bad limit values specified for cropping");
	}
	
	this->croppedXSize = endX - startX + 1;
	this->croppedYSize = endY - startY + 1;
	
	boost::shared_ptr<Eigen::MatrixXd > croppedImage;
	boost::shared_ptr<Eigen::MatrixXd > loadedImage;
	
	this->progressReporter->CalculationStarted();
	
	for (size_t n = 0; n < total_number_of_images; ++n) {
		// test for user abort and provide a progress update every 10 frames
		if (n % 10 == 0) {
			abortStatus = this->progressReporter->UpdateCalculationProgress((double)n, (double)total_number_of_images);
			if (abortStatus != 0) {
				// user abort
				this->progressReporter->CalculationAborted();
				throw USER_ABORTED("aborted by user request");
			}
		}
		
		loadedImage = image_loader->readImage(n);
		croppedImage = boost::shared_ptr<Eigen::MatrixXd> (new Eigen::MatrixXd((int)croppedXSize, (int)croppedYSize));
		
		for (size_t y = startY; y <= endY; ++y) {
			for (size_t x = startX; x <= endX; ++x) {
				(*croppedImage)(x - startX, y - startY) = (*loadedImage)(x, y);
			}
		}
		
		output_writer->write_image(croppedImage);
	}
	this->progressReporter->CalculationDone();
}

void CCDImagesProcessorConvertToPhotons::convert_images(boost::shared_ptr<ImageLoader> image_loader, boost::shared_ptr<ImageOutputWriter> output_writer) {
	assert(this->multiplicationFactor > 0.0);
	
	size_t total_number_of_images = image_loader->getNImages();
	size_t x_size = image_loader->getXSize();
	size_t y_size = image_loader->getYSize();
	
	boost::shared_ptr<Eigen::MatrixXd > loadedImage;
	boost::shared_ptr<Eigen::MatrixXd > convertedImage (new Eigen::MatrixXd((int)x_size, (int)y_size));;
	
	int abortStatus;
	this->progressReporter->CalculationStarted();
	
	for (size_t n = 0; n < total_number_of_images; ++n) {
		loadedImage = image_loader->readImage(n);
		for (size_t j = 0; j < y_size; ++j) {
			for (size_t i = 0; i < x_size; ++i) {
				(*convertedImage)(i, j) = ((*loadedImage)(i, j) - this->offset >= 0) ? (((*loadedImage)(i, j) - this->offset) / this->multiplicationFactor) : 0.0;
			}
		}
		output_writer->write_image(convertedImage);
		
		abortStatus = this->progressReporter->UpdateCalculationProgress((double)n / (double)(total_number_of_images) * 100.0, 100.0);
		if (abortStatus != 0) {
			progressReporter->CalculationAborted();
			throw USER_ABORTED("");
		}
	}
	this->progressReporter->CalculationDone();
}

