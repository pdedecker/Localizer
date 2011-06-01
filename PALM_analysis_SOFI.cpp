/*
 *  PALM_analysis_SOFI.cpp
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 01/06/2011.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "PALM_analysis_SOFI.h"

SOFICalculator_Order2_auto::SOFICalculator_Order2_auto(int lagTime_rhs) {
	this->lagTime = lagTime_rhs;
}

void SOFICalculator_Order2_auto::addNewImage(ImagePtr newImage) {
	// a new image is available
	// this can be the start of the calculation, or we can be
	// somewhere in the middle of the calculation
	
	this->imageQueue.push(newImage);
	
	// if the length of the queue is not equal to the lagTime + 1
	// then it's impossible to produce output
	if (this->imageQueue.size() < this->lagTime + 1)
		return;
	
	// if we're here then the queue is long enough
	size_t nRows = newImage->rows();
	size_t nCols = newImage->cols();
	
	// do the necessary images already exist?
	if (this->nEvaluations == 0) {
		this->outputImage = ImagePtr(new Image((int)nRows, (int)nCols));
		this->outputImage->setConstant(0.0);
		this->averageImage = ImagePtr(new Image((int)nRows, (int)nCols));
		this->averageImage->setConstant(0.0);
	}
	
	ImagePtr previousImage = this->imageQueue.front();
	ImagePtr currentImage = this->imageQueue.back();
	
	*this->outputImage += (*previousImage).cwise() * (*currentImage);
	*this->averageImage += *previousImage;
	
	this->imageQueue.pop();
	this->nEvaluations += 1;
}

ImagePtr SOFICalculator_Order2_auto::getResult() {
	if (this->nEvaluations == 0)
		throw std::runtime_error("Requested a SOFI image even though there are no evaluations");
	
	*this->averageImage /= (double)(this->nEvaluations);
	*this->outputImage /= (double)(this->nEvaluations);
	
	// correct for the fact that we have been multiplying intensities
	// instead of fluctuations
	*this->outputImage -= (*this->averageImage).cwise().square();
	
	ImagePtr imageToBeReturned = this->outputImage;
	
	// now reset everything for the next calculation
	// before returning
	this->nEvaluations = 0;
	this->averageImage.reset();
	this->outputImage.reset();
	
	while (this->imageQueue.size() > 0)
		this->imageQueue.pop();
	
	return imageToBeReturned;
}

void DoSOFIAnalysis(boost::shared_ptr<ImageLoader> imageLoader, boost::shared_ptr<ImageOutputWriter> outputWriter,
									 int lagTime, int order, int crossCorrelate) {
	size_t nImages = imageLoader->getNImages();
	if (nImages <= lagTime)
		throw std::runtime_error("Not enough images for the requested lagtime");
	
	boost::scoped_ptr<SOFICalculator> sofiCalculator(new SOFICalculator_Order2_auto(lagTime));
	
	ImagePtr currentImage;
	imageLoader->rewind();
	size_t dummyIndex;
	for (size_t i = 0; i < nImages; ++i) {
		currentImage = imageLoader->readNextImage(dummyIndex);
		sofiCalculator->addNewImage(currentImage);
	}
	
	ImagePtr outputImage = sofiCalculator->getResult();
	outputWriter->write_image(outputImage);
}
