/*
 *  PALM_analysis_SOFI.cpp
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 01/06/2011.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "PALM_analysis_SOFI.h"

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

SOFICalculator_Order2_auto::SOFICalculator_Order2_auto(int lagTime_rhs) {
	this->lagTime = lagTime_rhs;
	this->nEvaluations = 0;
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
	
	// remove the last image from the queue
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

SOFICalculator_Order2_cross::SOFICalculator_Order2_cross(int lagTime_rhs) {
	this->lagTime = lagTime_rhs;
	this->nEvaluations = 0;
}

void SOFICalculator_Order2_cross::addNewImage(ImagePtr newImage) {
	
	this->imageQueue.push(newImage);
	
	// if the length of the queue is not equal to the lagTime + 1
	// then it's impossible to produce output
	if (this->imageQueue.size() < this->lagTime + 1)
		return;
	
	// if we're here then the queue is long enough
	size_t nRows = newImage->rows();
	size_t nCols = newImage->cols();
	
	size_t nRowsOutputImage = 2 * nRows - 1;
	size_t nColsOutputImage = 2 * nCols - 1;
	
	if (this->nEvaluations == 0) {
		this->outputImage = ImagePtr(new Image((int)(nRowsOutputImage), (int)(nColsOutputImage)));
		this->outputImage->setConstant(0.0);
		this->averageImage = ImagePtr(new Image((int)nRows, (int)nCols));
		this->averageImage->setConstant(0.0);
	}
	
	ImagePtr previousImage = this->imageQueue.front();
	ImagePtr currentImage = this->imageQueue.back();
	
	*this->averageImage += *previousImage;
	
	// the cross correlation image has autocorrelation contributions and crosscorrelation contributions
	// treat all of these separately
	for (size_t j = 0; j < nColsOutputImage - 2; j+=2) {
		for (size_t i = 0; i < nRowsOutputImage - 2; i+=2) {
			// first handle the autocorrelation pixel
			// in the original image this pixel is at (i / 2, j / 2)
			(*outputImage)(i, j) = (*previousImage)(i / 2, j / 2) * (*currentImage)(i / 2, j / 2);
			
			// next handle the horizontal and vertical pixel
			(*outputImage)(i + 1, j) = (*previousImage)(i / 2, j / 2) * (*currentImage)(i / 2 + 1, j / 2);
			(*outputImage)(i, j + 1) = (*previousImage)(i / 2, j / 2) * (*currentImage)(i / 2, j / 2 + 1);
			
			// and the diagonal pixel
			(*outputImage)(i + 1, j + 1) = (*previousImage)(i / 2, j / 2) * (*currentImage)(i / 2 + 1, j / 2 + 1);
		}
	}
	
	// the previous loops have handled all pixels except those at the edges of the image
	// (the row and column with the highest indices), so handle those separately
	size_t index = nColsOutputImage - 1;
	for (size_t i = 0; i < nRowsOutputImage - 2; i += 2) {
		(*outputImage)(i, index) = (*previousImage)(i / 2, index / 2) * (*currentImage)(i / 2, index / 2);
		(*outputImage)(i + 1, index) = (*previousImage)(i / 2, index / 2) * (*currentImage)(i / 2 + 1, index / 2);
	}
	
	index = nRowsOutputImage - 1;
	for (size_t j = 0; j < nColsOutputImage - 2; j += 2) {
		(*outputImage)(index, j) = (*previousImage)(index / 2, j / 2) * (*currentImage)(index / 2, j / 2);
		
		(*outputImage)(index, j + 1) = (*previousImage)(index / 2, j / 2) * (*currentImage)(index / 2, j / 2 + 1);
	}
	
	// now there's just a single autocorrelation pixel left in the farthest corner of the image
	(*outputImage)(nRowsOutputImage - 1, nColsOutputImage - 1) = (*previousImage)((nRowsOutputImage - 1) / 2, (nColsOutputImage - 1) / 2) 
																	* (*currentImage)((nRowsOutputImage - 1) / 2, (nColsOutputImage - 1) / 2);
	
	// remove the last image from the queue
	this->imageQueue.pop();
	this->nEvaluations += 1;
}

ImagePtr SOFICalculator_Order2_cross::getResult() {
	if (this->nEvaluations == 0)
		throw std::runtime_error("Requested a SOFI image even though there are no evaluations");
	
	*this->averageImage /= (double)(this->nEvaluations);
	*this->outputImage /= (double)(this->nEvaluations);
	
	size_t nRowsOutputImage = this->outputImage->rows();
	size_t nColsOutputImage = this->outputImage->cols();
	// normalize the correlated values to get fluctuations
	for (size_t j = 0; j < nColsOutputImage; ++j) {
		for (size_t i = 0; i < nRowsOutputImage; ++i) {
			if ((i % 2 == 0) && (j & 2 == 0)) {
				// this is an autocorrelation pixel
				(*outputImage)(i, j) /= (*averageImage)(i / 2, j / 2) * (*averageImage)(i / 2, j / 2);
				continue;
			}
			
			if ((i % 2 == 1) && (j % 2 == 1)) {
				// this is a diagonal crosscorrelation pixel
				(*outputImage)(i, j) /= (*averageImage)(i / 2, j / 2) * (*averageImage)(i / 2 + 1, j / 2 + 1);
				continue;
			}
			
			if (i % 2 == 0) {
				(*outputImage)(i, j) /= (*averageImage)(i / 2, j / 2) * (*averageImage)(i / 2 + 1, j / 2);
			} else {	// i is even, j is odd
				(*outputImage)(i, j) /= (*averageImage)(i / 2, j / 2) * (*averageImage)(i / 2, j / 2 + 1);
			}
		}
	}
	
	// now we need to find the size of the psf and correct for that
	double psfStdDev = determinePSFStdDev(outputImage);
	ImagePtr correctedImage = performPSFCorrection(outputImage.get(), psfStdDev);
	return correctedImage;
}

double SOFICalculator_Order2_cross::determinePSFStdDev(ImagePtr image) {
	
	const gsl_min_fminimizer_type * minizerType = gsl_min_fminimizer_brent;
	gsl_min_fminimizer *minimizer = gsl_min_fminimizer_alloc(minizerType);
	if (minimizer == NULL) {
		throw std::bad_alloc();
	}
	
	gsl_function func;
	func.function = SOFICalculator_Order2_cross::functionToMinimize;
	func.params = (void *)(&(*image));
	
	gsl_min_fminimizer_set(minimizer, &func, 1.6, 0.5, 10.0);
	
	int status;
	int maxIterations = 100, iterations = 0;
	for (;;) {
		status = gsl_min_fminimizer_iterate(minimizer);
		if (status != GSL_SUCCESS) {
			gsl_min_fminimizer_free(minimizer);
			throw std::runtime_error("gsl_min_fminimizer_iterate reported an error");
		}
		
		status = gsl_min_test_interval(gsl_min_fminimizer_x_lower(minimizer), gsl_min_fminimizer_x_upper(minimizer), 0.001, 0.0);
		if (status == GSL_SUCCESS)
			break;
		
		iterations += 1;
		if (iterations > maxIterations)
			break;
	}
	
	double minimum = gsl_min_fminimizer_minimum(minimizer);
	gsl_min_fminimizer_free(minimizer);
	
	return minimum;
}

ImagePtr SOFICalculator_Order2_cross::performPSFCorrection(Image* image, double psfStdDev) {
	size_t nRows = image->rows();
	size_t nCols = image->cols();
	
	double horizontalFactor = exp((1 / sqrt(2.0)) / (2 * psfStdDev * psfStdDev));
	double diagonalFactor = exp(1.0 / (2 * psfStdDev * psfStdDev));	// the 1.0 comes from sqrt(2.0) / sqrt(2.0)
	
	ImagePtr correctedImage(new Image(*image));
	
	// only loop over the crosscorrelation pixels
	for (size_t j = 1; j < nCols; j+=2) {
		for (size_t i = 1; i < nRows; i+=2) {
			if ((i % 2 == 1) && (j % 2 == 1)) {
				// this is a diagonal crosscorrelation pixel
				(*correctedImage)(i, j) /= diagonalFactor;
			} else {
				(*correctedImage)(i, j) /= horizontalFactor;
			}
		}
	}
	
	return correctedImage;
}

double SOFICalculator_Order2_cross::functionToMinimize(double psfStdDev, void *params) {
	Image* image = (Image *)params;
	
	size_t nRows = image->rows();
	size_t nCols = image->cols();
	
	ImagePtr correctedImage = SOFICalculator_Order2_cross::performPSFCorrection(image, psfStdDev);
	
	// calculate the mean of the cross and autocorrelation pixels
	double nPixelsAuto = 0.0, nPixelsCross = 0.0;
	double sumOfAuto = 0, sumOfCross = 0;
	for (size_t j = 0; j < nCols; ++j) {
		for (size_t i = 0; i < nRows; ++i) {
			if ((i % 2 == 0) && (j % 2 == 0)) {
				sumOfAuto += (*image)(i, j);
				nPixelsAuto += 1.0;
				continue;
			} else {
				sumOfCross += (*image)(i, j);
				nPixelsCross += 1.0;
			}
		}
	}
	
	sumOfAuto /= nPixelsAuto;
	sumOfCross /= nPixelsCross;
	
	return (sumOfAuto - sumOfCross) * (sumOfAuto - sumOfCross);
}
