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
					size_t nFramesToSkip, int lagTime, int order, int crossCorrelate, int nFramesToGroup, double psfWidth) {
	size_t nImages = imageLoader->getNImages();
	size_t blockSize = 50;
	
	if (nImages <= lagTime)
		throw std::runtime_error("Not enough images for the requested lagtime");
	
	if (nFramesToSkip > nImages - lagTime)
		throw std::runtime_error("Too many frames are being skipped (/SKIP flag)");
	
	if ((nFramesToGroup <= lagTime) && (nFramesToGroup != 0))
		throw std::runtime_error("Cannot group less frames than the lag time");
	
	boost::shared_ptr<SOFICalculator> sofiCalculator;
	if (crossCorrelate == 0) {
		sofiCalculator = boost::shared_ptr<SOFICalculator>(new SOFICalculator_Order2_auto(lagTime));
	} else {
		sofiCalculator = boost::shared_ptr<SOFICalculator>(new SOFICalculator_Order2_cross(lagTime, psfWidth));
	}
	
	SOFICorrector_Order2 sofiCorrector;	// will only be used when crosscorrelating
	
	size_t nImagesToProcess = nImages - nFramesToSkip;
	int nGroups;
	if (nFramesToGroup != 0) {
		nGroups = ceil((double)(nImagesToProcess) / (double)(nFramesToGroup));
	} else {
		nGroups = 1;
	}
	
	ImagePtr currentImage;
	imageLoader->spoolTo(nFramesToSkip);
	
	size_t nFramesInThisGroup;
	size_t nFramesProcessedInThisRun;
	ImagePtr outputImage;
	ImagePtr sofiImage;
	size_t nImagesInBlock;
	
	for (size_t currentGroup = 0; currentGroup < nGroups; currentGroup += 1) {
		if (nFramesToGroup > 0) {
			nFramesInThisGroup = std::min(nFramesToGroup, (int)(nImagesToProcess - (currentGroup * nFramesToGroup)));
		} else {
			nFramesInThisGroup = nImagesToProcess;
		}
		
		nFramesProcessedInThisRun = 0;
		nImagesInBlock = 0;
		while (nFramesProcessedInThisRun < nFramesInThisGroup) {
			nImagesInBlock = std::min(nFramesInThisGroup - nFramesProcessedInThisRun, blockSize);
			if (nImagesInBlock < lagTime + 1)
				break;
			
			for (size_t i = 0; i < nImagesInBlock; ++i) {
				currentImage = imageLoader->readNextImage();
				sofiCalculator->addNewImage(currentImage);
				nFramesProcessedInThisRun += 1;
			}
			
			outputImage = sofiCalculator->getResult();
				
			if (sofiImage.get() == NULL) {
				sofiImage = outputImage;
			} else {
				*sofiImage += *outputImage;
			}
		}
		
		if (nImagesInBlock > 0) {
			*sofiImage /= (double)(nImagesInBlock);
			
			//if (crossCorrelate != 0)
			//	sofiImage = sofiCorrector.doImageCorrection(sofiImage);
			outputWriter->write_image(sofiImage);
		}
	}
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
	
	// if the length of the queue is less than lagTime + 1
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

SOFICalculator_Order2_cross::SOFICalculator_Order2_cross(int lagTime_rhs, double psfWidth_rhs) {
	this->lagTime = lagTime_rhs;
	this->psfWidth = psfWidth_rhs;
	this->nEvaluations = 0;
}

void SOFICalculator_Order2_cross::addNewImage(ImagePtr newImage) {
	
	this->imageQueue.push(newImage);
	
	// if the length of the queue is less than lagTime + 1
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
	
	for (size_t j = 0; j < nCols - 1; ++j) {
		for (size_t i = 0; i < nRows - 1; ++i) {
			// first handle the autocorrelation
			if ((i != 0) && (j != 0)) {
				// these autocorrelation pixels are now calculated using crosscorrelation
				(*outputImage)(2*i, 2*j) += (*previousImage)(i, j - 1) * (*currentImage)(i, j + 1);
			}
			
			// handle the horizontal and vertical pixel
			(*outputImage)(2*i + 1, 2*j) += (*previousImage)(i, j) * (*currentImage)(i + 1, j);
			(*outputImage)(2*i, 2*j + 1) += (*previousImage)(i, j) * (*currentImage)(i, j + 1);
			
			// handle the diagonal pixel
			(*outputImage)(2*i + 1, 2*j + 1) += (*previousImage)(i, j) * (*currentImage)(i + 1, j + 1);
		}
	}
	
	// the previous loops have handled all pixels except those at the edges of the image
	// (the row and column with the highest indices), so handle those separately
	size_t index = nCols - 1;
	for (size_t i = 0; i < nRows - 1; ++i) {
		// autocorrelation
		(*outputImage)(2 * i, 2 * index) += (*previousImage)(i, index) * (*currentImage)(i, index);
		
		// crosscorrelation
		(*outputImage)(2 * i + 1, 2 * index) += (*previousImage)(i, index) * (*currentImage)(i + 1, index);
	}
	
	index = nRows - 1;
	for (size_t j = 0; j < nCols - 1; ++j) {
		// autocorrelation
		(*outputImage)(2 * index, 2 * j) += (*previousImage)(index, j) * (*currentImage)(index, j);
		
		// crosscorrelation
		(*outputImage)(2 * index, 2 * j + 1) += (*previousImage)(index, j) * (*currentImage)(index, j + 1);
	}
	
	// now there's just a single autocorrelation pixel left in the farthest corner of the image
	(*outputImage)(2 * (nRows - 1), 2 * (nCols - 1)) += (*previousImage)(nRows - 1, nCols - 1) * (*currentImage)(nRows - 1, nCols - 1);
	
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
			if ((i % 2 == 0) && (j % 2 == 0)) {
				// this is an autocorrelation pixel
				if ((i == 0) || (i == nRowsOutputImage - 1) || (j == 0) || (j == nColsOutputImage - 1)) {
					// these pixels are on the edge of the image, it's impossible to crosscorrelate them
					continue;
				}
				
				(*outputImage)(i, j) -= (*averageImage)(i / 2, j / 2 - 1) * (*averageImage)(i / 2, j / 2 + 1);
				continue;
			}
			
			if ((i % 2 == 1) && (j % 2 == 1)) {
				// this is a diagonal crosscorrelation pixel
				(*outputImage)(i, j) -= (*averageImage)(i / 2, j / 2) * (*averageImage)(i / 2 + 1, j / 2 + 1);
				continue;
			}
			
			if (i % 2 == 1) {
				(*outputImage)(i, j) -= (*averageImage)(i / 2, j / 2) * (*averageImage)(i / 2 + 1, j / 2);
				continue;
			}
			
			if (j % 2 == 1) {
				(*outputImage)(i, j) -= (*averageImage)(i / 2, j / 2) * (*averageImage)(i / 2, j / 2 + 1);
				continue;
			}
		}
	}
	
	ImagePtr imageToReturn (new Image(*outputImage));
	
	// now reset everything for the next calculation
	// before returning
	this->nEvaluations = 0;
	this->averageImage.reset();
	this->outputImage.reset();
	
	while (this->imageQueue.size() > 0)
		this->imageQueue.pop();
	
	return imageToReturn;
}

/*ImagePtr SOFICalculator_Order2_cross::performCorrection_Averages(ImagePtr image) {
	
	size_t nRows = image->rows();
	size_t nCols = image->cols();
	
	ImagePtr correctedImage (new Image(*image));
	
	size_t nAutoPixels = (nRows / 2 + 1) * (nCols / 2 + 1);
	size_t nHorizontalPixels = (nRows / 2 + 1) * (nCols / 2);
	size_t nVerticalPixels = (nRows / 2) * (nCols / 2 + 1);
	size_t nDiagonalPixels = (nRows / 2) * (nCols / 2);
	
	double sumAuto = 0.0;
	double sumHorizontal = 0.0;
	double sumVertical = 0.0;
	double sumDiagonal = 0.0;
	
	for (size_t j = 0; j < nCols; j+=1) {
		for (size_t i = 0; i < nRows; i+=1) {
			if ((i % 2 == 0) && (j % 2 == 0)) {
				sumAuto += (*image)(i, j);
				continue;
			}
			
			if ((i % 2 == 1) && (j % 2 == 1)) {
				sumDiagonal += (*image)(i, j);
				continue;
			}
			
			if (i % 2 == 1) {
				sumVertical += (*image)(i, j);
			} else {
				sumHorizontal += (*image)(i, j);
			}
		}
	}
	
	sumAuto /= (double)nAutoPixels;
	sumHorizontal /= (double)nHorizontalPixels;
	sumVertical /= (double)nVerticalPixels;
	sumDiagonal /= (double)nDiagonalPixels;
	
	double horizontalFactor = sumAuto / sumHorizontal;
	double verticalFactor = sumAuto / sumVertical;
	double diagonalFactor = sumAuto / sumDiagonal;
	
	for (size_t j = 0; j < nCols; j+=1) {
		for (size_t i = 0; i < nRows; i+=1) {
			if ((i % 2 == 0) && (j % 2 == 0)) {
				continue;
			}
			
			if ((i % 2 == 1) && (j % 2 == 1)) {
				(*correctedImage)(i, j) *= diagonalFactor;
			}
			
			if (i % 2 == 1) {
				(*correctedImage)(i, j) *= verticalFactor;
			} else {
				(*correctedImage)(i, j) *= horizontalFactor;
			}
		}
	}
	
	return correctedImage;
}*/

ImagePtr SOFICorrector_Order2::doImageCorrection(ImagePtr imageToCorrect) {
	// estimate the PSF standard deviation
	double optimalPSFStdDev = determinePSFStdDev(imageToCorrect);
	
	ImagePtr correctedImage = performPSFCorrection(imageToCorrect.get(), optimalPSFStdDev);
	return correctedImage;
}

double SOFICorrector_Order2::determinePSFStdDev(ImagePtr imageToCorrect) {
	const gsl_min_fminimizer_type * minizerType = gsl_min_fminimizer_brent;
	gsl_min_fminimizer *minimizer = gsl_min_fminimizer_alloc(minizerType);
	if (minimizer == NULL) {
		throw std::bad_alloc();
	}
	
	gsl_function func;
	func.function = SOFICorrector_Order2::functionToMinimize;
	func.params = (void *)(&(*imageToCorrect));
	
	int status;
	
	// try to initialize the minimizer
	// somewhat annoyingly, the solver refuses to do anything
	// unless the boundaries and starting value indicate a minimum
	// so if we don't get lucky right away then search some more
	// for a starting guess
	double initialPSFWidth = 2.0;	// TODO: provide this initial value in a more sensible way
	status = gsl_min_fminimizer_set(minimizer, &func, initialPSFWidth, 0.1, 10.0);
	if (status == GSL_EINVAL) {
		double updatedInitialGuess = 0.1 + 0.5;
		do {
			status = gsl_min_fminimizer_set(minimizer, &func, updatedInitialGuess, 0.1, 10);
			updatedInitialGuess += 0.5;
		} while ((status == GSL_EINVAL) && (updatedInitialGuess < 10.0));
	}
	
	if (status == GSL_EINVAL) {
		// we didn't find any acceptable interval
		throw std::runtime_error("Unable to find acceptable starting guess for the PSF correction in 2nd order XC");
	}
	
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

ImagePtr SOFICorrector_Order2::performPSFCorrection(Image *image, double psfStdDev) {
	size_t nRows = image->rows();
	size_t nCols = image->cols();
	
	double autoPixelFactor = exp(- sqrt(2.0) / (2.0 * psfStdDev * psfStdDev));
	double horizontalFactor = exp(- (1.0 / 2.0) / (2.0 * psfStdDev * psfStdDev));
	double diagonalFactor = exp(- 1.0 / (2.0 * psfStdDev * psfStdDev));	// the 1.0 comes from sqrt(2.0) / sqrt(2.0)
	
	ImagePtr correctedImage(new Image(*image));
	
	// only loop over the crosscorrelation pixels
	for (size_t j = 0; j < nCols; ++j) {
		for (size_t i = 0; i < nRows; ++i) {
			if ((i % 2 == 0) && (j % 2 == 0)) {
				// autocorrelation pixel
				(*correctedImage)(i, j) /= autoPixelFactor;
			}
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

double SOFICorrector_Order2::functionToMinimize(double psfStdDev, void *params) {
	Image* image = (Image *)params;
	
	size_t nRows = image->rows();
	size_t nCols = image->cols();
	
	ImagePtr correctedImage = performPSFCorrection(image, psfStdDev);
	
	// calculate the mean of the cross and autocorrelation pixels
	double nPixelsAuto = 0.0, nPixelsCross = 0.0;
	double sumOfAuto = 0, sumOfCross = 0;
	for (size_t j = 0; j < nCols; ++j) {
		for (size_t i = 0; i < nRows; ++i) {
			if ((i % 2 == 0) && (j % 2 == 0)) {
				sumOfAuto += (*correctedImage)(i, j);
				nPixelsAuto += 1.0;
				continue;
			} else {
				sumOfCross += (*correctedImage)(i, j);
				nPixelsCross += 1.0;
			}
		}
	}
	
	sumOfAuto /= nPixelsAuto;
	sumOfCross /= nPixelsCross;
	
	return (sumOfAuto - sumOfCross) * (sumOfAuto - sumOfCross);
}
