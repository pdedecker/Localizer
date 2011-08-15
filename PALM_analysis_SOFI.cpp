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
					std::vector<boost::shared_ptr<SOFIFrameVerifier> > frameVerifiers, boost::shared_ptr<ProgressReporter> progressReporter,
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
	size_t nFramesProcessedTotal = 0;
	ImagePtr outputImage;
	ImagePtr sofiImage;
	size_t nImagesInBlock;
	int status;
	int isValidFrame;
	
	double weightOfThisBlock, sumOfBlockWeights;
	
	progressReporter->CalculationStarted();
	
	for (size_t currentGroup = 0; currentGroup < nGroups; currentGroup += 1) {
		if (nFramesToGroup > 0) {
			nFramesInThisGroup = std::min(nFramesToGroup, (int)(nImagesToProcess - (currentGroup * nFramesToGroup)));
		} else {
			nFramesInThisGroup = nImagesToProcess;
		}
		
		nFramesProcessedInThisRun = 0;
		nImagesInBlock = 0;
		sumOfBlockWeights = 0.0;
		while (nFramesProcessedInThisRun < nFramesInThisGroup) {
			nImagesInBlock = std::min(nFramesInThisGroup - nFramesProcessedInThisRun, blockSize);
			if (nImagesInBlock < lagTime + 1)
				break;
			
			for (size_t i = 0; i < nImagesInBlock; ++i) {
				currentImage = imageLoader->readNextImage();
				
				// verify if this is an acceptable image
				isValidFrame = 1;
				for (size_t j = 0; j < frameVerifiers.size(); j+=1) {
					if (frameVerifiers[j]->isValidFrame(currentImage) != 1) {
						isValidFrame = 0;
						break;
					}
				}
				if (isValidFrame == 1) {
					sofiCalculator->addNewImage(currentImage);
				}
				nFramesProcessedInThisRun += 1;
				nFramesProcessedTotal += 1;
			}
			
			status = progressReporter->UpdateCalculationProgress(nFramesProcessedTotal, nImagesToProcess);
			if (status != 0) {
				progressReporter->CalculationAborted();
				throw USER_ABORTED("Abort requested by user");
			}
			
			outputImage = sofiCalculator->getResult();
			weightOfThisBlock = (double)(nImagesInBlock) / (double)(blockSize);
			sumOfBlockWeights += weightOfThisBlock;
				
			if (sofiImage.get() == NULL) {
				sofiImage = outputImage;
			} else {
				*sofiImage += *outputImage * weightOfThisBlock;
			}
		}
		
		if (nImagesInBlock > 0) {
			*sofiImage /= sumOfBlockWeights;
			
			if (crossCorrelate != 0)
				sofiImage = sofiCorrector.doImageCorrection(sofiImage);
			outputWriter->write_image(sofiImage);
		}
	}
	progressReporter->CalculationDone();
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
	if (this->outputImage.get() == NULL) {
		this->outputImage = ImagePtr(new Image((int)nRows, (int)nCols));
		this->outputImage->setConstant(0.0);
		this->averageImage = ImagePtr(new Image((int)nRows, (int)nCols));
		this->averageImage->setConstant(0.0);
	}
	
	ImagePtr previousImage = this->imageQueue.front();
	ImagePtr currentImage = this->imageQueue.back();
	
	*this->outputImage += ((*previousImage).array() * (*currentImage).array()).matrix();
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
	*this->outputImage -= (*this->averageImage).array().square().matrix();
	
	ImagePtr imageToBeReturned(new Image(*this->outputImage));
	
	// now reset everything for the next calculation
	// before returning
	this->nEvaluations = 0;
	this->averageImage->setConstant(0.0);
	this->outputImage->setConstant(0.0);
	
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
	
	size_t nRowsOutputImage = 2 * nRows - 4;
	size_t nColsOutputImage = 2 * nCols - 4;
	
	if (this->nEvaluations == 0) {
		this->outputImageCrossCorrelation = ImagePtr(new Image((int)(nRowsOutputImage), (int)(nColsOutputImage)));
		this->outputImageHorizontalAutoCorrelation = ImagePtr(new Image((int)(nRowsOutputImage), (int)(nColsOutputImage)));
		this->outputImageVerticalAutoCorrelation = ImagePtr(new Image((int)(nRowsOutputImage), (int)(nColsOutputImage)));
		
		this->outputImageCrossCorrelation->setConstant(0.0);
		this->outputImageHorizontalAutoCorrelation->setConstant(0.0);
		this->outputImageVerticalAutoCorrelation->setConstant(0.0);
		
		this->averageImage = ImagePtr(new Image((int)nRows, (int)nCols));
		this->averageImage->setConstant(0.0);
	}
	
	ImagePtr previousImage = this->imageQueue.front();
	ImagePtr currentImage = this->imageQueue.back();
	
	*this->averageImage += *previousImage;
	
	for (size_t j = 1; j < nCols - 1; ++j) {
		for (size_t i = 1; i < nRows - 1; ++i) {
			// first handle the autocorrelation
			// these autocorrelation pixels are now calculated using crosscorrelation
			// since these require the correlation of two pixels with an intermediate pixel
			// in between, they might be more sensitive to noise. Therefore do two different
			// correlations, and at the end average between them
			(*outputImageHorizontalAutoCorrelation)(2*i - 2, 2*j - 2) += (*previousImage)(i - 1, j) * (*currentImage)(i + 1, j);
			(*outputImageVerticalAutoCorrelation)(2*i - 2, 2*j - 2) += (*previousImage)(i, j - 1) * (*currentImage)(i, j + 1);
			
			// handle the horizontal and vertical pixel
			(*outputImageCrossCorrelation)(2*i - 2 + 1, 2*j - 2) += (*previousImage)(i, j) * (*currentImage)(i + 1, j);
			(*outputImageCrossCorrelation)(2*i - 2, 2*j - 2 + 1) += (*previousImage)(i, j) * (*currentImage)(i, j + 1);
			
			// handle the diagonal pixel
			(*outputImageCrossCorrelation)(2*i - 2 + 1, 2*j - 2 + 1) += (*previousImage)(i, j) * (*currentImage)(i + 1, j + 1);
		}
	}
	
	// remove the last image from the queue
	this->imageQueue.pop();
	this->nEvaluations += 1;
}

ImagePtr SOFICalculator_Order2_cross::getResult() {
	if (this->nEvaluations == 0)
		throw std::runtime_error("Requested a SOFI image even though there are no evaluations");
	
	*this->averageImage /= (double)(this->nEvaluations);
	*this->outputImageHorizontalAutoCorrelation /= (double)(this->nEvaluations);
	*this->outputImageVerticalAutoCorrelation /= (double)(this->nEvaluations);
	*this->outputImageCrossCorrelation /= (double)(this->nEvaluations);
	
	size_t nRowsOutputImage = this->outputImageHorizontalAutoCorrelation->rows();
	size_t nColsOutputImage = this->outputImageHorizontalAutoCorrelation->cols();
	
	// normalize the correlated values to get fluctuations
	for (size_t j = 0; j < nColsOutputImage; ++j) {
		for (size_t i = 0; i < nRowsOutputImage; ++i) {
			
			if ((i % 2 == 0) && (j % 2 == 0)) {
				// this is an autocorrelation pixel
				(*outputImageHorizontalAutoCorrelation)(i, j) -= (*averageImage)(i / 2 - 1 + 1, j / 2 + 1) * (*averageImage)(i / 2 + 1 + 1, j / 2 + 1);
				(*outputImageVerticalAutoCorrelation)(i, j) -= (*averageImage)(i / 2 + 1, j / 2 - 1 + 1) * (*averageImage)(i / 2 + 1, j / 2 + 1 + 1);
				continue;
			}
			
			if ((i % 2 == 1) && (j % 2 == 1)) {
				// this is a diagonal crosscorrelation pixel
				(*outputImageCrossCorrelation)(i, j) -= (*averageImage)(i / 2 + 1, j / 2 + 1) * (*averageImage)(i / 2 + 1 + 1, j / 2 + 1 + 1);
				continue;
			}
			
			if (i % 2 == 1) {
				(*outputImageCrossCorrelation)(i, j) -= (*averageImage)(i / 2 + 1, j / 2 + 1) * (*averageImage)(i / 2 + 1 + 1, j / 2 + 1);
				continue;
			}
			
			if (j % 2 == 1) {
				(*outputImageCrossCorrelation)(i, j) -= (*averageImage)(i / 2 + 1, j / 2 + 1) * (*averageImage)(i / 2 + 1, j / 2 + 1 + 1);
				continue;
			}
		}
	}
	
	ImagePtr imageToReturn (new Image(*outputImageHorizontalAutoCorrelation));
	*imageToReturn = ((*outputImageHorizontalAutoCorrelation) + (*outputImageVerticalAutoCorrelation) / 2.0) + (*outputImageCrossCorrelation);
	
	// now reset everything for the next calculation
	// before returning
	this->nEvaluations = 0;
	this->averageImage->setConstant(0.0);
	this->outputImageHorizontalAutoCorrelation->setConstant(0.0);
	this->outputImageVerticalAutoCorrelation->setConstant(0.0);
	this->outputImageCrossCorrelation->setConstant(0.0);
	
	while (this->imageQueue.size() > 0)
		this->imageQueue.pop();
	
	return imageToReturn;
}

ImagePtr SOFICorrector_Order2::doImageCorrection(ImagePtr imageToCorrect) {
	size_t nRows = imageToCorrect->rows();
	size_t nCols = imageToCorrect->cols();
	
	double avgOfAuto = 0.0, avgOfHoriz = 0.0, avgOfDiag = 0.0;
	double nAuto = 0.0, nHoriz = 0.0, nDiag = 0.0;
	
	ImagePtr correctedImage(new Image(*imageToCorrect));
	
	for (size_t j = 0; j < nCols; j+=1) {
		for (size_t i = 0; i < nRows; i+=1) {
			if ((i % 2 == 0) && (j % 2 == 0)) {
				avgOfAuto += (*imageToCorrect)(i, j);
				nAuto += 1.0;
				continue;
			}
			if ((i % 2 == 1) && (j % 2 == 1)) {
				avgOfDiag += (*imageToCorrect)(i, j);
				nDiag += 1.0;
				continue;
			}
			
			avgOfHoriz += (*imageToCorrect)(i, j);
			nHoriz += 1.0;
		}
	}
	
	avgOfAuto /= nAuto;
	avgOfDiag /= nDiag;
	avgOfHoriz /= nHoriz;
	
	double horizFactor = avgOfHoriz / avgOfAuto;
	double diagFactor = avgOfDiag / avgOfAuto;
	
	for (size_t j = 0; j < nCols; j+=1) {
		for (size_t i = 0; i < nRows; i+=1) {
			if ((i % 2 == 0) && (j % 2 == 0)) {
				continue;
			}
			if ((i % 2 == 1) && (j % 2 == 1)) {
				(*correctedImage)(i, j) /= diagFactor;
				continue;
			}
			
			(*correctedImage)(i, j) /= horizFactor;
		}
	}
	
	return correctedImage;
}

/*ImagePtr SOFICorrector_Order2::doImageCorrection(ImagePtr imageToCorrect) {
	// estimate the PSF standard deviation
	double optimalPSFStdDev = determinePSFStdDev(imageToCorrect);
	
	ImagePtr correctedImage = performPSFCorrection(imageToCorrect.get(), optimalPSFStdDev);
	return correctedImage;
}*/

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
				continue;
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
	
	// calculate the relative variance of the mean
	double avg, avgOfSquares, variance;
	
	avg = (*image).sum() / (double)(nRows * nCols);
	avgOfSquares = (*image).array().square().sum() / (double)(nRows * nCols);
	variance = avgOfSquares - avg * avg;
	
	return variance / (avg * avg);
}

SOFIFrameVerifier_NoSaturation::SOFIFrameVerifier_NoSaturation(int storageType_rhs) {
	switch (storageType_rhs) {
            // some systems (e.g. the Nikon TIRF) never report absolutely saturated pixel
            // values, but instead report a different value close to it
            // so add in some arbitrary margin for the saturation detection
		case STORAGE_TYPE_INT8:
			this->saturationValue = 127.0 - 1.0;
			break;
		case STORAGE_TYPE_UINT8:
			this->saturationValue = 255.0 - 2.0;
			break;
		case STORAGE_TYPE_INT16:
			this->saturationValue = 32767.0 - 200.0;
			break;
		case STORAGE_TYPE_UINT16:
			this->saturationValue = 65535.0 - 200.0;
			break;
		case STORAGE_TYPE_INT32:
			this->saturationValue = 2147483647.0 - 500.0;
			break;
		case STORAGE_TYPE_UINT32:
			this->saturationValue = 4294967295.0 - 500.0;
			break;
		case STORAGE_TYPE_INT64:
			this->saturationValue = 9223372036854775807.0 - 1000.0;
			break;
		case STORAGE_TYPE_UINT64:
			this->saturationValue = 18446744073709551615.0 - 1000.0;
			break;
		case STORAGE_TYPE_INT4:
		case STORAGE_TYPE_UINT4:
		case STORAGE_TYPE_FP32:
		case STORAGE_TYPE_FP64:
			throw std::runtime_error("invalid storage type for saturation verification");
			break;
		default:
			throw std::runtime_error("unknown storage type for saturation verification");
			break;
	}
}

int SOFIFrameVerifier_NoSaturation::isValidFrame(ImagePtr frame) {
	size_t nRows = frame->rows();
	size_t nCols = frame->cols();
	
	size_t nPixels = nRows * nCols;
	double localSaturationValue = this->saturationValue;
	
	double *dataPtr = frame->data();
	for (size_t i = 0; i < nPixels; i+=1) {
		if (dataPtr[i] >= localSaturationValue) {
			return 0;
		}
	}
	
	return 1;
}

int SOFIFrameVerifier_MaxPixelValue::isValidFrame(ImagePtr image) {
	size_t nRows = image->rows();
	size_t nCols = image->cols();
	
	size_t nPixels = nRows * nCols;
	double localMaxPixelValue = this->maxPixelValue;
	
	double *dataPtr = image->data();
	for (size_t i = 0; i < nPixels; i+=1) {
		if (dataPtr[i] > localMaxPixelValue) {
			return 0;
		}
	}
	
	return 1;
}

