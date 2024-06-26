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

#include "PALMAnalysis.h"

#include "Errors.h"
#include "FileIO.h"
#include "Segmentation.h"
#include "ParticleFinding.h"
#include "Storage.h"


std::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > do_processing_and_thresholding(ImagePtr image, std::shared_ptr<ThresholdImage_Preprocessor>preprocessor, 
																									  std::shared_ptr<ThresholdImage> thresholder, std::shared_ptr<ThresholdImage_Postprocessor> postprocessor) {
	// this function takes care of the thresholding and the associated pre- and postprocessing
	
	ImagePtr preprocessed_image;
	std::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > thresholded_image;
	std::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > postprocessed_image;
	
	preprocessed_image = preprocessor->do_preprocessing(image);
	thresholded_image = thresholder->do_thresholding(preprocessed_image);
	postprocessed_image = postprocessor->do_postprocessing(thresholded_image, image);
	
	return postprocessed_image;
}

PALMAnalysisController::PALMAnalysisController (std::shared_ptr<ThresholdImage> thresholder_rhs,
												std::shared_ptr<ThresholdImage_Preprocessor> thresholdImagePreprocessor_rhs,
												std::shared_ptr<ThresholdImage_Postprocessor> thresholdImagePostprocessor_rhs,
												std::shared_ptr<ParticleFinder> particleFinder_rhs, 
												std::vector<std::shared_ptr<ParticleVerifier> > particleVerifiers_rhs,
												std::shared_ptr<FitPositions> fitPositions_rhs,
												std::shared_ptr<ProgressReporter> progressReporter_rhs,
												size_t firstFrame, size_t lastFrame) {
	thresholder = thresholder_rhs;
	thresholdImagePreprocessor = thresholdImagePreprocessor_rhs;
	thresholdImagePostprocessor = thresholdImagePostprocessor_rhs;
	particleFinder = particleFinder_rhs;
	particleVerifiers = particleVerifiers_rhs;
	fitPositions = fitPositions_rhs;
	progressReporter = progressReporter_rhs;
	
	this->firstFrameToAnalyze = firstFrame;
	this->lastFrameToAnalyze = lastFrame;
	if (lastFrame < firstFrame && (firstFrame != (size_t)-1))
		throw std::runtime_error("invalid range of frames to analyze specified");
	
	this->errorMessage.assign("");
	this->shouldAbort = false;
}

std::shared_ptr<LocalizedPositionsContainer> PALMAnalysisController::DoPALMAnalysis(std::shared_ptr<ImageLoader> imageLoader_rhs) {
	this->imageLoader = imageLoader_rhs;
	this->nImages = imageLoader->getNImages();
	
	size_t numberOfThreads;
	size_t firstFrame, lastFrame;
	size_t nFramesToBeAnalyzed;
	std::vector<std::future<void>> threads;
	int spinProcessStatus, progressStatus;
	size_t nFramesAnalyzed;
	
	this->localizedPositions = std::shared_ptr<LocalizedPositionsContainer>();
	
	if (this->nImages == 0) {	// if there are no images to load, do not do any processing
		return this->localizedPositions;
	}
	
	numberOfThreads = std::thread::hardware_concurrency();
	if (numberOfThreads == 0) {
		numberOfThreads = 1;
	}
	if (numberOfThreads > this->nImages) {
		numberOfThreads = nImages;
	}
	
	firstFrame = this->firstFrameToAnalyze;
	lastFrame = this->lastFrameToAnalyze;
	if (firstFrame == (size_t)-1) {
		firstFrame = 0;
	} else if (firstFrame >= this->nImages) {
		throw std::runtime_error("Invalid first frame to analyze specified");
	}
	if (lastFrame == (size_t)-1) {
		lastFrame = this->nImages - 1;
	} else if ((lastFrame >= this->nImages) || (lastFrame < firstFrame)) {
		throw std::runtime_error("Invalid last frame to analyze specified");
	}
	
	nFramesToBeAnalyzed = lastFrame - firstFrame + 1;
	
	this->nFramesRemainingToBeProcessed = lastFrame - firstFrame + 1;
	
	// spool the image loader so that the correct first frame will be returned
	this->imageLoader->spoolTo(firstFrame);
	
	progressReporter->CalculationStarted();
	
	// start the thread pool
	this->shouldAbort = false;
	for (size_t j = 0; j < numberOfThreads; ++j) {
		threads.push_back(std::async(std::launch::async, [=]() {ThreadPoolWorker(this); }));
	}
	
	// test if the threads have finished
	for (;;) {
		std::future_status firstThreadStatus = threads.at(0).wait_for(std::chrono::milliseconds(250));
		if (firstThreadStatus != std::future_status::ready) {	// the thread is not done yet, we're just waiting
			// while we wait we check for various things and give the interface the chance to update
			
			// did one of the threads run into an error?
			this->errorReportingMutex.lock();
			if (this->errorMessage.length() != 0) {
				errorReportingMutex.unlock();
				// an error occured, time to abort this analysis
				this->shouldAbort = true;
				// wait until the threads have completed
				for (size_t j = 0; j < numberOfThreads; ++j) {
					threads.at(j).get();
				}
				
				throw ERROR_RUNNING_THREADED_ANALYSIS(this->errorMessage);
			}
			this->errorReportingMutex.unlock();
			
			// allow the reporter to update with new progress
			{
				std::lock_guard<std::mutex> locker(this->acquireFrameForProcessingMutex);
				nFramesAnalyzed = nFramesToBeAnalyzed - this->nFramesRemainingToBeProcessed;
			}
			progressStatus = progressReporter->UpdateCalculationProgress((double)nFramesAnalyzed, (double)nFramesToBeAnalyzed);
			
#ifdef WITH_IGOR
			// does the user want to abort?
			spinProcessStatus = SpinProcess();
#else
			spinProcessStatus = 0;
#endif	// WITH_IGOR
			if ((spinProcessStatus != 0) || (progressStatus != 0)) {
				this->shouldAbort = true;
				// wait until the threads have completed
				for (size_t j = 0; j < numberOfThreads; ++j) {
					threads.at(j).get();
				}
				progressReporter->CalculationAborted();
				this->localizedPositions->sortPositionsByFrameNumber();
				return this->localizedPositions;
			}
			continue;
		} else {
			// the first thread has finished
			break;
		}
	}
	
	for (size_t j = 0; j < numberOfThreads; ++j) {
		threads.at(j).get();
	}
	
	// the processing itself is now done, but the results will not have been returned in the correct order
	this->localizedPositions->sortPositionsByFrameNumber();
	
	progressReporter->CalculationDone();
	
	// return the results
	return this->localizedPositions;
}

void ThreadPoolWorker(PALMAnalysisController* controller) {
	int currentImageToProcess;
	ImagePtr currentImage;
	std::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > thresholdedImage;
	std::shared_ptr<std::list<Particle> > locatedParticles;
	std::shared_ptr<LocalizedPositionsContainer> localizedPositions;
	
	try {
		for (;;) {	// loop continuously looking for more images until there are none left
			// if the main thread wants us to interrupt, then give it an opportunity to do so
			if (controller->shouldAbort) {
				return;
			}
			
			// start by a acquiring an image to process
			// see if there are still frames to be processed
			{
				std::lock_guard<std::mutex> locker(controller->acquireFrameForProcessingMutex);
				if (controller->nFramesRemainingToBeProcessed == 0) {
					// no more frames to be processed
					return;
				}
				controller->nFramesRemainingToBeProcessed -= 1;
			}
			
			// we need to process the image with index currentImageToProcess
			currentImage = controller->imageLoader->readNextImage(currentImageToProcess);
			
			thresholdedImage = do_processing_and_thresholding(currentImage, controller->thresholdImagePreprocessor, controller->thresholder,
															  controller->thresholdImagePostprocessor);
			
			locatedParticles = controller->particleFinder->findPositions(currentImage, thresholdedImage);
			
			// if the located particles are to be verified before fitting then do so
			for (std::vector<std::shared_ptr<ParticleVerifier> >::iterator it = controller->particleVerifiers.begin(); it != controller->particleVerifiers.end(); ++it) {
				(*it)->VerifyParticles(currentImage, locatedParticles);
			}
			
			localizedPositions = controller->fitPositions->fit_positions(currentImage, locatedParticles);
			
			// the localization routines have no idea what frame number the processed frame was
			// so we need to add the appropriate frame number to the fitted positions here
			localizedPositions->setFrameNumbers(currentImageToProcess);
			
			// pass the result to the output queue
			{
				std::lock_guard<std::mutex> locker(controller->addLocalizedPositionsMutex);
				// if this is the first time that positions are being returned then controller will contain a NULL pointer
				// set it to the positions we are now returning
				// this way ThreadPoolWorker does not need to know what the type of positions is
				if (controller->localizedPositions.get() == NULL) {
					controller->localizedPositions = localizedPositions;
				} else {
					controller->localizedPositions->addPositions(localizedPositions);
				}
			}
		}
	}
	catch (std::runtime_error &e) {
		// an error has appeared somewhere in this thread during the calculation
		// since there seems to be no easy way to communicate the exception to the
		// main thread, set an error message in the analysis controller.
		// The controller will periodically check this message and handle the error
		std::lock_guard<std::mutex> locker(controller->errorReportingMutex);
		controller->errorMessage.assign(e.what());
		
		// no point in continuing this thread
		return;
	}
	catch (std::exception &e) {
		std::lock_guard<std::mutex> locker(controller->errorReportingMutex);
		controller->errorMessage.assign(e.what());
		
		return;
	}
	catch (...) {
		// catch any other exception not handled by the above block
		std::lock_guard<std::mutex> locker(controller->errorReportingMutex);
		controller->errorMessage.assign("Encountered an unspecified exception while doing the PALM analysis");
		
		// no point in continuing this thread
		return;
	}
	
}

std::shared_ptr<LocalizedPositionsContainer> FitPositionsDeflate::fit_positions(const ImagePtr image, std::shared_ptr<std::list<Particle> > positions) {
	// TODO: for now we ignore the starting positions and ending position provided as arguments
	
	std::shared_ptr<LocalizedPositionsContainer> positionsFittedThusFar;
	std::shared_ptr<LocalizedPositionsContainer> positionsLocalizedThisFrame;
	ImagePtr subtractedImage;
	std::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > segmentedImage;
	std::shared_ptr<std::list<Particle> > locatedParticles;
	
	// get an initial set of positions to start from
	positionsFittedThusFar = this->positionsFitter->fit_positions(image, positions);
	
	// if no positions were localized then don't bother
	if (positionsFittedThusFar->getNPositions() == 0)
		return positionsFittedThusFar;
	
	// check if the type of FitPositions provided returns the required information to reconstruct
	// a Gaussian emitter
	if ((positionsFittedThusFar->getIntegral(0) == 0) || (positionsFittedThusFar->getXWidth(0) == 0) || (positionsFittedThusFar->getYWidth(0) == 0) || (positionsFittedThusFar->getBackground(0) == 0))
		throw std::runtime_error("The selected fitting algorithm does not provide sufficient information for deflation analysis");
	
	// set up for iteration
	positionsLocalizedThisFrame = positionsFittedThusFar;
	subtractedImage = image;
	
	while (positionsLocalizedThisFrame->getNPositions() != 0) {
		// continue while new positions are still being localized
		subtractedImage = this->subtractLocalizedPositions(subtractedImage, positionsLocalizedThisFrame);
		
		segmentedImage = do_processing_and_thresholding(subtractedImage, this->preprocessor, 
														this->thresholder, this->postprocessor);
		locatedParticles = this->particleFinder->findPositions(subtractedImage, segmentedImage);
		positionsLocalizedThisFrame = this->positionsFitter->fit_positions(subtractedImage, locatedParticles);
		
		// if we found new particles then append them
		if (positionsLocalizedThisFrame->getNPositions() != 0) {
			positionsFittedThusFar->addPositions(positionsLocalizedThisFrame);
		}
	}
	
	return positionsFittedThusFar;
	
}


ImagePtr FitPositionsDeflate::subtractLocalizedPositions(ImagePtr image, std::shared_ptr<LocalizedPositionsContainer> positions) {
	double fittedXPos, fittedYPos, fittedIntegral, fittedXWidth, fittedYWidth;
	double centerX, centerY, calculatedAmplitude;
	long startX, endX, startY, endY;
	double distanceXSquared, distanceYSquared, currentIntensity;
	
	size_t nPositions = positions->getNPositions();
	int xSize = image->rows();
	int ySize = image->cols();
	ImagePtr outputImage(new Image((int)xSize, (int)ySize));
	outputImage = image;
	
	for (size_t n = 0; n < nPositions; ++n) {
		fittedIntegral = positions->getIntegral(n);
		fittedXWidth = positions->getXWidth(n);
		fittedYWidth = positions->getYWidth(n);
		fittedXPos = positions->getXPosition(n);
		fittedYPos = positions->getYPosition(n);
		
		calculatedAmplitude = fittedIntegral / (2 * PI * fittedXWidth * fittedYWidth);
		
		centerX = fittedXPos;
		centerY = fittedYPos;
		
		startX = floor((double)centerX - 4.0 * fittedXWidth);	// only run the calculation over a subset of the image surrounding the position
		startY = floor((double)centerY - 4.0 * fittedYWidth);
		endX = ceil((double)centerX + 4.0 * fittedXWidth);
		endY = ceil((double)centerY + 4.0 * fittedYWidth);
		
		if (startX < 0)
			startX = 0;
		if (endX >= xSize)
			endX = xSize - 1;
		if (startY < 0)
			startY = 0;
		if (endY >= ySize)
			endY = ySize - 1;
		
		for (int j = startY; j < endY; ++j) {
			for (int i = startX; i <= endX; ++i) {
				distanceXSquared = ((double)i - (double)centerX) * ((double)i - (double)centerX);
				distanceYSquared = ((double)j - (double)centerY) * ((double)j - (double)centerY);
				currentIntensity = calculatedAmplitude * exp(- (distanceXSquared + distanceYSquared) / (2 * fittedXWidth * fittedYWidth));
				
				(*outputImage)(i, j) = (*outputImage)(i, j) - currentIntensity;
			}
		}
	}
	return outputImage;
}
