#include "PALM_analysis.h"


boost::shared_ptr<ublas::matrix <unsigned char> > do_processing_and_thresholding(boost::shared_ptr<ublas::matrix<double> > image, boost::shared_ptr<ThresholdImage_Preprocessor>preprocessor, 
																		 boost::shared_ptr<ThresholdImage> thresholder, boost::shared_ptr<ThresholdImage_Postprocessor> postprocessor) {
	// this function takes care of the thresholding and the associated pre- and postprocessing
	
	boost::shared_ptr<ublas::matrix<double> > preprocessed_image;
	boost::shared_ptr<ublas::matrix <unsigned char> > thresholded_image;
	boost::shared_ptr<ublas::matrix <unsigned char> > postprocessed_image;
	
	if (preprocessor.get() != NULL) {
		preprocessed_image = preprocessor->do_preprocessing(image);
		thresholded_image = thresholder->do_thresholding(preprocessed_image);
		
	} else {
		thresholded_image = thresholder->do_thresholding(image);
	}
	
	if (postprocessor.get() != NULL) {
		postprocessed_image = postprocessor->do_postprocessing(thresholded_image, image);
		thresholded_image = postprocessed_image;
	}
	
	return thresholded_image;
}


PALMAnalysisController::PALMAnalysisController (boost::shared_ptr<ThresholdImage> thresholder_rhs,
												boost::shared_ptr<ThresholdImage_Preprocessor> thresholdImagePreprocessor_rhs,
												boost::shared_ptr<ThresholdImage_Postprocessor> thresholdImagePostprocessor_rhs,
												boost::shared_ptr<ParticleFinder> particleFinder_rhs, 
												std::vector<boost::shared_ptr<ParticleVerifier> > particleVerifiers_rhs,
												boost::shared_ptr<FitPositions> fitPositions_rhs,
												boost::shared_ptr<PALMAnalysisProgressReporter> progressReporter_rhs) {
	thresholder = thresholder_rhs;
	thresholdImagePreprocessor = thresholdImagePreprocessor_rhs;
	thresholdImagePostprocessor = thresholdImagePostprocessor_rhs;
	particleFinder = particleFinder_rhs;
	particleVerifiers = particleVerifiers_rhs;
	fitPositions = fitPositions_rhs;
	progressReporter = progressReporter_rhs;
	
	this->errorMessage.assign("");
}

boost::shared_ptr<LocalizedPositionsContainer> PALMAnalysisController::DoPALMAnalysis(boost::shared_ptr<ImageLoader> imageLoader_rhs) {
	this->imageLoader = imageLoader_rhs;
	this->nImages = imageLoader->get_total_number_of_images();
	
	size_t numberOfProcessors = boost::thread::hardware_concurrency();
	size_t numberOfThreads;
	std::vector<boost::shared_ptr<boost::thread> > threads;
	boost::shared_ptr<boost::thread> singleThreadPtr;
	int firstThreadHasFinished, status;
	double percentDone;
	
	this->localizedPositions = boost::shared_ptr<LocalizedPositionsContainer>();
	
	if (this->nImages == 0) {	// if there are no images to load, do not do any processing
		return this->localizedPositions;
	}
	
	numberOfThreads = numberOfProcessors + 1;	// take one extra thread since every thread will be blocked on I/O sooner or later
	if (numberOfThreads == 0) {
		numberOfThreads = 1;
	}
	if (numberOfThreads > this->nImages) {
		numberOfThreads = nImages;
	}
	
	// TODO: REMOVE
	numberOfThreads = 1;
	
	// fill the queue holding the frames to be processed with the frames in the sequence
	for (size_t i = 0; i < this->nImages; ++i) {
		this->framesToBeProcessed.push(i);
	}
	
	progressReporter->CalculationStarted();
	
	// start the thread pool
	threads.clear();
	for (size_t j = 0; j < numberOfThreads; ++j) {
		singleThreadPtr = boost::shared_ptr<boost::thread>(new boost::thread(boost::bind(ThreadPoolWorker, this)));
		threads.push_back(singleThreadPtr);
	}
	
	// test if the threads have finished
	for (;;) {
		firstThreadHasFinished = threads.at(0)->timed_join(boost::posix_time::milliseconds(500));
		if (firstThreadHasFinished == 0) {	// the thread is not done yet, we're just waiting
											// while we wait we check for various things and give the interface the chance to update
			
			// did one of the threads run into an error?
			this->errorReportingMutex.lock();
			if (this->errorMessage.length() != 0) {
				errorReportingMutex.unlock();
				// an error occured, time to abort this analysis
				for (size_t j = 0; j < numberOfThreads; ++j) {
					threads.at(j)->interrupt();
				}
				// wait until the threads have completed
				for (size_t j = 0; j < numberOfThreads; ++j) {
				threads.at(j)->join();
				}
				
				throw ERROR_RUNNING_THREADED_ANALYSIS(this->errorMessage);
			}
			this->errorReportingMutex.unlock();
			
			// does the user want to abort?
#ifdef WITH_IGOR
			status = CheckAbort(0);
			if (status == -1) {
				for (size_t j = 0; j < numberOfThreads; ++j) {
					threads.at(j)->interrupt();
				}
				// wait until the threads have completed
				for (size_t j = 0; j < numberOfThreads; ++j) {
					threads.at(j)->join();
				}
				progressReporter->CalculationAborted();
				this->localizedPositions->sortPositionsByFrameNumber();
				return this->localizedPositions;
			}
#endif // WITH_IGOR
			
			// allow the reporter to update with new progress
			this->acquireFrameForProcessingMutex.lock();
			percentDone = (double)(nImages - framesToBeProcessed.size()) / (double)(this->nImages) * 100.0;
			progressReporter->UpdateCalculationProgress(percentDone);
			this->acquireFrameForProcessingMutex.unlock();
			continue;
		} else {
			break;
		}
	}
	
	for (size_t j = 1; j < numberOfThreads; ++j) {
		threads.at(j)->join();
	}
	
	// the processing itself is now done, but the results will not have been returned in the correct order
	this->localizedPositions->sortPositionsByFrameNumber();
	
	progressReporter->CalculationDone();
	
	// return the results
	return this->localizedPositions;
}

void ThreadPoolWorker(PALMAnalysisController* controller) {
	size_t currentImageToProcess;
	boost::shared_ptr<ublas::matrix<double> > currentImage;
	boost::shared_ptr<ublas::matrix <unsigned char> > thresholdedImage;
	boost::shared_ptr<std::list<position> > locatedParticles;
	boost::shared_ptr<LocalizedPositionsContainer> localizedPositions;
	
	try {
		for (;;) {	// loop continuously looking for more images until there are none left
			// if the main thread wants us to interrupt, then give it an opportunity to do so
			if (boost::this_thread::interruption_requested()) {
				return;
			}
			
			// start by a acquiring an image to process
			controller->acquireFrameForProcessingMutex.lock();
			if (controller->framesToBeProcessed.size() == 0) {
				// no more frames to be processed
				controller->acquireFrameForProcessingMutex.unlock();
				return;
			}
			currentImageToProcess = controller->framesToBeProcessed.front();
			controller->framesToBeProcessed.pop();
			controller->acquireFrameForProcessingMutex.unlock();
			
			// we need to process the image with index currentImageToProcess
			currentImage = controller->imageLoader->get_nth_image(currentImageToProcess);
			
			thresholdedImage = do_processing_and_thresholding(currentImage, controller->thresholdImagePreprocessor, controller->thresholder,
															  controller->thresholdImagePostprocessor);
			
			locatedParticles = controller->particleFinder->findPositions(currentImage, thresholdedImage);
			
			// if the located particles are to be verified before fitting then do so
			for (std::vector<boost::shared_ptr<ParticleVerifier> >::iterator it = controller->particleVerifiers.begin(); it != controller->particleVerifiers.end(); ++it) {
				(*it)->VerifyParticles(currentImage, locatedParticles);
			}
			
			localizedPositions = controller->fitPositions->fit_positions(currentImage, locatedParticles);
			
			// the localization routines have no idea what frame number the processed frame was
			// so we need to add the appropriate frame number to the fitted positions here
			localizedPositions->setFrameNumbers(currentImageToProcess);
			
			// pass the result to the output queue
			controller->addLocalizedPositionsMutex.lock();
			// if this is the first time that positions are being returned then controller will contain a NULL pointer
			// set it to the positions we are now returning
			// this way ThreadPoolWorker does not need to know what the type of positions is
			if (controller->localizedPositions.get() == NULL) {
				controller->localizedPositions = localizedPositions;
			} else {
				controller->localizedPositions->addPositions(localizedPositions);
			}
			
			controller->addLocalizedPositionsMutex.unlock();
		}
	}
	catch (std::runtime_error &e) {
		// an error has appeared somewhere in this thread during the calculation
		// since there seems to be no easy way to communicate the exception to the
		// main thread, set an error message in the analysis controller.
		// The controller will periodically check this message and handle the error
		controller->errorReportingMutex.lock();
		controller->errorMessage.assign(e.what());
		controller->errorReportingMutex.unlock();
		
		// no point in continuing this thread
		return;
	}
	catch (boost::thread_interrupted) {
		// the main thread wants us to stop
		return;
	}
	catch (...) {
		// catch any other exception not handled by the above block
		controller->errorReportingMutex.lock();
		controller->errorMessage.assign("Encountered an unspecified exception while doing the PALM analysis");
		controller->errorReportingMutex.unlock();
		
		// no point in continuing this thread
		return;
	}
	
}

#ifdef WITH_IGOR
void PALMAnalysisProgressReporter_IgorCommandLine::UpdateCalculationProgress(double percentDone) {
	char XOPOut[10];
	if (percentDone - previousPercentage > 10.0) {
		previousPercentage = floor(percentDone / 10.0) * 10.0;
		sprintf(XOPOut, "%.0lf%% ", previousPercentage);
		XOPNotice(XOPOut);
	}
}
#endif // WITH_IGOR

boost::shared_ptr<LocalizedPositionsContainer> FitPositionsDeflate::fit_positions(const boost::shared_ptr<ublas::matrix<double> > image, boost::shared_ptr<std::list<position> > positions) {
	// TODO: for now we ignore the starting positions and ending position provided as arguments
	
	boost::shared_ptr<LocalizedPositionsContainer> positionsFittedThusFar;
	boost::shared_ptr<LocalizedPositionsContainer> positionsLocalizedThisFrame;
	boost::shared_ptr<ublas::matrix<double> > subtractedImage;
	boost::shared_ptr<ublas::matrix <unsigned char> > segmentedImage;
	boost::shared_ptr<std::list<position> > locatedParticles;
	
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


boost::shared_ptr<ublas::matrix<double> > FitPositionsDeflate::subtractLocalizedPositions(boost::shared_ptr<ublas::matrix<double> > image, boost::shared_ptr<LocalizedPositionsContainer> positions) {
	double fittedXPos, fittedYPos, fittedIntegral, fittedXWidth, fittedYWidth;
	double centerX, centerY, calculatedAmplitude;
	size_t startX, endX, startY, endY;
	double distanceXSquared, distanceYSquared, currentIntensity;
	
	size_t nPositions = positions->getNPositions();
	size_t xSize = image->size1();
	size_t ySize = image->size2();
	boost::shared_ptr<ublas::matrix<double> > outputImage(new ublas::matrix<double> (xSize, ySize));
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
		
		for (size_t i = startX; i <= endX; ++i) {
			for (size_t j = startY; j < endY; ++j) {
				distanceXSquared = ((double)i - (double)centerX) * ((double)i - (double)centerX);
				distanceYSquared = ((double)j - (double)centerY) * ((double)j - (double)centerY);
				currentIntensity = calculatedAmplitude * exp(- (distanceXSquared + distanceYSquared) / (2 * fittedXWidth * fittedYWidth));
				
				(*outputImage)(i, j) = (*outputImage)(i, j) - currentIntensity;
			}
		}
	}
	return outputImage;
}
