#include "PALM_analysis.h"


int load_partial_ccd_image(ImageLoader *image_loader, size_t n_start, size_t n_end) {
	size_t total_n_images = image_loader->get_total_number_of_images();
	size_t n_images_to_load;
	size_t x_size, y_size;
	waveHndl output_wave;
	long dimension_sizes[MAX_DIMENSIONS + 1];
	long indices[MAX_DIMENSIONS];
	
	int result;
	boost::shared_ptr<PALMMatrix<double> > current_image;
	double current_value;
	double current_value_array[2];
	
	if (n_start > n_end)
		throw END_SHOULD_BE_LARGER_THAN_START();
	
	// how many images should we load?
	if (n_end <= total_n_images) {
		n_images_to_load = n_end - n_start + 1;
	} else {
		n_images_to_load = total_n_images - n_start;
		n_end = total_n_images - 1;
	}
	
	x_size = image_loader->getXSize();
	y_size = image_loader->getYSize();
	
	SetIgorIntVar("V_numberOfImages", (long)total_n_images, 1);
	SetIgorIntVar("V_xSize", (long)x_size, 1);
	SetIgorIntVar("V_ySize", (long)y_size, 1);
	SetIgorIntVar("V_firstImageLoaded", (long)n_start, 1);
	SetIgorIntVar("V_lastImageLoaded", (long)n_end, 1);
	
	
	// make the output wave that will store the data
	dimension_sizes[0] = (long)x_size;
	dimension_sizes[1] = (long)y_size;
	dimension_sizes[2] = (long)n_images_to_load;
	dimension_sizes[3] = 0;
	
	
	result = MDMakeWave(&output_wave, "M_CCDFrames", NULL, dimension_sizes, NT_FP32, 1);
	
	if (result != 0) {
		throw result;
	}
	
	// load the data
	for (size_t i = n_start; i <= n_end; i++) {
		current_image = image_loader->get_nth_image(i);
		
		indices[2] = i - n_start;
		
		// store the data in the output wave
		for (size_t k = 0; k < y_size; k++) {
			for (size_t j = 0; j < x_size; j++) {
				current_value = (*current_image)(j, k);
				
				indices[0] = j;
				indices[1] = k;
				
				current_value_array[0] = current_value;
				
				result = MDSetNumericWavePointValue(output_wave, indices, current_value_array);
				if (result != 0) {
					throw result;
				}
			}
		}
		
	}
	
	
	return 0;
}


int parse_ccd_headers(ImageLoader *image_loader) {
	size_t total_n_images = image_loader->get_total_number_of_images();
	size_t x_size = image_loader->getXSize();
	size_t y_size = image_loader->getYSize();
	
	int result;
	
	result = SetIgorIntVar("V_numberOfImages", (long)total_n_images, 1);
	if (result != 0)
		return result;
	result = SetIgorIntVar("V_xSize", (long)x_size, 1);
	if (result != 0)
		return result;
	result = SetIgorIntVar("V_ySize", (long)y_size, 1);
	if (result != 0)
		return result;
	return 0;
}


boost::shared_ptr<PALMMatrix <unsigned char> > do_processing_and_thresholding(boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<ThresholdImage_Preprocessor>preprocessor, 
																		 boost::shared_ptr<ThresholdImage> thresholder, boost::shared_ptr<ThresholdImage_Postprocessor> postprocessor) {
	// this function takes care of the thresholding and the associated pre- and postprocessing
	
	boost::shared_ptr<PALMMatrix<double> > preprocessed_image;
	boost::shared_ptr<PALMMatrix <unsigned char> > thresholded_image;
	boost::shared_ptr<PALMMatrix <unsigned char> > postprocessed_image;
	
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

int construct_summed_intensity_trace(ImageLoader *image_loader, string output_wave_name, long startX, long startY, long endX, long endY) {
	size_t n_images = image_loader->get_total_number_of_images();
	size_t x_size = image_loader->getXSize();
	size_t y_size = image_loader->getYSize();
	
	boost::shared_ptr<PALMMatrix<double> > current_image;
	double summed_intensity;
	
	waveHndl output_wave;
	long dimension_sizes[MAX_DIMENSIONS + 1];
	dimension_sizes[0] = (long)n_images;
	dimension_sizes[1] = 0;
	int result;
	
	// check if we want the full frame
	if (startX == -1) {
		startX = 0;
		startY = 0;
		endX = x_size - 1;
		endY = y_size - 1;
	} else {
		if ((startX >= x_size) || (endX >= x_size) || (startY >= y_size) || (endY >= y_size)) {
			throw kBadROIDimensions;
		}
	}
	
	// try to allocate a buffer that will hold the intensity trace
	boost::scoped_array<double> intensity_trace_buffer(new double[n_images]);
	
	for (size_t i = 0; i < n_images; i++) {
		summed_intensity = 0;
		current_image = image_loader->get_nth_image(i);
		
		// calculate the total sum of the image
		for (size_t k = startY; k <= endY; k++) {
			for (size_t j = startX; j <= endX; j++) {
				summed_intensity += (*current_image)(j, k);
			}
		}
		// store the contents in the buffer
		intensity_trace_buffer[i] = summed_intensity;
	}
	
	// try to create the output wave
	result = MDMakeWave(&output_wave, output_wave_name.c_str(), NULL, dimension_sizes, NT_FP64, 1);
	if (result != 0) {
		throw result;
	}
	
	// write the output data to the wave
	result = MDStoreDPDataInNumericWave(output_wave, intensity_trace_buffer.get());
	if (result != 0) {
		throw result;
	}
	
	return 0;
}

PALMAnalysisController::PALMAnalysisController(boost::shared_ptr<ImageLoader> imageLoader_rhs, boost::shared_ptr<ThresholdImage> thresholder_rhs,
											   boost::shared_ptr<ThresholdImage_Preprocessor> thresholdImagePreprocessor_rhs,
											   boost::shared_ptr<ThresholdImage_Postprocessor> thresholdImagePostprocessor_rhs,
											   boost::shared_ptr<ParticleFinder> particleFinder_rhs, boost::shared_ptr<FitPositions> fitPositions_rhs,
											   boost::shared_ptr<PALMResultsWriter> resultsWriter_rhs,
											   boost::shared_ptr<PALMAnalysisProgressReporter> progressReporter_rhs) {
	imageLoader = imageLoader_rhs;
	thresholder = thresholder_rhs;
	thresholdImagePreprocessor = thresholdImagePreprocessor_rhs;
	thresholdImagePostprocessor = thresholdImagePostprocessor_rhs;
	particleFinder = particleFinder_rhs;
	fitPositions = fitPositions_rhs;
	resultsWriter = resultsWriter_rhs;
	progressReporter = progressReporter_rhs;
	
	nImages = imageLoader->get_total_number_of_images();
}

void PALMAnalysisController::DoPALMAnalysis() {
	size_t numberOfProcessors = boost::thread::hardware_concurrency();
	size_t numberOfThreads;
	vector<boost::shared_ptr<boost::thread> > threads;
	boost::shared_ptr<boost::thread> singleThreadPtr;
	std::list<boost::shared_ptr<PALMResults> >::iterator it;
	int firstThreadHasFinished, status;
	double percentDone;
	
	numberOfThreads = numberOfProcessors * 2;	// take two threads for every processor since every thread will be blocked on I/O sooner or later
	if (numberOfThreads == 0) {
		numberOfThreads = 1;
	}
	if (numberOfThreads > this->nImages) {
		numberOfThreads = nImages;
	}
	
	this->fittedPositionsList.clear();
	
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
											// while we wait we check for a user abort and give the interface the chance to update
			status = CheckAbort(0);
			if (status == -1) {
				for (size_t j = 0; j < numberOfThreads; ++j) {
					threads.at(j)->interrupt();
				}
				progressReporter->CalculationAborted();
				return;
			}
			this->addPALMResultsMutex.lock();
			percentDone = (double)(fittedPositionsList.size()) / (double)(this->nImages) * 100.0;
			progressReporter->UpdateCalculationProgress(percentDone);
			this->addPALMResultsMutex.unlock();
			continue;
		} else {
			break;
		}
	}
	
	for (size_t j = 1; j < numberOfThreads; ++j) {
		threads.at(j)->join();
	}
	
	// the processing itself is now done, but the results will not have been returned in the correct order
	fittedPositionsList.sort(ComparePALMResults);
	
	// store the results
	for (it = this->fittedPositionsList.begin(); it != this->fittedPositionsList.end(); ++it) {
		this->resultsWriter->AppendNewResult(*it);
	}
	
	progressReporter->CalculationDone();
}

void ThreadPoolWorker(PALMAnalysisController* controller) {
	size_t currentImageToProcess;
	boost::shared_ptr<PALMMatrix<double> > currentImage;
	boost::shared_ptr<PALMMatrix <unsigned char> > thresholdedImage;
	boost::shared_ptr<PALMMatrix<double> > locatedParticles;
	boost::shared_ptr<PALMMatrix<double> > fittedPositions;
	boost::shared_ptr<PALMResults> analysisResult;
	
	for (;;) {	// loop continuously looking for more images until there are none left
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
		fittedPositions = controller->fitPositions->fit_positions(currentImage, locatedParticles);
		
		analysisResult = boost::shared_ptr<PALMResults> (new PALMResults(currentImageToProcess, fittedPositions));
		
		// pass the result to the output queue
		controller->addPALMResultsMutex.lock();
		controller->fittedPositionsList.push_back(analysisResult);
		controller->addPALMResultsMutex.unlock();
	}
	
}

int ComparePALMResults(boost::shared_ptr<PALMResults> result1, boost::shared_ptr<PALMResults> result2) {
	assert(result1->getFrameIndex() != result2->getFrameIndex());
	if (result1->getFrameIndex() < result2->getFrameIndex()) {
		return 1;
	} else {
		return 0;
	}
}

void PALMAnalysisProgressReporter_IgorCommandLine::UpdateCalculationProgress(double percentDone) {
	char XOPOut[10];
	if (percentDone - previousPercentage > 10.0) {
		previousPercentage = floor(percentDone / 10.0) * 10.0;
		sprintf(XOPOut, "%.0lf%% ", previousPercentage);
		XOPNotice(XOPOut);
	}
}


boost::shared_ptr<PALMVolume <unsigned short> > calculate_PALM_bitmap_image(boost::shared_ptr<PALMMatrix<double> > positions, boost::shared_ptr<PALMMatrix<double> > colors, boost::shared_ptr<PALMBitmapImageDeviationCalculator> deviationCalculator,
																size_t xSize, size_t ySize, size_t imageWidth, size_t imageHeight, int normalizeColors) {
	size_t nPositions = positions->getXSize();
	size_t nColors = colors->getXSize();
	size_t nFrames = (size_t)(positions->get(nPositions - 1, 0)) + 1;
	size_t currentFrame;
	
	double imageScaleFactor = (double)(imageWidth - 1) / (double)xSize;
	double maxAmplitude = 0;
	
	double currentX, currentY, currentAmplitude;
	size_t centerX, centerY;
	long startX, endX, startY, endY;
	double deviation, currentIntensity;
	double distanceXSquared, distanceYSquared;
	size_t colorIndex;
	double currentColors[3];
	double summedIntensity;
	
	boost::shared_ptr<PALMVolume <unsigned short> > outputImage(new PALMVolume <unsigned short>(imageWidth, imageHeight, 3));	// 3 layers because it will be a direct color image
	boost::shared_ptr<PALMMatrix<double> > totalIntensities(new PALMMatrix<double>(imageWidth, imageHeight));	// keep track of the total intensities in each pixel
	
	outputImage->set_all(0);
	totalIntensities->set_all(0);
	
	if (normalizeColors != 0) {
		// get the position with the maximum amplitude
		// those positions will have the 'full' colors, the colors of the other positions will be scaled relative to it
		for (size_t n = 0; n < nPositions; ++n) {
			if ((*positions)(n, 1) > maxAmplitude) {
				maxAmplitude = positions->get(n,1);
			}
		}
	}
	
	for (size_t n = 0; n < nPositions; ++n) {
		currentFrame = (size_t)((*positions)(n, 0) + 0.5);
		currentAmplitude = positions->get(n, 1);
		currentX = positions->get(n, 3);
		currentY = positions->get(n, 4);
		
		if ((currentAmplitude < 0) || (currentX < 0) || (currentX >= xSize) || (currentY < 0) || (currentY >= ySize)) {
			continue;
		}
		
		if (normalizeColors == 0) {
			currentAmplitude = 1.0;	// every position is equally important when we don't do scaling
		}
		
		centerX = (size_t)(currentX * imageScaleFactor + 0.5);
		centerY = (size_t)(currentY * imageScaleFactor + 0.5);
		deviation = (deviationCalculator->getDeviation(positions, n) * imageScaleFactor);
		
		startX = floor((double)centerX - 4.0 * deviation);
		startY = floor((double)centerY - 4.0 * deviation);
		endX = ceil((double)centerX + 4.0 * deviation);
		endY = ceil((double)centerY + 4.0 * deviation);
		
		if (startX < 0)
			startX = 0;
		if (endX >= imageWidth)
			endX = imageWidth - 1;
		if (startY < 0)
			startY = 0;
		if (endY >= imageHeight)
			endY = imageHeight - 1;
		
		colorIndex = (size_t)((double)currentFrame / (double)nFrames * (double)(nColors - 1) + 0.5);
		
		for (size_t i = startX; i <= endX; ++i) {
			for (size_t j = startY; j < endY; ++j) {
				distanceXSquared = ((double)i - (double)centerX) * ((double)i - (double)centerX);
				distanceYSquared = ((double)j - (double)centerY) * ((double)j - (double)centerY);
				currentIntensity = currentAmplitude * exp(- (distanceXSquared + distanceYSquared) / (2 * deviation * deviation));
				
				if (normalizeColors != 0) {
					currentColors[0] = (*colors)(colorIndex, 0) * currentIntensity / maxAmplitude;	// Simplification of colors->get(colorIndex, 0) * currentIntensity / currentAmplitude * currentAmplitude / maxAmplitude
					currentColors[1] = (*colors)(colorIndex, 1) * currentIntensity / maxAmplitude;
					currentColors[2] = (*colors)(colorIndex, 2) * currentIntensity / maxAmplitude;
				} else {
					currentColors[0] = (*colors)(colorIndex, 0) * currentIntensity / currentAmplitude;
					currentColors[1] = (*colors)(colorIndex, 1) * currentIntensity / currentAmplitude;
					currentColors[2] = (*colors)(colorIndex, 2) * currentIntensity / currentAmplitude;
				}
				
				summedIntensity = currentIntensity + (*totalIntensities)(i, j);
				
				(*outputImage)(i, j, 0) = (unsigned short)(currentIntensity / summedIntensity * currentColors[0] + totalIntensities->get(i, j) / summedIntensity * outputImage->get(i, j, 0));
				(*outputImage)(i, j, 1) = (unsigned short)(currentIntensity / summedIntensity * currentColors[1] + totalIntensities->get(i, j) / summedIntensity * outputImage->get(i, j, 1));
				(*outputImage)(i, j, 2) = (unsigned short)(currentIntensity / summedIntensity * currentColors[2] + totalIntensities->get(i, j) / summedIntensity * outputImage->get(i, j, 2));
				
				(*totalIntensities)(i,j) = totalIntensities->get(i, j) + currentIntensity;
			}
		}
	}
	
	return outputImage;
}

boost::shared_ptr<PALMVolume <unsigned short> > calculate_PALM_bitmap_image_parallel(boost::shared_ptr<PALMMatrix<double> > positions, boost::shared_ptr<PALMMatrix<double> > colors, boost::shared_ptr<PALMBitmapImageDeviationCalculator> deviationCalculator,
																		 size_t xSize, size_t ySize, size_t imageWidth, size_t imageHeight, int normalizeColors) {
	vector<boost::shared_ptr<boost::thread> > threads;
	boost::shared_ptr<boost::thread> singleThreadPtr;
	vector<size_t> startPositions;
	vector<size_t> endPositions;
	vector<boost::shared_ptr<calculate_PALM_bitmap_image_ThreadStartParameters> > threadData;
	boost::shared_ptr<PALMVolume <unsigned short> > outputImage;
	boost::shared_ptr<PALMMatrix<double> > totalIntensities;
	
	size_t nPositions = positions->getXSize();
	size_t numberOfProcessors = boost::thread::hardware_concurrency();
	if (numberOfProcessors == 0) {
		numberOfProcessors = 1;
	}
	size_t nThreads = numberOfProcessors;
	size_t nPositionsPerThread;
	size_t currentThreadStart;
	size_t currentThreadEnd;
	size_t nFrames = (size_t)(positions->get(nPositions - 1, 0)) + 1;
	
	double maxAmplitude = 0;
	double imageScaleFactor = (double)(imageWidth - 1) / (double)xSize;
	double summedIntensity;
	
	nThreads = nPositions / 5;
	if (nThreads > numberOfProcessors) {
		nThreads = numberOfProcessors;
	}
	
	if (nThreads == 0) {
		nThreads = 1;
	}
	
	nPositionsPerThread = nPositions / nThreads;
	
	if (normalizeColors != 0) {
		// get the position with the maximum amplitude
		// those positions will have the 'full' colors, the colors of the other positions will be scaled relative to it
		for (size_t n = 0; n < nPositions; ++n) {
			if ((*positions)(n,1) > maxAmplitude) {
				maxAmplitude = (*positions)(n,1);
			}
		}
	}
	
	// set up the vector containing the data for the threads
	threadData.clear();
	currentThreadEnd = (size_t)-1;
	for (size_t i = 0; i < nThreads; ++i) {
		if (i == nThreads - 1) {	// the last thread is special
			currentThreadStart = currentThreadEnd + 1;
			currentThreadEnd = nPositions - 1;
		} else {
			currentThreadStart = currentThreadEnd + 1;
			currentThreadEnd = currentThreadStart + nPositionsPerThread;
		}
		
		boost::shared_ptr<calculate_PALM_bitmap_image_ThreadStartParameters> singleThreadStartParameter(new calculate_PALM_bitmap_image_ThreadStartParameters);
		singleThreadStartParameter->positions = positions;
		singleThreadStartParameter->image = boost::shared_ptr<PALMVolume <unsigned short> > (new PALMVolume <unsigned short>(imageWidth, imageHeight, 3));
		singleThreadStartParameter->totalIntensities = boost::shared_ptr<PALMMatrix<double> > (new PALMMatrix<double>(imageWidth, imageHeight));
		singleThreadStartParameter->colors = colors;
		singleThreadStartParameter->deviationCalculator = deviationCalculator;
		singleThreadStartParameter->normalizeColors = normalizeColors;
		singleThreadStartParameter->nFrames = nFrames;
		singleThreadStartParameter->startIndex = currentThreadStart;
		singleThreadStartParameter->endIndex = currentThreadEnd;
		singleThreadStartParameter->imageWidth = imageWidth;
		singleThreadStartParameter->imageHeight = imageHeight;
		singleThreadStartParameter->xSize = xSize;
		singleThreadStartParameter->ySize = ySize;
		singleThreadStartParameter->maxAmplitude = maxAmplitude;
		singleThreadStartParameter->scaleFactor = imageScaleFactor;
		
		threadData.push_back(singleThreadStartParameter);
	}
	
	// now start the threads
	threads.clear();
	for (size_t j = 0; j < nThreads; ++j) {
		singleThreadPtr = boost::shared_ptr<boost::thread>(new boost::thread(&calculate_PALM_bitmap_image_ThreadStart, threadData.at(j)));
		threads.push_back(singleThreadPtr);
	}
	
	// wait for the threads to finish
	for (size_t j = 0; j < nThreads; ++j) {
		threads.at(j)->join();
	}
	
	// combine the individual images into the final output image
	outputImage = threadData.at(0)->image;
	totalIntensities = threadData.at(0)->totalIntensities;
	
	for (size_t n = 1; n < nThreads; ++n) {
		for (size_t i = 0; i < imageWidth; ++i) {
			for (size_t j = 0; j < imageHeight; ++j) {
				summedIntensity = totalIntensities->get(i,j) + threadData[n]->totalIntensities->get(i,j);
				(*outputImage)(i, j, 0) = (unsigned short)((*totalIntensities)(i,j) / summedIntensity * (*outputImage)(i, j, 0) + threadData[n]->totalIntensities->get(i,j) / summedIntensity * threadData[n]->image->get(i, j, 0));
				(*outputImage)(i, j, 1) = (unsigned short)((*totalIntensities)(i,j) / summedIntensity * (*outputImage)(i, j, 1) + threadData[n]->totalIntensities->get(i,j) / summedIntensity * threadData[n]->image->get(i, j, 0));
				(*outputImage)(i, j, 2) = (unsigned short)((*totalIntensities)(i,j) / summedIntensity * (*outputImage)(i, j, 2) + threadData[n]->totalIntensities->get(i,j) / summedIntensity * threadData[n]->image->get(i, j, 0));
				(*totalIntensities)(i,j) = summedIntensity;
			}
		}
	}
	
	return outputImage;
}


void calculate_PALM_bitmap_image_ThreadStart(boost::shared_ptr<calculate_PALM_bitmap_image_ThreadStartParameters> startParameters) {
	boost::shared_ptr<PALMMatrix<double> > positions = startParameters->positions;
	boost::shared_ptr<PALMVolume <unsigned short> > image = startParameters->image;
	boost::shared_ptr<PALMMatrix<double> > totalIntensities = startParameters->totalIntensities;
	boost::shared_ptr<PALMMatrix<double> > colors = startParameters->colors;
	
	boost::shared_ptr<PALMBitmapImageDeviationCalculator> deviationCalculator = startParameters->deviationCalculator;
	
	size_t nColors = colors->getXSize();
	size_t startIndex = startParameters->startIndex;
	size_t endIndex = startParameters->endIndex;
	size_t nFrames = startParameters->nFrames;
	size_t currentFrame;
	int normalizeColors = startParameters->normalizeColors;
	
	size_t imageWidth = startParameters->imageWidth;
	size_t imageHeight = startParameters->imageHeight;
	size_t xSize = startParameters->xSize;
	size_t ySize = startParameters->ySize;
	double imageScaleFactor = startParameters->scaleFactor;
	double maxAmplitude = startParameters->maxAmplitude;
	
	double currentX, currentY, currentAmplitude;
	size_t centerX, centerY;
	long startX, endX, startY, endY;
	double deviation, currentIntensity;
	double distanceXSquared, distanceYSquared;
	size_t colorIndex;
	double currentColors[3];
	double summedIntensity;
	
	image->set_all(0);
	totalIntensities->set_all(0);
	
	for (size_t n = startIndex; n <= endIndex; ++n) {
		currentFrame = (size_t)((*positions)(n, 0) + 0.5);
		currentAmplitude = (*positions)(n, 1);
		currentX = (*positions)(n, 3);
		currentY = (*positions)(n, 4);
		
		if ((currentAmplitude < 0) || (currentX < 0) || (currentX >= xSize) || (currentY < 0) || (currentY >= ySize)) {
			continue;
		}
		
		if (normalizeColors == 0) {
			currentAmplitude = 1.0;	// every position is equally important when we don't do scaling
		}
		
		centerX = (size_t)(currentX * imageScaleFactor + 0.5);
		centerY = (size_t)(currentY * imageScaleFactor + 0.5);
		deviation = (deviationCalculator->getDeviation(positions, n) * imageScaleFactor);
		
		startX = floor((double)centerX - 4.0 * deviation);
		startY = floor((double)centerY - 4.0 * deviation);
		endX = ceil((double)centerX + 4.0 * deviation);
		endY = ceil((double)centerY + 4.0 * deviation);
		
		if (startX < 0)
			startX = 0;
		if (endX >= imageWidth)
			endX = imageWidth - 1;
		if (startY < 0)
			startY = 0;
		if (endY >= imageHeight)
			endY = imageHeight - 1;
		
		colorIndex = (size_t)((double)currentFrame / (double)nFrames * (double)(nColors - 1) + 0.5);
		
		for (size_t i = startX; i <= endX; ++i) {
			for (size_t j = startY; j < endY; ++j) {
				distanceXSquared = ((double)i - (double)centerX) * ((double)i - (double)centerX);
				distanceYSquared = ((double)j - (double)centerY) * ((double)j - (double)centerY);
				currentIntensity = currentAmplitude * exp(- (distanceXSquared + distanceYSquared) / (2 * deviation * deviation));
				
				if (normalizeColors != 0) {
					currentColors[0] = colors->get(colorIndex, 0) * currentIntensity / maxAmplitude;	// Simplification of colors->get(colorIndex, 0) * currentIntensity / currentAmplitude * currentAmplitude / maxAmplitude
					currentColors[1] = colors->get(colorIndex, 1) * currentIntensity / maxAmplitude;
					currentColors[2] = colors->get(colorIndex, 2) * currentIntensity / maxAmplitude;
				} else {
					currentColors[0] = colors->get(colorIndex, 0) * currentIntensity / currentAmplitude;
					currentColors[1] = colors->get(colorIndex, 1) * currentIntensity / currentAmplitude;
					currentColors[2] = colors->get(colorIndex, 2) * currentIntensity / currentAmplitude;
				}
				
				summedIntensity = currentIntensity + totalIntensities->get(i, j);
				
				image->set(i, j, 0, (unsigned short)(currentIntensity / summedIntensity * currentColors[0] + totalIntensities->get(i, j) / summedIntensity * image->get(i, j, 0)));
				image->set(i, j, 1, (unsigned short)(currentIntensity / summedIntensity * currentColors[1] + totalIntensities->get(i, j) / summedIntensity * image->get(i, j, 1)));
				image->set(i, j, 2, (unsigned short)(currentIntensity / summedIntensity * currentColors[2] + totalIntensities->get(i, j) / summedIntensity * image->get(i, j, 2)));
				
				totalIntensities->set(i,j, (totalIntensities->get(i, j) + currentIntensity));
			}
		}
	}
}


int construct_average_image(ImageLoader *image_loader, string output_wave_name, long startX, long startY, long endX, long endY) {
	size_t n_images = image_loader->get_total_number_of_images();
	size_t x_size = image_loader->getXSize();
	size_t y_size = image_loader->getYSize();
	
	long xRange, yRange;
	
	if (endX > x_size - 1)
		endX = x_size - 1;
	if (endY > y_size - 1)
		endY = y_size - 1;
	
	if ((startX < 0) || (startY < 0) || (startY < 0) || (startY < 0)) {	// we want the full frame
		startX = 0;
		startY = 0;
		endX = x_size - 1;
		endY = y_size - 1;
	}
	
	xRange = endX - startX + 1;
	yRange = endY - startY + 1;
	
	boost::shared_ptr<PALMMatrix<double> > current_image;
	boost::shared_ptr<PALMMatrix<double> > average_image(new PALMMatrix<double>(xRange, yRange));
	
	waveHndl output_wave;
	long dimension_sizes[MAX_DIMENSIONS + 1];
	dimension_sizes[0] = xRange;
	dimension_sizes[1] = yRange;
	dimension_sizes[2] = 0;
	long indices[MAX_DIMENSIONS];
	double current_value[2];
	int result;
	
	average_image->set_all(0);
	
	
	for (size_t i = 0; i < n_images; i++) {
		current_image = image_loader->get_nth_image(i);
		
		// add the values of the newly loaded image to the average image
		(*average_image) += (*current_image);
	}
	
	// divide by the number of images
	(*average_image) = (*average_image).DivideByScalar(n_images);
	
	// try to create the output wave
	result = MDMakeWave(&output_wave, output_wave_name.c_str(), NULL, dimension_sizes, NT_FP64, 1);
	if (result != 0) {
		throw result;
	}
	
	// write the output data to the wave
	for (size_t i = startX; i <= endX; ++i) {
		for (size_t j = startY; j <= endY; ++j) {
			current_value[0] = (*average_image)(i, j);
			indices[0] = i;
			indices[1] = j;
			result = MDSetNumericWavePointValue(output_wave, indices, current_value);
			if (result != 0) {
				throw result;
			}
		}
	}
	
	return 0;
}


void calculateStandardDeviationImage(ImageLoader *image_loader, string output_wave_name, long startX, long startY, long endX, long endY) {
	size_t n_images = image_loader->get_total_number_of_images();
	size_t x_size = image_loader->getXSize();
	size_t y_size = image_loader->getYSize();
	int result;
	waveHndl output_wave;
	long dimension_sizes[MAX_DIMENSIONS + 1];
	long indices[MAX_DIMENSIONS];
	
	double current_value[2];
	
	long xRange, yRange;
	
	if (endX > x_size - 1)
		endX = x_size - 1;
	if (endY > y_size - 1)
		endY = y_size - 1;
	
	if ((startX < 0) || (startY < 0) || (startY < 0) || (startY < 0)) {	// we want the full frame
		startX = 0;
		startY = 0;
		endX = x_size - 1;
		endY = y_size - 1;
	}
	
	xRange = endX - startX + 1;
	yRange = endY - startY + 1;
	
	boost::scoped_ptr<PALMMatrix<double> > stdDevImage(new PALMMatrix<double>(xRange, yRange));
	boost::scoped_ptr<PALMMatrix<double> > average_image(new PALMMatrix<double>(xRange, yRange));
	boost::shared_ptr<PALMMatrix<double> > current_image;
	
	average_image->set_all(0);
	stdDevImage->set_all(0);
	
	// construct an average image
	for (size_t i = 0; i < n_images; i++) {
		current_image = image_loader->get_nth_image(i);
		
		// add the values of the newly loaded image to the average image
		(*average_image) += (*current_image);
	}
	
	// divide by the number of images
	(*average_image) = (*average_image).DivideByScalar(n_images);
	
	// now loop over the images again, calculating the standard deviation of each pixel
	for (size_t i = 0; i < n_images; i++) {
		current_image = image_loader->get_nth_image(i);
		
		// add the deviation of the newly loaded image from the mean to the stddev image
		(*stdDevImage) += (((*current_image) - (*average_image)) * ((*current_image) - (*average_image)));
	}
	
	// divide by the number of images to get the average deviation, and take the square root
	(*stdDevImage) = (*stdDevImage).DivideByScalar(n_images);
	(*stdDevImage) = (*stdDevImage).RaiseToPower(0.5);
	
	// try to create the output wave
	dimension_sizes[0] = xRange;
	dimension_sizes[1] = yRange;
	dimension_sizes[2] = 0;
	result = MDMakeWave(&output_wave, output_wave_name.c_str(), NULL, dimension_sizes, NT_FP64, 1);
	if (result != 0) {
		throw result;
	}
	
	// write the output data to the wave
	for (size_t i = startX; i <= endX; ++i) {
		for (size_t j = startY; j <= endX; ++j) {
			current_value[0] = (*stdDevImage)(i, j);
			indices[0] = i;
			indices[1] = j;
			result = MDSetNumericWavePointValue(output_wave, indices, current_value);
			if (result != 0) {
				throw result;
			}
		}
	}
}
	