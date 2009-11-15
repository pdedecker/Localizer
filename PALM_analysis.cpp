#include "PALM_analysis.h"


int load_partial_ccd_image(ImageLoader *image_loader, size_t n_start, size_t n_end) {
	size_t total_n_images = image_loader->get_total_number_of_images();
	size_t n_images_to_load;
	size_t x_size, y_size;
	int storage_type, waveType;
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
	storage_type = image_loader->getStorageType();
	
	result = SetOperationNumVar("V_numberOfImages", total_n_images);
	if (result != 0)
		throw result;
	result = SetOperationNumVar("V_xSize", x_size);
	if (result != 0)
		throw result;
	result = SetOperationNumVar("V_ySize", y_size);
	if (result != 0)
		throw result;
	result = SetOperationNumVar("V_firstImageLoaded", n_start);
	if (result != 0)
		throw result;
	result = SetOperationNumVar("V_lastImageLoaded", n_end);
	if (result != 0)
		throw result;
	
	
	// make the output wave that will store the data
	dimension_sizes[0] = (long)x_size;
	dimension_sizes[1] = (long)y_size;
	dimension_sizes[2] = (long)n_images_to_load;
	dimension_sizes[3] = 0;
	
	// decide on the storage type to use
	// we use the storage type that matches that of the original frames
	switch (storage_type) {
		case STORAGE_TYPE_INT4:
		case STORAGE_TYPE_UINT4:
		case STORAGE_TYPE_INT8:
			waveType = NT_I8;
			break;
		case STORAGE_TYPE_UINT8:
			waveType = NT_I8 | NT_UNSIGNED;
			break;
		case STORAGE_TYPE_INT16:
			waveType = NT_I16;
			break;
		case STORAGE_TYPE_UINT16:
			waveType = NT_I16 | NT_UNSIGNED;
			break;
		case STORAGE_TYPE_INT32:
			waveType = NT_I32;
			break;
		case STORAGE_TYPE_UINT32:
			waveType = NT_I32 | NT_UNSIGNED;
			break;
		case STORAGE_TYPE_FP32:
			waveType = NT_FP32;
			break;
		case STORAGE_TYPE_FP64:
			waveType = NT_FP64;
			break;
		default:
			waveType = NT_FP64;
	}
	
	result = MDMakeWave(&output_wave, "M_CCDFrames", NULL, dimension_sizes, waveType, 1);
	
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
	
	result = SetOperationNumVar("V_numberOfImages", total_n_images);
	if (result != 0)
		throw result;
	result = SetOperationNumVar("V_xSize", x_size);
	if (result != 0)
		throw result;
	result = SetOperationNumVar("V_ySize", y_size);
	if (result != 0)
		throw result;
	result = SetOperationNumVar("V_firstImageLoaded", -1);
	if (result != 0)
		throw result;
	result = SetOperationNumVar("V_lastImageLoaded", -1);
	if (result != 0)
		throw result;
	
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
	
	if (this->nImages == 0) {	// if there are no images to load, do not do any processing
		return;
	}
	
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
	boost::shared_ptr<std::vector<LocalizedPosition> > fittedPositions;
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

	