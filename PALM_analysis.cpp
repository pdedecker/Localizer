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

int do_analyze_images_operation_no_positions_finding(boost::shared_ptr<ImageLoader> image_loader, const string output_wave_name, waveHndl fitting_positions, 
													 boost::shared_ptr<FitPositions> positions_fitter, int quiet) {
	
	// the format for passing the positions to the fitting routines is intensity, x, y
	
	size_t number_of_images;
	int status;
	long indices[MAX_DIMENSIONS + 1];
	long n_dimensions;
	size_t offset_in_positions_wave = 0;
	size_t total_n_positions;
	size_t n_positions_in_current_frame;
	double double_image_index;
	double current_value;
	double intensity;
	size_t ulong_image_index, counter;
	size_t ulong_current_value;
	size_t x_size, y_size;
	size_t current_x, current_y;
	
	size_t progress_indices[9];		// keeps track of the indices of the images that correspond to 10% done, 20% done, and so on
	int current_progress_index = 0;
	
	number_of_images = image_loader->get_total_number_of_images();
	x_size = image_loader->getXSize();
	y_size = image_loader->getYSize();
	
	boost::shared_ptr<PALMMatrix<double> > supplied_fitting_positions;
	boost::shared_ptr<PALMMatrix<double> > positions;
	boost::shared_ptr<PALMMatrix<double> > current_image;
	boost::shared_ptr<PALMMatrix<double> > fitted_positions;
	
	IgorOutputWriter output_writer(output_wave_name);
	
	for (int i = 0; i < 9; i++) {
		progress_indices[i] = floor((double)number_of_images / 100.0 * (double)(i + 1));
	}
	
	if (quiet == 0) {
		XOPNotice("Calculating");
	}
	
	// we will start by copying the positions into a gsl_matrix
	// even though this is quite unnecessary, it is easier to work with
	
	/*** BY CONVENTION THE POSITIONS HAVE TO BE FORMATTED IN THE FOLLOWING WAY: ***/
	// INDEX	X	Y
	
	indices[1] = 0;
	
	status = MDGetWaveDimensions(fitting_positions, &n_dimensions, indices);
	if (status != 0) {
		throw status;
	}
	
	total_n_positions = indices[0];
	supplied_fitting_positions = boost::shared_ptr<PALMMatrix<double> >(new PALMMatrix<double>(total_n_positions, 3));
	
	// copy the data in the newly allocated array
	
	for (size_t i = 0; i < total_n_positions; i++) {
		indices[0] = i;
		indices[1] = 0;
		status = MDGetNumericWavePointValue(fitting_positions, indices, &double_image_index);
		
		if (status != 0) {
			throw status;
		}
		
		if (double_image_index < 0) {
			throw EXPECT_POS_NUM;
		}
		
		ulong_image_index = (size_t)(double_image_index + 0.5);
		
		if (ulong_image_index > (number_of_images - 1)) {
			throw INDEX_OUT_OF_RANGE;
		}
		
		supplied_fitting_positions->set(i, 0, double_image_index);
		
		
		// get the x position
		indices[1] = 1;
		status = MDGetNumericWavePointValue(fitting_positions, indices, &current_value);
		
		if (status != 0) {
			throw status;
		}
		
		if (double_image_index < 0) {
			throw EXPECT_POS_NUM;
		}
		
		ulong_current_value = (size_t)(current_value + 0.5);
		
		if (ulong_current_value > (x_size - 1)) {
			throw INDEX_OUT_OF_RANGE;
		}
		
		supplied_fitting_positions->set(i, 1, current_value);
		
		// get the y position
		indices[1] = 2;
		status = MDGetNumericWavePointValue(fitting_positions, indices, &current_value);
		
		if (status != 0) {
			throw status;
		}
		
		if (double_image_index < 0) {
			throw EXPECT_POS_NUM;
		}
		
		ulong_current_value = (size_t)(current_value + 0.5);
		
		if (ulong_current_value > (y_size - 1)) {
			throw INDEX_OUT_OF_RANGE;
		}
		
		supplied_fitting_positions->set(i, 2, current_value);
	}
	
	// a full copy of the wave data is now stored in supplied_fitting_positions
	offset_in_positions_wave = 0;
	
	for (size_t i = 0; i < number_of_images; i++) {
		
		n_positions_in_current_frame = 0;
		
		// we start by checking if we have some positions to fit in the current_frame
		counter = 0;
		
		while ((offset_in_positions_wave + counter < total_n_positions) && ((size_t)(supplied_fitting_positions->get(offset_in_positions_wave + counter, 0) + 0.5) == i)) {
			n_positions_in_current_frame++;
			counter++;
		}
		
		if (n_positions_in_current_frame == 0) {	// we found no positions
			positions = boost::shared_ptr<PALMMatrix<double> >();
		} else {	// we found some positions
			// we need to load the image
			try {
				current_image = image_loader->get_nth_image(i);
			}
			catch (OUT_OF_MEMORY) {
				throw NOMEM;
			}
			
			positions = boost::shared_ptr<PALMMatrix<double> >(new PALMMatrix<double>(n_positions_in_current_frame, 3));
			
			for (size_t j = offset_in_positions_wave; j < (offset_in_positions_wave + n_positions_in_current_frame); j++) {
				current_x = (size_t)(supplied_fitting_positions->get(j, 1) + 0.5);
				current_y = (size_t)(supplied_fitting_positions->get(j, 2) + 0.5);
				
				intensity = current_image->get(current_x, current_y);
				
				(*positions)((j - offset_in_positions_wave), 0) = intensity;
				(*positions)((j - offset_in_positions_wave), 1) = (*supplied_fitting_positions)(j, 1);
				(*positions)((j - offset_in_positions_wave), 2) = (*supplied_fitting_positions)(j, 2);
			}
			
			offset_in_positions_wave += n_positions_in_current_frame;
			n_positions_in_current_frame = 0;
		}
		
		
		if (positions.get() == NULL) {
			// if we get back NULL and no exception has been thrown then this means that we found no positions
			output_writer.append_new_positions(positions);	// we append a NULL pointer anyway because it keeps track of the image sequence
			
			if (i == progress_indices[current_progress_index]) {
				XOPNotice(".");
				current_progress_index++;
			}
			
			continue;
		}
		
		
		fitted_positions = positions_fitter->fit_positions(current_image, positions);
		
		// now we need to pass the fitted images to the output routine
		output_writer.append_new_positions(fitted_positions);
		
		if (i == progress_indices[current_progress_index]) {
			if (quiet == 0) {
				XOPNotice(".");
			}
			current_progress_index++;
		}
		
	}
	
	// time to write the positions to an Igor wave
	if (quiet == 0) {
		XOPNotice(" Writing output... ");
	}
	status = output_writer.write_positions_to_wave();
	if (status != 0) {
		return status;
	}
	
	if (quiet == 0) {
		XOPNotice(" Done!\r");
	}
	
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
	


gsl_histogram * make_histogram_from_matrix(boost::shared_ptr<PALMMatrix<double> > image, size_t number_of_bins) {
	size_t x_size, y_size;
	gsl_histogram *hist;
	double min = 1e100;
	double max = -1e100;
	double current_value;
	
	x_size = image->getXSize();
	y_size = image->getYSize();
	
	string error;
	error = "Unable to allocate a gsl_histogram in make_histogram_from_matrix()\r";
	
	hist = gsl_histogram_alloc(number_of_bins);
	if (hist == NULL) {
		throw OUT_OF_MEMORY(error);
		return NULL;
	}
	
	for (size_t j = 0; j < y_size; j++) {
		for (size_t i = 0; i < x_size; i++) {
			current_value = (*image)(i, j);
			if (current_value < min)
				min = current_value;
			if (current_value > max)
				max = current_value;
		}
	}
	
	// adjust the histogram bins so that they range uniformly from min to max and set the values to zero
	gsl_histogram_set_ranges_uniform(hist, min, max + 1e-10);
	
	// now populate the histogram
	for (size_t j = 0; j < y_size; j++) {
		for (size_t i = 0; i < x_size; i++) {
			current_value = (*image)(i, j);
			
			gsl_histogram_increment(hist, current_value);
			
		}
	}
	
	return hist;
}
