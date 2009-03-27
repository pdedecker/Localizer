#include "PALM_analysis.h"


int load_partial_ccd_image(ImageLoader *image_loader, unsigned long n_start, unsigned long n_end) {
	unsigned long total_n_images = image_loader->get_total_number_of_images();
	unsigned long n_images_to_load;
	unsigned long x_size, y_size;
	waveHndl output_wave;
	long dimension_sizes[MAX_DIMENSIONS + 1];
	long indices[MAX_DIMENSIONS];
	
	int result;
	boost::shared_ptr<encap_gsl_matrix> current_image;
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
	
	x_size = image_loader->get_x_size();
	y_size = image_loader->get_y_size();
	
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
	for (unsigned long i = n_start; i <= n_end; i++) {
		current_image = image_loader->get_nth_image(i);
		
		indices[2] = i - n_start;
		
		// store the data in the output wave
		for (unsigned long k = 0; k < y_size; k++) {
			for (unsigned long j = 0; j < x_size; j++) {
				current_value = current_image->get(j, k);
				
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
	unsigned long total_n_images = image_loader->get_total_number_of_images();
	unsigned long x_size = image_loader->get_x_size();
	unsigned long y_size = image_loader->get_y_size();
	
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


int do_analyze_images_operation(boost::shared_ptr<ImageLoader> image_loader, const string output_wave_name, boost::shared_ptr<FitPositions> positions_fitter, 
								boost::shared_ptr<ParticleFinder> particle_finder, boost::shared_ptr<ThresholdImage_Preprocessor> preprocessor, 
								boost::shared_ptr<ThresholdImage> thresholder, boost::shared_ptr<ThresholdImage_Postprocessor> postprocessor) {
	
	unsigned long number_of_images;
	int status;
	
	number_of_images = image_loader->get_total_number_of_images();
	
	boost::shared_ptr<encap_gsl_matrix_uchar> thresholded_image;
	boost::shared_ptr<encap_gsl_matrix> positions;
	boost::shared_ptr<encap_gsl_matrix> current_image;
	boost::shared_ptr<encap_gsl_matrix> fitted_positions;
	
	unsigned long progress_indices[9];		// keeps track of the indices of the images that correspond to 10% done, 20% done, and so on
	int current_progress_index = 0;
	
	IgorOutputWriter output_writer(output_wave_name);
	
	for (int i = 0; i < 9; ++i) {
		progress_indices[i] = floor((double)(i + 1) / 10.0 * (double)number_of_images);
	}
	
	XOPNotice("Calculating");
	
	for (unsigned long i = 0; i < number_of_images; i++) {
		// check if the user wants to cancel the calculation
		status = CheckAbort(0);
		if (status == -1) {
			XOPNotice(" Abort requested by user\r");
			return 0;
		}
		
		current_image = image_loader->get_nth_image(i);
		
		thresholded_image = do_processing_and_thresholding(current_image, preprocessor, thresholder, postprocessor);
		
		// positions = find_positions_no_sort(current_image, thresholded_image, radius, min_distance_from_edge);
		positions = particle_finder->findPositions(current_image, thresholded_image);
		
		
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
			XOPNotice(".");
			current_progress_index++;
		}
		
	}
	
	// time to write the positions to an Igor wave
	XOPNotice(" Writing output... ");
	status = output_writer.write_positions_to_wave();
	if (status != 0) {
		return status;
	}
	
	XOPNotice(" Done!\r");
	
    return 0;
	
}

int do_analyze_images_operation_parallel(boost::shared_ptr<ImageLoader> image_loader, const string output_wave_name, boost::shared_ptr<FitPositions> positions_fitter, 
								boost::shared_ptr<ParticleFinder> particle_finder, boost::shared_ptr<ThresholdImage_Preprocessor> preprocessor, 
								boost::shared_ptr<ThresholdImage> thresholder, boost::shared_ptr<ThresholdImage_Postprocessor> postprocessor) {
	
	unsigned long number_of_images;
	int status;
	unsigned int numberOfProcessors = boost::thread::hardware_concurrency();
	unsigned int numberOfThreads;
	vector<boost::shared_ptr<boost::thread> > threads;
	boost::shared_ptr<boost::thread> singleThreadPtr;
	// vector<pthread_t> threads;
	// threads.resize(numberOfProcessors, pthread_t());
	// vector <boost::shared_ptr<encap_gsl_matrix> > localFittedPositions;
	// localFittedPositions.resize(numberOfProcessors, boost::shared_ptr<encap_gsl_matrix>());
	vector<unsigned long> startPositions;
	vector<unsigned long> endPositions;
	vector<threadStartData> threadData;
	unsigned long positionsPerThread, startPosition, endPosition;
	unsigned long nPositions, nLocalPositions;
	
	
	
	number_of_images = image_loader->get_total_number_of_images();
	
	boost::shared_ptr<encap_gsl_matrix_uchar> thresholded_image;
	boost::shared_ptr<encap_gsl_matrix> positions;
	boost::shared_ptr<encap_gsl_matrix> current_image;
	boost::shared_ptr<encap_gsl_matrix> fitted_positions;
	
	unsigned long progress_indices[9];		// keeps track of the indices of the images that correspond to 10% done, 20% done, and so on
	int current_progress_index = 0;
	
	IgorOutputWriter output_writer(output_wave_name);
	
	for (int i = 0; i < 9; ++i) {
		progress_indices[i] = floor((double)(i + 1) / 10.0 * (double)number_of_images);
	}
	
	XOPNotice("Calculating");
	
	for (unsigned long i = 0; i < number_of_images; i++) {
		// check if the user wants to cancel the calculation
		status = CheckAbort(0);
		if (status == -1) {
			XOPNotice( " Abort requested by user\r");
			return 0;
		}
		
		current_image = image_loader->get_nth_image(i);
		
		thresholded_image = do_processing_and_thresholding(current_image, preprocessor, thresholder, postprocessor);
		
		// positions = find_positions_no_sort(current_image, thresholded_image, radius, min_distance_from_edge);
		positions = particle_finder->findPositions(current_image, thresholded_image);
		
		if (positions.get() == NULL) {
			// if we get back NULL and no exception has been thrown then this means that we found no positions
			output_writer.append_new_positions(positions);	// we append a NULL pointer anyway because it keeps track of the image sequence
			
			if (i == progress_indices[current_progress_index]) {
				XOPNotice(".");
				current_progress_index++;
			}
			
			continue;
		}
		
		nPositions = positions->get_x_size();
		numberOfThreads = nPositions / 4;	// assume that each thread needs about 4 positions in order for the parallel processing to be useful
		if (numberOfThreads > numberOfProcessors) {
			numberOfThreads = numberOfProcessors;	// never spawn more threads than there are processors
		}
		if (numberOfThreads == 0) {
			numberOfThreads = 1;
		}
		
		positionsPerThread = nPositions / numberOfThreads;
		
		if (numberOfThreads == 1) {	// there is no real reason to bother with the increased overhead of the parallel approach, or the computer has only a single core
									// use the single-threaded approach (in this frame) and move on to the next one
			fitted_positions = positions_fitter->fit_positions(current_image, positions);
			
			// now we need to pass the fitted images to the output routine
			output_writer.append_new_positions(fitted_positions);
			
			if (i == progress_indices[current_progress_index]) {
				XOPNotice(".");
				current_progress_index++;
			}
			
			continue;
		}
		
		// allocate storage for the fitted results
		fitted_positions = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(nPositions, N_OUTPUT_PARAMS_PER_FITTED_POSITION));
		
		
		// divide the positions where we have to fit so that we can pass each thread some part of it
		startPositions.clear();
		endPositions.clear();
		threadData.clear();
		threads.clear();
		
		endPosition = (unsigned long)-1;
		for (unsigned long l = 0; l < numberOfThreads - 1; ++l) {
			startPosition = endPosition + 1;
			endPosition = startPosition + positionsPerThread - 1;
			nLocalPositions = endPosition - startPosition + 1;
			
			startPositions.push_back(startPosition);
			endPositions.push_back(endPosition);
		}
		
		// the last thread is special because we need to make sure that we include all positions
		startPosition = endPosition + 1;
		endPosition = nPositions - 1;
		nLocalPositions = endPosition - startPosition + 1;
		
		startPositions.push_back(startPosition);
		endPositions.push_back(endPosition);
		
		// now set up the data that will be passed to each thread
		for (unsigned long j = 0; j < numberOfThreads; ++j) {
			threadData.push_back(threadStartData(current_image, positions_fitter, positions, fitted_positions, startPositions.at(j), endPositions.at(j)));
		}
		
		// now start the threads
		for (unsigned long j = 0; j < numberOfThreads; ++j) {
			singleThreadPtr = boost::shared_ptr<boost::thread>(new boost::thread(&fitPositionsThreadStart, threadData.at(j)));
			threads.push_back(singleThreadPtr);
		}
		
		// wait for the treads to finish
		for (unsigned long j = 0; j < numberOfThreads; ++j) {
			threads.at(j)->join();
		}
			
		
		
		// now we need to pass the fitted positions to the output routine
		output_writer.append_new_positions(fitted_positions);
		
		if (i == progress_indices[current_progress_index]) {
			XOPNotice(".");
			current_progress_index++;
		}
		
	}
	
	// time to write the positions to an Igor wave
	XOPNotice(" Writing output... ");
	status = output_writer.write_positions_to_wave();
	if (status != 0) {
		return status;
	}
	
	XOPNotice(" Done!\r");
	
    return 0;
	
}


int do_analyze_images_operation_parallel2(boost::shared_ptr<ImageLoader> image_loader, const string output_wave_name, boost::shared_ptr<FitPositions> positions_fitter, 
										 boost::shared_ptr<ParticleFinder> particle_finder, boost::shared_ptr<ThresholdImage_Preprocessor> preprocessor, 
										 boost::shared_ptr<ThresholdImage> thresholder, boost::shared_ptr<ThresholdImage_Postprocessor> postprocessor) {
	
	unsigned long number_of_images;
	int status;
	unsigned int numberOfProcessors = boost::thread::hardware_concurrency();
	unsigned int numberOfThreads;
	vector<boost::shared_ptr<boost::thread> > threads;
	boost::shared_ptr<boost::thread> singleThreadPtr;
	vector<unsigned long> startPositions;
	vector<unsigned long> endPositions;
	vector<boost::shared_ptr<threadStartData2> > threadData;
	unsigned long nFramesRemaining;
	
	
	number_of_images = image_loader->get_total_number_of_images();
	boost::shared_ptr<encap_gsl_matrix> current_image;
	boost::shared_ptr<encap_gsl_matrix> fitted_positions;
	
	unsigned long progress_indices[9];		// keeps track of the indices of the images that correspond to 10% done, 20% done, and so on
	int current_progress_index = 0;
	
	IgorOutputWriter output_writer(output_wave_name);
	
	for (int i = 0; i < 9; ++i) {
		progress_indices[i] = floor((double)(i + 1) / 10.0 * (double)number_of_images);
	}
	
	XOPNotice("Calculating");
	
	// set up the vector containing the data for the threads
	for (unsigned int i = 0; i < numberOfProcessors; ++i) {
		threadData.push_back(boost::shared_ptr<threadStartData2> (new threadStartData2(thresholder, preprocessor, postprocessor, particle_finder, positions_fitter)));
	}
	
	for (unsigned long i = 0; i < number_of_images; ) {	// i is incremented when assigning the threads
		// check if the user wants to cancel the calculation
		status = CheckAbort(0);
		if (status == -1) {
			XOPNotice( " Abort requested by user\r");
			return 0;
		}
		
		nFramesRemaining = number_of_images - i;
		if (nFramesRemaining < numberOfProcessors) {
			numberOfThreads = nFramesRemaining;
		} else {
			numberOfThreads = numberOfProcessors;
		}
		
		// provide the starting image for the threads
		for (unsigned long j = 0; j < numberOfThreads; ++j) {
			current_image = image_loader->get_nth_image(i);
			++i;
			threadData.at(j)->image = current_image;
			
			if (i == progress_indices[current_progress_index]) {
				XOPNotice(".");
				DoUpdate();
				current_progress_index++;
			}
		}
		
		 // now start the threads
		threads.clear();
		for (unsigned long j = 0; j < numberOfThreads; ++j) {
			singleThreadPtr = boost::shared_ptr<boost::thread>(new boost::thread(&fitPositionsThreadStart2, threadData.at(j)));
			threads.push_back(singleThreadPtr);
		}
		
		// wait for the threads to finish
		for (unsigned long j = 0; j < numberOfThreads; ++j) {
			threads.at(j)->join();
		} 
		
		// output the fitted positions
		for (unsigned long j = 0; j < numberOfThreads; ++j) {
			output_writer.append_new_positions(threadData.at(j)->fittedPositions);
		}
	}
	
	// write the positions to an Igor wave
	XOPNotice(" Writing output... ");
	status = output_writer.write_positions_to_wave();
	if (status != 0) {
		return status;
	}
	
	XOPNotice(" Done!\r");
	
	return 0;
}


void fitPositionsThreadStart(threadStartData startData) {
	
	boost::shared_ptr<encap_gsl_matrix> positions = startData.positions;
	boost::shared_ptr<encap_gsl_matrix> image = startData.image;
	boost::shared_ptr<encap_gsl_matrix> outputPositions = startData.outputPositions;
	boost::shared_ptr<encap_gsl_matrix> localFittedPositions;
	
	boost::shared_ptr<FitPositions> positionsFitter = startData.positionsFitter;
	
	unsigned long start = startData.start;
	unsigned long end = startData.end;
	unsigned long nPositions = end - start + 1;
	
	// do the fit
	localFittedPositions = positionsFitter->fit_positions(image, positions, start, end);
	
	// copy the result into the output matrix
	for (unsigned long j = 0; j < nPositions; ++j) {
		for (unsigned long l = 0; l < N_OUTPUT_PARAMS_PER_FITTED_POSITION; ++l) {
			outputPositions->set(j + start, l, localFittedPositions->get(j, l));
		}
	}
	
}

void fitPositionsThreadStart2(boost::shared_ptr<threadStartData2> data) {
	
	boost::shared_ptr<encap_gsl_matrix_uchar> thresholded_image;
	boost::shared_ptr<encap_gsl_matrix> positions;
	boost::shared_ptr<encap_gsl_matrix> fitted_positions;
		
	thresholded_image = do_processing_and_thresholding(data->image, data->preprocessor, data->thresholder, data->postprocessor);
		
	positions = data->particleFinder->findPositions(data->image, thresholded_image);
	
	if (positions.get() == NULL) {
		// if we get back NULL and no exception has been thrown then this means that we found no positions
		data->fittedPositions = positions;
		
		return;
	}
	
	fitted_positions = data->positionsFitter->fit_positions(data->image, positions);
	data->fittedPositions = fitted_positions;
	
	return;
	
}
	


int do_analyze_images_operation_no_positions_finding(boost::shared_ptr<ImageLoader> image_loader, const string output_wave_name, waveHndl fitting_positions, 
													 boost::shared_ptr<FitPositions> positions_fitter) {
	
	// the format for passing the positions to the fitting routines is intensity, x, y
	
	unsigned long number_of_images;
	int status;
	long indices[MAX_DIMENSIONS + 1];
	long n_dimensions;
	unsigned long offset_in_positions_wave = 0;
	unsigned long total_n_positions;
	unsigned long n_positions_in_current_frame;
	double double_image_index;
	double current_value;
	double intensity;
	unsigned long ulong_image_index, counter;
	unsigned long ulong_current_value;
	unsigned long x_size, y_size;
	unsigned long current_x, current_y;
	
	unsigned long progress_indices[9];		// keeps track of the indices of the images that correspond to 10% done, 20% done, and so on
	int current_progress_index = 0;
	
	number_of_images = image_loader->get_total_number_of_images();
	x_size = image_loader->get_x_size();
	y_size = image_loader->get_y_size();
	
	boost::shared_ptr<encap_gsl_matrix> supplied_fitting_positions;
	boost::shared_ptr<encap_gsl_matrix> positions;
	boost::shared_ptr<encap_gsl_matrix> current_image;
	boost::shared_ptr<encap_gsl_matrix> fitted_positions;
	
	IgorOutputWriter output_writer(output_wave_name);
	
	for (int i = 0; i < 9; i++) {
		progress_indices[i] = floor((double)number_of_images / 100.0 * (double)(i + 1));
	}
	
	
	XOPNotice("Calculating");
	
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
	supplied_fitting_positions = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(total_n_positions, 3));
	
	// copy the data in the newly allocated array
	
	for (unsigned long i = 0; i < total_n_positions; i++) {
		indices[0] = i;
		indices[1] = 0;
		status = MDGetNumericWavePointValue(fitting_positions, indices, &double_image_index);
		
		if (status != 0) {
			throw status;
		}
		
		if (double_image_index < 0) {
			throw EXPECT_POS_NUM;
		}
		
		ulong_image_index = (unsigned long)(double_image_index + 0.5);
		
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
		
		ulong_current_value = (unsigned long)(current_value + 0.5);
		
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
		
		ulong_current_value = (unsigned long)(current_value + 0.5);
		
		if (ulong_current_value > (y_size - 1)) {
			throw INDEX_OUT_OF_RANGE;
		}
		
		supplied_fitting_positions->set(i, 2, current_value);
	}
	
	// a full copy of the wave data is now stored in supplied_fitting_positions
	offset_in_positions_wave = 0;
	
	for (unsigned long i = 0; i < number_of_images; i++) {
		
		n_positions_in_current_frame = 0;
		
		// we start by checking if we have some positions to fit in the current_frame
		counter = 0;
		
		while ((offset_in_positions_wave + counter < total_n_positions) && ((unsigned long)(supplied_fitting_positions->get(offset_in_positions_wave + counter, 0) + 0.5) == i)) {
			n_positions_in_current_frame++;
			counter++;
		}
		
		if (n_positions_in_current_frame == 0) {	// we found no positions
			// positions = NULL;
			positions = boost::shared_ptr<encap_gsl_matrix>();
		} else {	// we found some positions
			// we need to load the image
			try {
				current_image = image_loader->get_nth_image(i);
			}
			catch (OUT_OF_MEMORY) {
				throw NOMEM;
			}
			
			positions = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(n_positions_in_current_frame, 3));
			
			for (unsigned long j = offset_in_positions_wave; j < (offset_in_positions_wave + n_positions_in_current_frame); j++) {
				current_x = (unsigned long)(supplied_fitting_positions->get(j, 1) + 0.5);
				current_y = (unsigned long)(supplied_fitting_positions->get(j, 2) + 0.5);
				
				intensity = current_image->get(current_x, current_y);
				
				positions->set((j - offset_in_positions_wave), 0, intensity);
				positions->set((j - offset_in_positions_wave), 1, supplied_fitting_positions->get(j, 1));
				positions->set((j - offset_in_positions_wave), 2, supplied_fitting_positions->get(j, 2));
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
		
		//		print_fitted_positions(fitted_positions);
		
		
		// now we need to pass the fitted images to the output routine
		// this will also take care of freeing the matrix
		output_writer.append_new_positions(fitted_positions);
		
		if (i == progress_indices[current_progress_index]) {
			XOPNotice(".");
			current_progress_index++;
		}
		
		// if (((int)((double)i / (double)number_of_images * 100.0) % 10) == 0) {
		// 	XOPNotice(".");
		// }
		
	}
	
	// time to write the positions to an Igor wave
	XOPNotice(" Writing output... ");
	status = output_writer.write_positions_to_wave();
	if (status != 0) {
		return status;
	}
	
	XOPNotice(" Done!\r");
	
    return 0;
	
}


boost::shared_ptr<encap_gsl_matrix_uchar> do_processing_and_thresholding(boost::shared_ptr<encap_gsl_matrix> image, boost::shared_ptr<ThresholdImage_Preprocessor>preprocessor, 
																		 boost::shared_ptr<ThresholdImage> thresholder, boost::shared_ptr<ThresholdImage_Postprocessor> postprocessor) {
	// this function takes care of the thresholding and the associated pre- and postprocessing
	
	boost::shared_ptr<encap_gsl_matrix> preprocessed_image;
	boost::shared_ptr<encap_gsl_matrix_uchar> thresholded_image;
	boost::shared_ptr<encap_gsl_matrix_uchar> postprocessed_image;
	
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
	unsigned long n_images = image_loader->get_total_number_of_images();
	unsigned long x_size = image_loader->get_x_size();
	unsigned long y_size = image_loader->get_y_size();
	
	boost::shared_ptr<encap_gsl_matrix> current_image;
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
	
	for (unsigned long i = 0; i < n_images; i++) {
		summed_intensity = 0;
			current_image = image_loader->get_nth_image(i);
		
		// calculate the total sum of the image
		for (unsigned long k = startY; k <= endY; k++) {
			for (unsigned long j = startX; j <= endX; j++) {
				summed_intensity += current_image->get(j, k);
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


int construct_average_image(ImageLoader *image_loader, string output_wave_name, long startX, long startY, long endX, long endY) {
	unsigned long n_images = image_loader->get_total_number_of_images();
	unsigned long x_size = image_loader->get_x_size();
	unsigned long y_size = image_loader->get_y_size();
	
	long xRange, yRange;
	
	if (endX > x_size - 1)
		endX = x_size - 1;
	if (endY > y_size - 1)
		endY = y_size - 1;
	
	if ((startX < 0) || (startY < 0) || (startY < 0) || (startY < 0)) {	// we want the full frame
		startX = 0;
		startY = 0;
		endX = x_size - 1;
		endY = x_size - 1;
	}
	
	xRange = endX - startX + 1;
	yRange = endY - startY + 1;
	
	boost::shared_ptr<encap_gsl_matrix> current_image;
	boost::shared_ptr<encap_gsl_matrix> average_image(new encap_gsl_matrix(xRange, yRange));
	
	waveHndl output_wave;
	long dimension_sizes[MAX_DIMENSIONS + 1];
	dimension_sizes[0] = xRange;
	dimension_sizes[1] = yRange;
	dimension_sizes[2] = 0;
	long indices[MAX_DIMENSIONS];
	double current_value[2];
	double value;
	int result;
	
	average_image->set_all(0);
	
	
	for (unsigned long i = 0; i < n_images; i++) {
		current_image = image_loader->get_nth_image(i);
		
		// add the values of the newly loaded image to the average image
		for (long i = startX; i <= endX; ++i) {
			for (long j = startY; j <= endY; ++j) {
				value = average_image->get(i, j);
				value += current_image->get(i, j);
				average_image->set(i, j, value);
			}
		}
	}
	
	// divide by the number of images
	for (long i = startX; i <= endX; ++i) {
		for (long j = startY; j <= endY; ++j) {
			current_value[0] = average_image->get(i, j);
			current_value[0] /= (double)n_images;
			average_image->set(i, j, current_value[0]);
		}
	}
	
	// try to create the output wave
	result = MDMakeWave(&output_wave, output_wave_name.c_str(), NULL, dimension_sizes, NT_FP64, 1);
	if (result != 0) {
		throw result;
	}
	
	// write the output data to the wave
	for (long i = startX; i <= endX; ++i) {
		for (long j = startY; j <= endY; ++j) {
			current_value[0] = average_image->get(i, j);
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
	unsigned long n_images = image_loader->get_total_number_of_images();
	unsigned long x_size = image_loader->get_x_size();
	unsigned long y_size = image_loader->get_y_size();
	int result;
	waveHndl output_wave;
	long dimension_sizes[MAX_DIMENSIONS + 1];
	long indices[MAX_DIMENSIONS];
	
	double current_value[2], deviation, value;
	
	long xRange, yRange;
	
	if (endX > x_size - 1)
		endX = x_size - 1;
	if (endY > y_size - 1)
		endY = y_size - 1;
	
	if ((startX < 0) || (startY < 0) || (startY < 0) || (startY < 0)) {	// we want the full frame
		startX = 0;
		startY = 0;
		endX = x_size - 1;
		endY = x_size - 1;
	}
	
	xRange = endX - startX + 1;
	yRange = endY - startY + 1;
	
	boost::scoped_ptr<encap_gsl_matrix> stdDevImage(new encap_gsl_matrix(xRange, yRange));
	boost::scoped_ptr<encap_gsl_matrix> average_image(new encap_gsl_matrix(xRange, yRange));
	boost::shared_ptr<encap_gsl_matrix> current_image;
	
	average_image->set_all(0);
	stdDevImage->set_all(0);
	
	// construct an average image
	for (unsigned long i = 0; i < n_images; i++) {
		current_image = image_loader->get_nth_image(i);
		
		// add the values of the newly loaded image to the average image
		for (long j = startX; j <= endX; ++j) {
			for (long k = startY; k <= endY; ++k) {
				value = average_image->get(j, k);
				value += current_image->get(j, k);
				average_image->set(j, k, value);
			}
		}
	}
	
	// divide by the number of images
	for (long i = startX; i <= endX; ++i) {
		for (long j = startY; j <= endY; ++j) {
			value = average_image->get(i, j);
			value /= (double)n_images;
			average_image->set(i, j, value);
		}
	}
	
	// now loop over the images again, calculating the standard deviation of each pixel
	for (unsigned long i = 0; i < n_images; i++) {
		current_image = image_loader->get_nth_image(i);
		
		// add the deviation of the newly loaded image from the mean to the stddev jmage
		for (long j = startX; j <= endX; ++j) {
			for (long k = startY; k <= endY; ++k) {
				value = stdDevImage->get(j, k);
				deviation = (current_image->get(j, k) - average_image->get(j, k)) * (current_image->get(j, k) - average_image->get(j, k));
				value += deviation;
				stdDevImage->set(j, k, value);
			}
		}
	}
	
	// divide by the number of images to get the average deviation, and take the square root
	for (long i = startX; i <= endX; ++i) {
		for (long j = startY; j <= endY; ++j) {
			value = stdDevImage->get(i, j);
			value /= (double)n_images;
			value = sqrt(value);
			stdDevImage->set(i, j, value);
		}
	}
	
	// try to create the output wave
	dimension_sizes[0] = xRange;
	dimension_sizes[1] = yRange;
	dimension_sizes[2] = 0;
	result = MDMakeWave(&output_wave, output_wave_name.c_str(), NULL, dimension_sizes, NT_FP64, 1);
	if (result != 0) {
		throw result;
	}
	
	// write the output data to the wave
	for (long i = startX; i <= endX; ++i) {
		for (long j = startY; j <= endX; ++j) {
			current_value[0] = stdDevImage->get(i, j);
			indices[0] = i;
			indices[1] = j;
			result = MDSetNumericWavePointValue(output_wave, indices, current_value);
			if (result != 0) {
				throw result;
			}
		}
	}
}
	


gsl_histogram * make_histogram_from_matrix(boost::shared_ptr<encap_gsl_matrix> image, unsigned long number_of_bins) {
	unsigned long x_size, y_size;
	gsl_histogram *hist;
	double min = 1e100;
	double max = -1e100;
	double current_value;
	
	x_size = image->get_x_size();
	y_size = image->get_y_size();
	
	string error;
	error = "Unable to allocate a gsl_histogram in make_histogram_from_matrix()\r";
	
	hist = gsl_histogram_alloc(number_of_bins);
	if (hist == NULL) {
		throw OUT_OF_MEMORY(error);
		return NULL;
	}
	
	for (unsigned long j = 0; j < y_size; j++) {
		for (unsigned long i = 0; i < x_size; i++) {
			current_value = image->get(i, j);
			if (current_value < min)
				min = current_value;
			if (current_value > max)
				max = current_value;
		}
	}
	
	// adjust the histogram bins so that they range uniformly from min to max and set the values to zero
	gsl_histogram_set_ranges_uniform(hist, min, max);
	
	// now populate the histogram
	for (unsigned long j = 0; j < y_size; j++) {
		for (unsigned long i = 0; i < x_size; i++) {
			current_value = image->get(i, j);
			
			gsl_histogram_increment(hist, current_value);
			
		}
	}
	
	return hist;
}
	
	
/*
void print_localized_positions(gsl_matrix *positions) {
	unsigned long number_of_positions = (unsigned long)positions->size1;
	
	cout << "Number of positions found:\t" << number_of_positions << "\n";
	
	for (unsigned long i = 0; i < number_of_positions; i++) {
		cout << "(" << gsl_matrix_get(positions, i, 1) << "," << gsl_matrix_get(positions, i, 2) << ")\n";
	}
}

void print_fitted_positions(gsl_matrix *positions) {
	unsigned long number_of_positions = (unsigned long)positions->size1;
	
	cout << "Number of positions found:\t" << number_of_positions << "\n";
	
	for (unsigned long i = 0; i < number_of_positions; i++) {
		cout << "\nPosition " << i << ":\n";
		cout << "Amplitude:\t\t" << gsl_matrix_get(positions, i, 0) << "\n";
		cout << "Radius:\t\t\t" << gsl_matrix_get(positions, i, 1) << "\n";
		cout << "X position:\t\t" << gsl_matrix_get(positions, i, 2) << "\n";
		cout << "Y positions:\t\t" << gsl_matrix_get(positions, i, 3) << "\n";
		cout << "Offset:\t\t\t" << gsl_matrix_get(positions, i, 4) << "\n";
	}
}

void print_image_to_file(gsl_matrix *image, string file_path) {
	ofstream file(file_path.c_str(), ios::out | ios::trunc);
	
	unsigned long x_size = (unsigned long)image->size1;
	unsigned long y_size = (unsigned long)image->size2;
	
	for (unsigned long j = 0; j < y_size; j++) {
		for (unsigned long i = 0; i < x_size; i++) {
			file << gsl_matrix_get(image, i, j) << "\t";
		}
		file << "\n";
	}
}

void print_fit_state(unsigned long iteration, gsl_multifit_fdfsolver *s) {
	cout << "iteration: " << iteration << "\t" << gsl_vector_get(s->x, 0) << "\t" << gsl_vector_get(s->x, 1) << "\t" << gsl_vector_get(s->x, 2)
<< "\t" << gsl_vector_get(s->x, 3) << "\t" << gsl_vector_get(s->x, 4) << "\t\t|f(x)| = " << gsl_blas_dnrm2(s->f) << "\n";
}

void print_initial_parameters(gsl_vector *params) {
	cout << "\t\t\t" << gsl_vector_get(params, 0) << "\t" << gsl_vector_get(params, 1) << "\t" << gsl_vector_get(params, 2) << "\t" << gsl_vector_get(params, 3) << "\t"
<< gsl_vector_get(params, 4) << "\n";
}

void print_image(gsl_matrix * image) {
	unsigned long x_size = image->size1;
	unsigned long y_size = image->size2;
	
	for (unsigned long j = 0; j < y_size; j++) {
		for (unsigned long i = 0; i < x_size; i++) {
			if (i != x_size - 1) {
				cout <<  gsl_matrix_get(image, i, j) << "\t";
			}
			else {
				cout <<  gsl_matrix_get(image, i, j);
			}
		}
		cout << "\n";
	}
}
*/