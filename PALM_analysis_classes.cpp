/*
 *  PALM_analysis_classes.cpp
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 25/01/2008.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */


#include "PALM_analysis_classes.h"
#include "PALM_analysis.h"
#include "PALM_analysis_IgorXOP.h"







CCDImagesProcessorAverageSubtraction::CCDImagesProcessorAverageSubtraction(ImageLoader *i_loader, OutputWriter *o_writer) {
	image_loader = i_loader;
	output_writer = o_writer;
	
	total_number_of_images = image_loader->get_total_number_of_images();
	x_size = image_loader->getXSize();
	y_size = image_loader->getYSize();
	
	n_frames_averaging = 0;
}


void CCDImagesProcessorAverageSubtraction::set_n_parameter(double n) {
	// check is we are indeed averaging over an odd number
	// if not then we throw an error
	// the exception is that a value of '0' means that we have to average over the entire sequence
	
	// by convention the calling function checks that n is positive
	n_frames_averaging = (size_t)(n + 0.5);
	
	if (((n_frames_averaging % 2) != 1) && (n_frames_averaging != 0)) {
		throw NUMBER_OF_AVERAGING_FRAMES_SHOULD_BE_ODD();
	}
	
	if (n_frames_averaging > total_number_of_images) {
		throw NUMBER_OF_AVERAGING_LARGER_THAN_N_FRAMES_IN_FILE();
	}
}


int CCDImagesProcessorAverageSubtraction::convert_images() {
	if (n_frames_averaging == 0) {	// we want to average over the entire trace
		subtract_average_of_entire_trace();
	} else {
		subtract_partial_average();
	}
	
	return 0;
}

void CCDImagesProcessorAverageSubtraction::subtract_average_of_entire_trace() {
	size_t n;
	boost::shared_ptr<PALMMatrix<double> > average_image;
	boost::shared_ptr<PALMMatrix<double> > loaded_image;
	boost::shared_ptr<PALMMatrix<double> > subtracted_image;
	double current_double;
	double value;
	
	// we pass through the images two times:
	// the first pass calculates the average,
	// the second pass subtracts it from the image
	// fortunately out intermediate format uses doubles to store the data!
	
	average_image = boost::shared_ptr<PALMMatrix<double> >(new PALMMatrix<double>(x_size, y_size));
	
	average_image->set_all(0);	// zero the matrix
	
	for (n = 0; n < total_number_of_images; n++) {
		loaded_image = image_loader->get_nth_image(n);
		
		// average_image->add(*loaded_image);
		for (size_t l = 0; l < y_size; ++l) {
			for (size_t k = 0; k < x_size; ++k) {
				value = average_image->get(k, l);
				value += loaded_image->get(k, l);
				average_image->set(k, l, value);
			}
		}
	}
	
	// now divide each point so that we get the average
	// gsl_matrix_scale(average_image, (1.0 / (double)n_frames_averaging));
	
	// the approach using the gsl functions seems off so we use a different one instead
	
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {
			current_double = average_image->get(i, j);
			current_double /= (double)total_number_of_images;
			average_image->set(i, j, current_double);
		}
	}
	
	// now subtract the average for each frame
	for (n = 0; n < total_number_of_images; n++) {
		loaded_image = image_loader->get_nth_image(n);
		
		// loaded_image->sub(*average_image);
		for (size_t k = 0; k < x_size; ++k) {
			for (size_t l = 0; l < y_size; ++l) {
				value = loaded_image->get(k, l);
				value -= average_image->get(k, l);
				loaded_image->set(k, l, value);
			}
		}
		
		subtracted_image = loaded_image;
		
		output_writer->write_image(subtracted_image);	// the output writer will take care of freeing the memory
	}
}



void CCDImagesProcessorAverageSubtraction::subtract_partial_average() {
	boost::shared_ptr<PALMMatrix<double> > current_image;
	boost::shared_ptr<PALMMatrix<double> > average_image;
	boost::shared_ptr<PALMMatrix<double> > subtracted_image;
	// gsl_matrix *averaging_buffer[n_frames_averaging];	// strictly speaking this isn't legal C++ code because declarations on the stack should be
															// a const size
															// g++ accepts this, but the Microsoft compiler throws an error on this
															// so we will convert it to an assignment on the heap instead
	vector<boost::shared_ptr<PALMMatrix<double> > > averaging_buffer;
	averaging_buffer.resize(n_frames_averaging, boost::shared_ptr<PALMMatrix<double> > ());
	
	long average_starting_index, average_ending_index;
	size_t cache_loading_offset = 0;
	
	double current_double;
	double value;
	
	average_image = boost::shared_ptr<PALMMatrix<double> >(new PALMMatrix<double>(x_size, y_size));
	
	/*// the buffer for the averaging will be allocated by the get_nth_image() routine
	for (size_t i = 0; i < n_frames_averaging; i++) {
		averaging_buffer[i] = NULL;
	}*/
	
	// we normally try to subtract the average of the frames surrounding the frame that we are interested in
	// this is why we demand that n_frames_averaging is odd
	
	// however, if we have a frame at the beginning or end of the image stack then we cannot do this
	// instead we construct an average of n_frames_averaging as close as possible to the frame we are interested in
	
	// loop over all the images
	for (size_t n = 0; n < total_number_of_images; n++) {
		average_starting_index = n - floor((double)total_number_of_images / 2.0);
		average_ending_index = n + floor((double)total_number_of_images / 2.0);
		
		// SPECIAL CASE: check if we are too close to the START of the image stack to calculate a normal average
		
		if (average_starting_index <= 0) {
			// all the points that fulfill this condition calculate their average using the same set of frames
			// so we only load these frames once, for the first image
			if (n == 0) {
				for (size_t i = 0; i < n_frames_averaging; i++) {
					
					averaging_buffer.at(i) = image_loader->get_nth_image(i);
					
				}
				
				// now calculate the average once
				average_image->set_all(0);
				for (size_t i = 0; i < n_frames_averaging; i++) {
					for (size_t k = 0; k < x_size; ++k) {
						for (size_t l = 0; l < y_size; ++l) {
							value = average_image->get(k, l);
							value += averaging_buffer.at(i)->get(k, l);
							average_image->set(k, l, value);
						}
					}
				}
				// gsl_matrix_scale(average_image, (1.0 / (double)n_frames_averaging));
				// the calculation using the gsl seems to be off so we use a direct one instead
				
				for (size_t i = 0; i < x_size; i++) {
					for (size_t j = 0; j < y_size; j++) {
						current_double = average_image->get(i, j);
						current_double /= n_frames_averaging;
						average_image->set(i, j, current_double);
					}
				}
				
			}
			// the frame that we are interested in (index n) can now also be found at
			// averaging_buffer[n]
			current_image = averaging_buffer.at(n);
			
			// now allocate the subtracted image
			subtracted_image = boost::shared_ptr<PALMMatrix<double> >(new PALMMatrix<double>(x_size, y_size));
			
			// subtract the average from the current frame
			// subtracted_image->copy(*current_image);
			for (size_t j = 0; j < y_size; ++j) {
				for (size_t i = 0; i < x_size; ++i) {
					subtracted_image->set(i, j, current_image->get(i, j));
				}
			}
			// subtracted_image->sub(*average_image);
			for (size_t k = 0; k < x_size; ++k) {
				for (size_t l = 0; l < y_size; ++l) {
					value = subtracted_image->get(k, l);
					value -= average_image->get(k, l);
					subtracted_image->set(k, l, value);
				}
			}
			
			// store the subtracted image
			output_writer->write_image(subtracted_image);	// output_writer will take care of freeing the memory
			
			continue;	// don't go through any of the other stuff in the main for loop
						// while we are too close to the edge of the image stack
		}
		
		
		// SPECIAL CASE: check if we are too close to the END of the image stack to calculate a normal average
		if ((size_t)average_ending_index >= total_number_of_images) {
			
			// we don't need to update the image buffer anymore, all the images we need to calculate the average are stored in memory
			// also, we don't need to calculate this average as it will have been calculated by the last image
			// that could calculate a normal average
			
			current_image = averaging_buffer[(size_t)floor((double)(n_frames_averaging) / 2.0) + cache_loading_offset];
			
			// now allocate the subtracted image
			subtracted_image = boost::shared_ptr<PALMMatrix<double> >(new PALMMatrix<double>(x_size, y_size));
			
			// subtracted_image->copy(*current_image);
			for (size_t j = 0; j < y_size; ++j) {
				for (size_t i = 0; i < x_size; ++i) {
					subtracted_image->set(i, j, current_image->get(i, j));
				}
			}
			
			
			// subtracted_image->sub(*average_image);
			for (size_t k = 0; k < x_size; ++k) {
				for (size_t l = 0; l < y_size; ++l) {
					value = subtracted_image->get(k, l);
					value -= average_image->get(k, l);
					subtracted_image->set(k, l, value);
				}
			}
			
			output_writer->write_image(subtracted_image);
			
			cache_loading_offset++;
			continue;
		}
				
		
		// NORMAL CASE: we are somewhere in the middle of the image stack
		
		averaging_buffer.erase(averaging_buffer.begin());
		// we shift the images to the right by one position
		for (size_t i = 0; i < (n_frames_averaging - 1); i++) {
			
			averaging_buffer.at(i) = averaging_buffer.at(i + 1);
		}
		
		averaging_buffer.at(n_frames_averaging - 1) = image_loader->get_nth_image(n + floor((double)(n_frames_averaging) / 2.0));
		
		// if we are here then we can calculate a new average
		average_image->set_all(0);
		
		for (size_t i = 0; i < n_frames_averaging; i++) {
			for (size_t i = 0; i < x_size; ++i) {
				for (size_t j = 0; j < y_size; ++j) {
					value = average_image->get(i, j);
					value += current_image->get(i, j);
					average_image->set(i, j, value);
				}
			}
		}
		
		// gsl_matrix_scale(average_image, (1.0 / (double)n_frames_averaging));
		// the calculation using the gsl seems to be off so we use a direct one instead
		
		for (size_t i = 0; i < x_size; i++) {
			for (size_t j = 0; j < y_size; j++) {
				current_double = average_image->get(i, j);
				current_double /= n_frames_averaging;
				average_image->set(i, j, current_double);
			}
		}
		
		// the image that we are interested in is at the center of the imaging buffer
		current_image = averaging_buffer.at((size_t)floor((double)(n_frames_averaging) / 2.0));
		
		// now allocate the subtracted image
		subtracted_image = boost::shared_ptr<PALMMatrix<double> >(new PALMMatrix<double>(x_size, y_size));
		
		// subtracted_image->copy(*current_image);
		for (size_t j = 0; j < y_size; ++j) {
			for (size_t i = 0; i < x_size; ++i) {
				subtracted_image->set(i, j, current_image->get(i, j));
			}
		}
		
		// subtracted_image->sub(*average_image);
		for (size_t k = 0; k < x_size; ++k) {
			for (size_t l = 0; l < y_size; ++l) {
				value = subtracted_image->get(k, l);
				value -= average_image->get(k, l);
				subtracted_image->set(k, l, value);
			}
		}
		
		output_writer->write_image(subtracted_image);
	}
	
	// phew! we're done
	// now we just have to make sure that we don't leave a mess behind when we close
	
}
	
	
CCDImagesProcessorDifferenceImage::CCDImagesProcessorDifferenceImage(ImageLoader *i_loader, OutputWriter *o_writer) {
	image_loader = i_loader;
	output_writer = o_writer;
	
	total_number_of_images = image_loader->get_total_number_of_images();
	x_size = image_loader->getXSize();
	y_size = image_loader->getYSize();
}

int CCDImagesProcessorDifferenceImage::convert_images() {
	
	boost::shared_ptr<PALMMatrix<double> > current_image;
	boost::shared_ptr<PALMMatrix<double> > next_image;
	
	// we start by loading the first image in next_image
	// this is required so the loop that follows can start properly
	
	next_image = image_loader->get_nth_image(0);
	
	
	for (size_t n = 0; n < (total_number_of_images - 1); n++) {
		
		// the previous image for this run of the loop is the image that was previously in current_image
		// so we have to shift it down
		current_image = next_image;
		
		next_image = image_loader->get_nth_image(n + 1);
		
		// now do the actual subtraction
		*current_image = *current_image - *next_image;
		
		// current_image now contains the subtracted image, we should write it to disk
		output_writer->write_image(current_image);
		
		// the output_writer also takes care of freeing current_image
	}
	
	// before we exit the function we need to free next_image
	// gsl_matrix_free(next_image);
	
	return 0;
}
		

CCDImagesProcessorConvertToSimpleFileFormat::CCDImagesProcessorConvertToSimpleFileFormat(ImageLoader *i_loader, OutputWriter *o_writer) {
	image_loader = i_loader;
	output_writer = o_writer;
	
	total_number_of_images = image_loader->get_total_number_of_images();
	x_size = image_loader->getXSize();
	y_size = image_loader->getYSize();
}

int CCDImagesProcessorConvertToSimpleFileFormat::convert_images() {
	
	boost::shared_ptr<PALMMatrix<double> > current_image;
	
	for (size_t n = 0; n < total_number_of_images; ++n) {
		current_image = image_loader->get_nth_image(n);
		output_writer->write_image(current_image);
	}
	return 0;
}


boost::shared_ptr<PALMMatrix<double> > ParticleFinder_radius::findPositions(boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<PALMMatrix <unsigned char> > threshold_image) {
	vector<position> positions;
	position position;
	// we store the pixels above the treshold as a vector containing x,y,intensity
	size_t x_size = image->getXSize(), y_size = image->getYSize();
	size_t number_of_positions;
	double current_intensity, previous_intensity, current_x, current_y;
	double distance_squared;
	double radius_squared = (double)(radius * radius);
	int skip;
	double x, y;
	boost::shared_ptr<PALMMatrix<double> > output_positions;
	
	// we run over all the points in the image to see if they are above the treshold
	for (size_t j = minDistanceFromEdge; j < y_size - minDistanceFromEdge; j++) {
		for (size_t i = minDistanceFromEdge; i < x_size - minDistanceFromEdge; i++) {
			
			if (threshold_image->get(i, j) < 128)
				continue;	// we don't care about this point, it's not included in the thresholded image
			
			current_intensity = image->get(i, j);
			
			current_x = (double)i;
			current_y = (double)j;
			
			// is this point too close to the edge of the image?
			if ((current_x < minDistanceFromEdge) || (current_x > (x_size - minDistanceFromEdge)))
				continue;
			if ((current_y < minDistanceFromEdge) || (current_y > (y_size - minDistanceFromEdge)))
				continue;
			
			// if we are still here then we need to take a closer look at this point
			// check if the current point overlaps with a previous point
			skip = 0;
			number_of_positions = positions.size();
			for (size_t k = 0; k < number_of_positions; k++) {
				x = positions[k].get_x();
				y = positions[k].get_y();
				distance_squared = (current_x - x) * (current_x - x) + (current_y - y) * (current_y - y);
				
				if (distance_squared < radius_squared) {
					// we need to skip one of the two pixels that we are comparing
					// we will keep the pixel with the largest intensity
					previous_intensity = positions[k].get_intensity();
					if (current_intensity > previous_intensity) {
						positions[k].set_intensity(current_intensity);
						positions[k].set_x(current_x);
						positions[k].set_y(current_y);
					}
					skip = 1;
					break;
				}
			}
			
			if (skip == 0) {	// we should store this point
				position.set_intensity(current_intensity);
				position.set_x(current_x);
				position.set_y(current_y);
				positions.push_back(position);
			}
		}
	}
	
	// now we need to store the data in the standard matrix format
	// this means that the columns are oriented as intensity, x, y
	
	number_of_positions = positions.size();
	
	if (number_of_positions == 0) {	// no positions were found
		return output_positions;	// is equal to NULL due to its initialization as a shared_ptr
	}
	
	output_positions = boost::shared_ptr<PALMMatrix<double> >(new PALMMatrix<double>(number_of_positions, 3));
	
	for (size_t i = 0; i < number_of_positions; i++) {
		output_positions->set(i, 0, positions[i].get_intensity());
		output_positions->set(i, 1, positions[i].get_x());
		output_positions->set(i, 2, positions[i].get_y());
	}
	
	return output_positions;
}


boost::shared_ptr<PALMMatrix<double> > ParticleFinder_adjacent4::findPositions(boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<PALMMatrix <unsigned char> > threshold_image) {
	
	boost::shared_ptr<PALMMatrix<double> > output_positions;
	list<position> positionsInCurrentParticleList;
	vector<position> positionsInCurrentParticle;
	position currentPosition;
	vector<position> particles;
	boost::shared_ptr<PALMMatrix<long> > mapped_image;	// keeps track of which pixels have already been mapped to a particle
	size_t x_size = image->getXSize();
	size_t y_size = image->getYSize();
	size_t x, y;
	long particleIndex = 0;
	double average_x, average_y;
	double maxIntensity;
	
	mapped_image = boost::shared_ptr<PALMMatrix<long> >(new PALMMatrix<long>(x_size, y_size));
	
	mapped_image->set_all(-1);
	
	for (size_t j = minDistanceFromEdge; j < y_size - minDistanceFromEdge; ++j) {
		for (size_t i = minDistanceFromEdge; i < x_size - minDistanceFromEdge; ++i) {	// loop over the entire image
			
			if (threshold_image->get(i, j) < 128) {
				continue;
			}
			
			if (mapped_image->get(i, j) != -1) {	// this point is already assigned to another particle
				continue;
			}
			
			// if we are still here then we found an active pixel that is not yet assigned to a new particle
			// we create a new particle at this position
			positionsInCurrentParticle.clear();
			positionsInCurrentParticleList.clear();
			
			mapped_image->set(i, j, particleIndex);
			
			// store this position
			currentPosition.set_x((double)i);
			currentPosition.set_y((double)j);
			currentPosition.set_intensity(image->get(i, j));
			
			positionsInCurrentParticleList.push_back(currentPosition);
			
			// growParticle(currentPosition, i, j, threshold_image, mapped_image);
			
			
			while (positionsInCurrentParticleList.size() > 0) {
				currentPosition = positionsInCurrentParticleList.front();
				x = (size_t)(currentPosition.get_x() + 0.5);
				y = (size_t)(currentPosition.get_y() + 0.5);
				
				growParticle(currentPosition, positionsInCurrentParticleList, image, threshold_image, mapped_image);
				// growParticle will update the list with new positions
				// this position has been checked, so we don't need to include it in future searches
				positionsInCurrentParticleList.pop_front();
				positionsInCurrentParticle.push_back(currentPosition);
			}
			
			// we have found all positions in this particle
			++particleIndex;
			
			// store the output positions
			// first calculate an average in x and y
			// and get an estimate for the intensity of the particle
			maxIntensity = 0;
			average_x = 0;
			average_y = 0;
			for (size_t k = 0; k < positionsInCurrentParticle.size(); ++k) {
				average_x += positionsInCurrentParticle[k].get_x();
				average_y += positionsInCurrentParticle[k].get_y();
				if (positionsInCurrentParticle[k].get_intensity() > maxIntensity) {
					maxIntensity = positionsInCurrentParticle[k].get_intensity();
				}
			}
			average_x /= positionsInCurrentParticle.size();
			average_y /= positionsInCurrentParticle.size();
			currentPosition.set_intensity(maxIntensity);
			currentPosition.set_x(average_x);
			currentPosition.set_y(average_y);
			
			// we store the particle only if it is not too close to the edge of the frame
			if ((average_x < minDistanceFromEdge) || (average_x > ((double)x_size - minDistanceFromEdge)))
				continue;
			if ((average_y < minDistanceFromEdge) || (average_y > ((double)y_size - minDistanceFromEdge)))
				continue;
			
			particles.push_back(currentPosition);
		}
	}
	
	// if we have found some particles then return them, else return a null pointer
	if (particles.size() == 0) {
		return output_positions;
	}
	
	output_positions = boost::shared_ptr<PALMMatrix<double> >(new PALMMatrix<double>(particles.size(), 3));
	
	// now copy the output to a gsl matrix
	for (size_t k = 0; k < particles.size(); ++k) {
		output_positions->set(k, 0, particles[k].get_intensity());
		output_positions->set(k, 1, particles[k].get_x());
		output_positions->set(k, 2, particles[k].get_y());
	}
	
	return output_positions;
	
}

void ParticleFinder_adjacent4::growParticle(position centerPosition, list<position> &positionsInCurrentParticle, boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<PALMMatrix <unsigned char> > threshold_image, boost::shared_ptr<PALMMatrix<long> > mapped_image) {
	// the pixel at position (x,y) belongs to a particle
	// do the surrounding pixels belong to the same particle?
	
	// the function checks which of the pixels surrounding pos are active
	// if one or more of these is active, then it checks if they are already assigned to the particle by checking mapped_image
	// if they are not known then they are added to to the list with positions of the current particle
	// and also added to mapped_image
	
	size_t x_size = image->getXSize();
	size_t y_size = image->getYSize();
	position currentPosition;
	
	size_t x = (size_t)(centerPosition.get_x() + 0.5);
	size_t y = (size_t)(centerPosition.get_y() + 0.5);
	
	long particleIndex = mapped_image->get(x, y);
	
	if ((x < minDistanceFromEdge) || (x > x_size - minDistanceFromEdge - 1))
		return;
	if ((y < minDistanceFromEdge) || (y > y_size - minDistanceFromEdge - 1))
		return;
	
	// is the pixel to the left of the current one active?
	if (x > 0) {
		if (threshold_image->get(x - 1, y) > 128) {
			// it's active
			// did we already include this pixel?
			if (mapped_image->get(x - 1, y) == -1) {
				mapped_image->set(x - 1, y, particleIndex);
				// add the point to the vector
				currentPosition.set_x((double)x - 1);
				currentPosition.set_y((double)y);
				currentPosition.set_intensity(image->get(x - 1, y));
				positionsInCurrentParticle.push_back(currentPosition);
				// now call this function itself (recursion) on the new pixel
				// growParticle(positions, x - 1, y, threshold_image, mapped_image);
			}
		}
	}
	// is the pixel to the right of the current one active?
	if (x < x_size - 1) {
		if (threshold_image->get(x + 1, y) > 128) {
			// it's active
			// did we already include this pixel?
			if (mapped_image->get(x + 1, y) == -1) {
				mapped_image->set(x + 1, y, particleIndex);
				// add the point to the vector
				currentPosition.set_x((double)x + 1);
				currentPosition.set_y((double)y);
				currentPosition.set_intensity(image->get(x + 1, y));
				positionsInCurrentParticle.push_back(currentPosition);
				// now call this function itself (recursion) on the new pixel
				// growParticle(positions, x + 1, y, threshold_image, mapped_image);
			}
		}
	}
	// is the pixel above the current one active?
	if (y > 0) {
		if (threshold_image->get(x, y - 1) > 128) {
			// it's active
			// did we already include this pixel?
			if (mapped_image->get(x, y - 1) == -1) {
				mapped_image->set(x, y - 1, particleIndex);
				// add the point to the vector
				currentPosition.set_x((double)x);
				currentPosition.set_y((double)y - 1);
				currentPosition.set_intensity(image->get(x, y - 1));
				positionsInCurrentParticle.push_back(currentPosition);
				// now call this function itself (recursion) on the new pixel
				// growParticle(positions, x, y - 1, threshold_image, mapped_image);
			}
		}
	}
	// is the pixel below the current one active?
	if (y < y_size - 1) {
		if (threshold_image->get(x, y + 1) > 128) {
			// it's active
			// did we already include this pixel?
			if (mapped_image->get(x, y + 1) == -1) {
				mapped_image->set(x, y + 1, particleIndex);
				// add the point to the vector
				currentPosition.set_x((double)x);
				currentPosition.set_y((double)y + 1);
				currentPosition.set_intensity(image->get(x, y + 1));
				positionsInCurrentParticle.push_back(currentPosition);
				// now call this function itself (recursion) on the new pixel
				// growParticle(positions, x, y + 1, threshold_image, mapped_image);
			}
		}
	}
}


boost::shared_ptr<PALMMatrix<double> > ParticleFinder_adjacent8::findPositions(boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<PALMMatrix <unsigned char> > threshold_image) {
	
	boost::shared_ptr<PALMMatrix<double> > output_positions;
	list<position> positionsInCurrentParticleList;
	vector<position> positionsInCurrentParticle;
	position currentPosition;
	vector<position> particles;
	boost::shared_ptr<PALMMatrix<long> > mapped_image;	// keeps track of which pixels have already been mapped to a particle
	size_t x_size = image->getXSize();
	size_t y_size = image->getYSize();
	size_t x, y;
	long particleIndex = 0;
	double average_x, average_y;
	double maxIntensity;
	
	mapped_image = boost::shared_ptr<PALMMatrix<long> >(new PALMMatrix<long>(x_size, y_size));
	
	mapped_image->set_all(-1);
	
	for (size_t j = minDistanceFromEdge; j < y_size - minDistanceFromEdge; ++j) {
		for (size_t i = minDistanceFromEdge; i < x_size - minDistanceFromEdge; ++i) {	// loop over the entire image
			
			if (threshold_image->get(i, j) < 128) {
				continue;
			}
			
			if (mapped_image->get(i, j) != -1) {	// this point is already assigned to another particle
				continue;
			}
			
			// if we are still here then we found an active pixel that is not yet assigned to a new particle
			// we create a new particle at this position
			positionsInCurrentParticle.clear();
			positionsInCurrentParticleList.clear();
			
			mapped_image->set(i, j, particleIndex);
			
			// store this position
			currentPosition.set_x((double)i);
			currentPosition.set_y((double)j);
			currentPosition.set_intensity(image->get(i, j));
			
			positionsInCurrentParticleList.push_back(currentPosition);
			
			// growParticle(currentPosition, i, j, threshold_image, mapped_image);
			
			
			while (positionsInCurrentParticleList.size() > 0) {
				currentPosition = positionsInCurrentParticleList.front();
				x = (size_t)(currentPosition.get_x() + 0.5);
				y = (size_t)(currentPosition.get_y() + 0.5);
				
				growParticle(currentPosition, positionsInCurrentParticleList, image, threshold_image, mapped_image);
				// growParticle will update the list with new positions
				// this position has been checked, so we don't need to include it in future searches
				positionsInCurrentParticleList.pop_front();
				positionsInCurrentParticle.push_back(currentPosition);
			}
			
			
			// we're done with this particle, time for the next one
			++particleIndex;
			
			// store the output positions
			// first calculate an average in x and y
			// and get an estimate for the intensity of the particle
			maxIntensity = 0;
			average_x = 0;
			average_y = 0;
			for (size_t k = 0; k < positionsInCurrentParticle.size(); ++k) {
				average_x += positionsInCurrentParticle[k].get_x();
				average_y += positionsInCurrentParticle[k].get_y();
				if (positionsInCurrentParticle[k].get_intensity() > maxIntensity) {
					maxIntensity = positionsInCurrentParticle[k].get_intensity();
				}
			}
			average_x /= positionsInCurrentParticle.size();
			average_y /= positionsInCurrentParticle.size();
			currentPosition.set_intensity(maxIntensity);
			currentPosition.set_x(average_x);
			currentPosition.set_y(average_y);
			
			// we store the particle only if it is not too close to the edge of the frame
			if ((average_x < minDistanceFromEdge) || (average_x > ((double)x_size - minDistanceFromEdge)))
				continue;
			if ((average_y < minDistanceFromEdge) || (average_y > ((double)y_size - minDistanceFromEdge)))
				continue;
			
			particles.push_back(currentPosition);
		}
	}
	
	// if we have found some particles then return them, else return a null pointer
	if (particles.size() == 0) {
		return output_positions;
	}
	
	output_positions = boost::shared_ptr<PALMMatrix<double> >(new PALMMatrix<double>(particles.size(), 3));
	
	// now copy the output to a gsl matrix
	for (size_t k = 0; k < particles.size(); ++k) {
		output_positions->set(k, 0, particles[k].get_intensity());
		output_positions->set(k, 1, particles[k].get_x());
		output_positions->set(k, 2, particles[k].get_y());
	}
	
	return output_positions;
	
}

void ParticleFinder_adjacent8::growParticle(position centerPosition, list<position> &positionsInCurrentParticle, boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<PALMMatrix <unsigned char> > threshold_image, boost::shared_ptr<PALMMatrix<long> > mapped_image) {
	// the pixel at position (x,y) belongs to a particle
	// do the surrounding pixels belong to the same particle?
	
	// the function checks which of the pixels surrounding pos are active
	// if one or more of these is active, then it checks if they are already assigned to the particle by checking mapped_image
	// if they are not known then they are added to to the list with positions of the current particle
	// and also added to mapped_image
	
	size_t x_size = image->getXSize();
	size_t y_size = image->getYSize();
	position currentPosition;
	
	size_t x = (size_t)(centerPosition.get_x() + 0.5);
	size_t y = (size_t)(centerPosition.get_y() + 0.5);
	
	long particleIndex = mapped_image->get(x, y);
	
	if ((x < minDistanceFromEdge) || (x > x_size - minDistanceFromEdge - 1))
		return;
	if ((y < minDistanceFromEdge) || (y > y_size - minDistanceFromEdge - 1))
		return;
	
	assert((x > 0) && (y > 0));
	
	for (size_t j = y - 1; j <= y + 1; ++j) {
		for (size_t i = x - 1; i <= x + 1; ++i) {
			
			if ((i == 0) && (j == 0))
				continue;
			
			if (threshold_image->get(i, j) < 128)
				continue;
			
			if (mapped_image->get(i, j) != -1)
				continue;
			
			// add the current position
			mapped_image->set(i, j, particleIndex);
			currentPosition.set_x((double)i);
			currentPosition.set_y((double)j);
			currentPosition.set_intensity(image->get(i, j));
			positionsInCurrentParticle.push_back(currentPosition);
		}
	}
}

boost::shared_ptr<PALMMatrix <unsigned char> > ThresholdImage_Direct::do_thresholding(boost::shared_ptr<PALMMatrix<double> > image) {
	
	size_t x_size, y_size;
	double current_value;
	
	x_size = image->getXSize();
	y_size = image->getYSize();
	
	boost::shared_ptr<PALMMatrix <unsigned char> > thresholded_image = boost::shared_ptr<PALMMatrix <unsigned char> >(new PALMMatrix<unsigned char>_uchar(x_size, y_size));
	
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {
			current_value = image->get(i, j);
			if (current_value >= threshold) {
				thresholded_image->set(i, j, 255);
			} else {
				thresholded_image->set(i, j, 0);
			}
		}
	}
	
	return thresholded_image;

}

boost::shared_ptr<PALMMatrix <unsigned char> > ThresholdImage_Igor_Iterative::do_thresholding(boost::shared_ptr<PALMMatrix<double> > image) {
	// we use a built-in Igor operation to handle this case
	// the approach is quite clumsy: we create a temporary wave in Igor that has the same values as the frame we are interested in
	// the igor routine then runs on that frame
	// and next we copy the data back out into a new image
	
	// this is wasteful, but in this way the interface is the same as the other classes
	size_t x_size, y_size;
	long dimensionSizes[MAX_DIMENSIONS+1];
	long indices[MAX_DIMENSIONS];
	double value[2];
	int result;
	unsigned char threshold_result;
	
	waveHndl tmp_storage_wave;
	
	x_size = image->getXSize();
	y_size = image->getYSize();
	
	dimensionSizes[0] = (long)x_size;
	dimensionSizes[1] = (long)y_size;
	dimensionSizes[2] = 0;
	
	/****** this function is not reentrant, only one thread can be past this stage. Make sure that this is so by locking the mutex ******/
	threadMutex.lock();
	
	
	// try to make the wave for the temporary storage
	result = MDMakeWave(&tmp_storage_wave, "tmp_thresh_storage", NULL, dimensionSizes, NT_FP64, 1);
	if (result != 0) {
		throw result;
	}
	
	boost::shared_ptr<PALMMatrix <unsigned char> > thresholded_image = boost::shared_ptr<PALMMatrix <unsigned char> >(new PALMMatrix<unsigned char>(x_size, y_size));
	
	// now copy the data
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {
			indices[0] = i;
			indices[1] = j;
			
			value[0] = image->get(i, j);
			
			result = MDSetNumericWavePointValue(tmp_storage_wave, indices, value) ;
			if (result != 0) {
				throw result;
			}
		}
	}
	
	// now let Igor handle the calculation of the threshold
	result = XOPSilentCommand("ImageThreshold /Q /O /M=1 tmp_thresh_storage");
	if (result != 0) {
		throw result;
	}
	
	// now copy the data back to the output matrix
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {
			indices[0] = i;
			indices[1] = j;
			
			result = MDGetNumericWavePointValue(tmp_storage_wave, indices, value) ;
			if (result != 0) {
				throw result;
			}
			
			threshold_result = (unsigned char)(value[0] + 0.5);
			
			thresholded_image->set(i, j, threshold_result);
		}
	}
	
	// now delete the temporary storage wave
	result = KillWave(tmp_storage_wave);
	if (result != 0) {
		throw result;
	}
	
	/****** release the mutex ******/
	threadMutex.unlock();
	
	return thresholded_image;
}

boost::shared_ptr<PALMMatrix <unsigned char> > ThresholdImage_Igor_Bimodal::do_thresholding(boost::shared_ptr<PALMMatrix<double> > image) {
	// we use a built-in Igor operation to handle this case
	// the approach is quite clumsy: we create a temporary wave in Igor that has the same values as the frame we are interested in
	// the igor routine then runs on that frame
	// and next we copy the data back out into a new image
	
	// this is wasteful, but the interface is the same as the other classes
	size_t x_size, y_size;
	long dimensionSizes[MAX_DIMENSIONS+1];
	long indices[MAX_DIMENSIONS];
	double value[2];
	int result;
	unsigned char threshold_result;
	
	waveHndl tmp_storage_wave;
	
	x_size = image->getXSize();
	y_size = image->getYSize();
	
	dimensionSizes[0] = (long)x_size;
	dimensionSizes[1] = (long)y_size;
	dimensionSizes[2] = 0;
	
	/****** this function is not reentrant, only one thread can be past this stage. Make sure that this is so by locking the mutex ******/
	threadMutex.lock();
	
	// try to make the wave for the temporary storage
	result = MDMakeWave(&tmp_storage_wave, "tmp_thresh_storage", NULL, dimensionSizes, NT_FP64, 1);
	if (result != 0) {
		throw result;
	}
	
	boost::shared_ptr<PALMMatrix <unsigned char> > thresholded_image = boost::shared_ptr<PALMMatrix <unsigned char> >(new PALMMatrix<unsigned char>(x_size, y_size));
	
	// now copy the data
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {
			indices[0] = i;
			indices[1] = j;
			
			value[0] = image->get(i, j);
			
			result = MDSetNumericWavePointValue(tmp_storage_wave, indices, value) ;
			if (result != 0) {
				throw result;
			}
		}
	}
	
	// now let Igor handle the calculation of the threshold
	result = XOPSilentCommand("ImageThreshold /Q /O /M=2 tmp_thresh_storage");
	if (result != 0) {
		throw result;
	}
	
	// now copy the data back to the output matrix
	for (size_t i = 0; i < x_size; ++i) {
		for (size_t j = 0; j < y_size; ++j) {
			indices[0] = i;
			indices[1] = j;
			
			result = MDGetNumericWavePointValue(tmp_storage_wave, indices, value) ;
			if (result != 0) {
				throw result;
			}
			
			threshold_result = (unsigned char)(value[0] + 0.5);
			
			thresholded_image->set(i, j, threshold_result);
		}
	}
	
	// now delete the temporary storage wave
	result = KillWave(tmp_storage_wave);
	if (result != 0) {
		throw result;
	}
	
	/****** release the mutex ******/
	threadMutex.unlock();
	
	return thresholded_image;
}

boost::shared_ptr<PALMMatrix <unsigned char> > ThresholdImage_Igor_Adaptive::do_thresholding(boost::shared_ptr<PALMMatrix<double> > image) {
	// we use a built-in Igor operation to handle this case
	// the approach is quite clumsy: we create a temporary wave in Igor that has the same values as the frame we are interested in
	// the igor routine then runs on that frame
	// and next we copy the data back out into a new image
	
	// this is wasteful, but the interface is the same as the other classes
	size_t x_size, y_size;
	long dimensionSizes[MAX_DIMENSIONS+1];
	long indices[MAX_DIMENSIONS];
	double value[2];
	int result;
	unsigned char threshold_result;
	boost::shared_ptr<PALMMatrix <unsigned char> > original_thresholded;
	boost::shared_ptr<PALMMatrix <unsigned char> > transposed_tresholded;
	
	waveHndl tmp_storage_wave;
	waveHndl thresholded_wave;
	
	x_size = image->getXSize();
	y_size = image->getYSize();
	
	// we make two images for the original and the transposed threshold
	original_thresholded = boost::shared_ptr<PALMMatrix <unsigned char> >(new PALMMatrix<unsigned char>(x_size, y_size));
	transposed_tresholded = boost::shared_ptr<PALMMatrix <unsigned char> >(new PALMMatrix<unsigned char>(x_size, y_size));
	
	dimensionSizes[0] = (long)x_size;
	dimensionSizes[1] = (long)y_size;
	dimensionSizes[2] = 0;
	
	/****** this function is not reentrant, only one thread can be past this stage. Make sure that this is so by locking the mutex ******/
	threadMutex.lock();
	
	// try to make the wave for the temporary storage
	result = MDMakeWave(&tmp_storage_wave, "tmp_thresh_storage", NULL, dimensionSizes, NT_FP64, 1);
	if (result != 0) {
		throw result;
	}
	
	// now copy the data
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {
			indices[0] = i;
			indices[1] = j;
			
			value[0] = image->get(i, j);
			
			result = MDSetNumericWavePointValue(tmp_storage_wave, indices, value) ;
			if (result != 0) {
				throw result;
			}
		}
	}
	
	// now let Igor handle the calculation of the threshold
	result = XOPSilentCommand("ImageThreshold /Q /M=3 tmp_thresh_storage");
	if (result != 0) {
		throw result;
	}
	
	
	thresholded_wave = FetchWave("M_ImageThresh");
	if (thresholded_wave == NULL) {
		throw NOWAV;
	}
	
	// now copy the thresholded image back to the first output wave
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {
			indices[0] = i;
			indices[1] = j;
			
			result = MDGetNumericWavePointValue(thresholded_wave, indices, value) ;
			if (result != 0) {
				throw result;
			}
			
			threshold_result = (unsigned char)(value[0] + 0.5);
			
			original_thresholded->set(i, j, threshold_result);
		}
	}
	
	// now transpose the image
	result = XOPSilentCommand("MatrixTranspose tmp_thresh_storage");
	if (result != 0) {
		throw result;
	}
	
	// now calculate the threshold again
	result = XOPSilentCommand("ImageThreshold /Q /M=3 tmp_thresh_storage");
	if (result != 0) {
		throw result;
	}
	
	// we transpose the thresholded image back to the original orientation
	result = XOPSilentCommand("MatrixTranspose M_ImageThresh");
	if (result != 0) {
		throw result;
	}
	
	thresholded_wave = FetchWave("M_ImageThresh");
	if (thresholded_wave == NULL) {
		throw NOWAV;
	}
	
	// now copy the thresholded image back to the second output wave
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {
			indices[0] = i;
			indices[1] = j;
			
			result = MDGetNumericWavePointValue(thresholded_wave, indices, value) ;
			if (result != 0) {
				throw result;
			}
			
			threshold_result = (unsigned char)(value[0] + 0.5);
			
			transposed_tresholded->set(i, j, threshold_result);
		}
	}
	
	// now construct the combined threshold image by AND'ing the original and transpose together
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {
			if (original_thresholded->get(i, j) < 128) {	// below the threshold, we skip it
				continue;
			} else {
				// is the transposed matrix also above the threshold?
				if (transposed_tresholded->get(i, j) < 128) {	// below the threshold. We should not include this point
					original_thresholded->set(i, j, 0);
				} else {
					continue;
				}
			}
		}
	}
	
	
	
	// now delete the temporary storage wave
	result = KillWave(tmp_storage_wave);
	if (result != 0) {
		throw result;
	}
	
	/****** release the mutex ******/
	threadMutex.unlock();
	
	return original_thresholded;
}

boost::shared_ptr<PALMMatrix <unsigned char> > ThresholdImage_Igor_Fuzzy1::do_thresholding(boost::shared_ptr<PALMMatrix<double> > image) {
	// we use a built-in Igor operation to handle this case
	// the approach is quite clumsy: we create a temporary wave in Igor that has the same values as the frame we are interested in
	// the igor routine then runs on that frame
	// and next we copy the data back out into a new image
	
	// this is wasteful, but the interface is the same as the other classes
	size_t x_size, y_size;
	long dimensionSizes[MAX_DIMENSIONS+1];
	long indices[MAX_DIMENSIONS];
	double value[2];
	int result;
	unsigned char threshold_result;
	
	waveHndl tmp_storage_wave;
	
	x_size = image->getXSize();
	y_size = image->getYSize();
	
	dimensionSizes[0] = (long)x_size;
	dimensionSizes[1] = (long)y_size;
	dimensionSizes[2] = 0;
	
	/****** this function is not reentrant, only one thread can be past this stage. Make sure that this is so by locking the mutex ******/
	threadMutex.lock();
	
	// try to make the wave for the temporary storage
	result = MDMakeWave(&tmp_storage_wave, "tmp_thresh_storage", NULL, dimensionSizes, NT_FP64, 1);
	if (result != 0) {
		throw result;
	}
	
	boost::shared_ptr<PALMMatrix <unsigned char> > thresholded_image = boost::shared_ptr<PALMMatrix <unsigned char> >(new PALMMatrix<unsigned char>(x_size, y_size));
	
	// now copy the data
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {
			indices[0] = i;
			indices[1] = j;
			
			value[0] = image->get(i, j);
			
			result = MDSetNumericWavePointValue(tmp_storage_wave, indices, value) ;
			if (result != 0) {
				throw result;
			}
		}
	}
	
	// now let Igor handle the calculation of the threshold
	result = XOPSilentCommand("ImageThreshold /Q /O /M=4 tmp_thresh_storage");
	if (result != 0) {
		throw result;
	}
	
	// now copy the data back to the output matrix
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {
			indices[0] = i;
			indices[1] = j;
			
			result = MDGetNumericWavePointValue(tmp_storage_wave, indices, value) ;
			if (result != 0) {
				throw result;
			}
			
			threshold_result = (unsigned char)(value[0] + 0.5);
			
			thresholded_image->set(i, j, threshold_result);
		}
	}
	
	// now delete the temporary storage wave
	result = KillWave(tmp_storage_wave);
	if (result != 0) {
		throw result;
	}
	
	/****** release the mutex ******/
	threadMutex.unlock();
	
	return thresholded_image;	
}

boost::shared_ptr<PALMMatrix <unsigned char> > ThresholdImage_Igor_Fuzzy2::do_thresholding(boost::shared_ptr<PALMMatrix<double> > image) {
	// we use a built-in Igor operation to handle this case
	// the approach is quite clumsy: we create a temporary wave in Igor that has the same values as the frame we are interested in
	// the igor routine then runs on that frame
	// and next we copy the data back out into a new image
	
	// this is wasteful, but the interface is the same as the other classes
	size_t x_size, y_size;
	long dimensionSizes[MAX_DIMENSIONS+1];
	long indices[MAX_DIMENSIONS];
	double value[2];
	int result;
	unsigned char threshold_result;
	
	waveHndl tmp_storage_wave;
	
	x_size = image->getXSize();
	y_size = image->getYSize();
	
	dimensionSizes[0] = (long)x_size;
	dimensionSizes[1] = (long)y_size;
	dimensionSizes[2] = 0;
	
	/****** this function is not reentrant, only one thread can be past this stage. Make sure that this is so by locking the mutex ******/
	threadMutex.lock();
	
	// try to make the wave for the temporary storage
	result = MDMakeWave(&tmp_storage_wave, "tmp_thresh_storage", NULL, dimensionSizes, NT_FP64, 1);
	if (result != 0) {
		throw result;
	}
	
	boost::shared_ptr<PALMMatrix <unsigned char> > thresholded_image = boost::shared_ptr<PALMMatrix <unsigned char> >(new PALMMatrix<unsigned char>(x_size, y_size));
	
	// now copy the data
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {
			indices[0] = i;
			indices[1] = j;
			
			value[0] = image->get(i, j);
			
			result = MDSetNumericWavePointValue(tmp_storage_wave, indices, value) ;
			if (result != 0) {
				throw result;
			}
		}
	}
	
	// now let Igor handle the calculation of the threshold
	result = XOPSilentCommand("ImageThreshold /Q /O /M=5 tmp_thresh_storage");
	if (result != 0) {
		throw result;
	}
	
	// now copy the data back to the output matrix
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {
			indices[0] = i;
			indices[1] = j;
			
			result = MDGetNumericWavePointValue(tmp_storage_wave, indices, value) ;
			if (result != 0) {
				throw result;
			}
			
			threshold_result = (unsigned char)(value[0] + 0.5);
			
			thresholded_image->set(i, j, threshold_result);
		}
	}
	
	// now delete the temporary storage wave
	result = KillWave(tmp_storage_wave);
	if (result != 0) {
		throw result;
	}
	
	/****** release the mutex ******/
	threadMutex.unlock();
	
	return thresholded_image;
}

boost::shared_ptr<PALMMatrix <unsigned char> > ThresholdImage_Isodata::do_thresholding(boost::shared_ptr<PALMMatrix<double> > image) {
	gsl_histogram *hist;
	boost::shared_ptr<PALMMatrix <unsigned char> > threshold_image;
	
	size_t x_size = image->getXSize();
	size_t y_size = image->getYSize();
	
	size_t number_of_bins = 256;
	size_t current_threshold_bin = 127;
	int n_iterations = 0;
	int max_iters = 50;
	double lower_mean, upper_mean;
	double sum, denominator_sum;
	double current_threshold = -1, previous_threshold;
	double lower_bin_limit, upper_bin_limit;
	size_t bin_threshold;
	double intensity_threshold;
	int converged = 0;
	
	// since this is a histogram-based approach we start by constructing the histogram
	hist = make_histogram_from_matrix(image, number_of_bins);
	
	// because this approach is based on thresholding it makes sense to only express the threshold in terms of bins, stored in "current_threshold_bin".
	// a value of 0 for "current_threshold_bin" means that all bins at index 0 and higher are considered to be 'signal', not 'background'.
	
	while ((converged == 0) && (n_iterations < max_iters)) {
		
		previous_threshold = current_threshold;
		
		// calculate the lower mean
		sum = 0;
		denominator_sum = 0;
		for (size_t i = 0; i < current_threshold_bin; i++) {
			sum += (double)i * gsl_histogram_get(hist, i);
			denominator_sum += gsl_histogram_get(hist, i);
		}
		
		lower_mean = sum / denominator_sum;
		
		// calculate the upper mean
		sum = 0;
		denominator_sum = 0;
		for (size_t i = current_threshold_bin; i < number_of_bins; i++) {
			sum += (double)i * gsl_histogram_get(hist, i);
			denominator_sum += gsl_histogram_get(hist, i);
		}
		
		upper_mean = sum / denominator_sum;
		
		current_threshold = (lower_mean + upper_mean) / 2.0;
		current_threshold_bin = floor(current_threshold);
		
		if (floor(current_threshold + 0.5) == floor(previous_threshold + 0.5)) {
			bin_threshold = floor(current_threshold + 0.5);
			converged = 1;
		}
		
	}
	
	if (converged == 0) {	// the iterations did not converge, there is no clear threshold
		// to indicate this we set everything to 'off' (0)
		threshold_image = boost::shared_ptr<PALMMatrix <unsigned char> >(new PALMMatrix<unsigned char>(x_size, y_size));
		threshold_image->set_all(0);
		return threshold_image;
	}
	
	// now translate the threshold value to an intensity instead of being in bins
	gsl_histogram_get_range(hist, bin_threshold, &lower_bin_limit, &upper_bin_limit);
	intensity_threshold = lower_bin_limit;
	
	// get another threshold class to do the work for us
	ThresholdImage_Direct thresholder(intensity_threshold);
	
	try {
		threshold_image = thresholder.do_thresholding(image);
	}
	catch (OUT_OF_MEMORY) {
		gsl_histogram_free(hist);
		string error;
		error = "unable to allocate threshold_image in ThresholdImage_Isodata::do_thresholding()\r";
		throw OUT_OF_MEMORY(error);
	}
	
	gsl_histogram_free(hist);
	
	return threshold_image;
}

boost::shared_ptr<PALMMatrix <unsigned char> > ThresholdImage_Triangle::do_thresholding(boost::shared_ptr<PALMMatrix<double> > image) {
	gsl_histogram *hist;
	size_t number_of_bins = 256;
	size_t maximum_bin;
	double max_val, max_bin_double;
	double end_val, end_bin_double;
	double slope, intercept, perpendicular_slope, perpendicular_intercept;
	double current_bin_value, double_i;
	double intercept_x, intercept_y;
	double distance;
	double max_distance = -1;
	size_t max_index;
	double lower_bin_limit, upper_bin_limit, intensity_threshold;
	
	boost::shared_ptr<PALMMatrix <unsigned char> > threshold_image;
	size_t x_size = image->getXSize();
	size_t y_size = image->getYSize();
	
	// since this is a histogram-based approach we start by constructing the histogram
	
	try {
		hist = make_histogram_from_matrix(image, number_of_bins);
	}
	catch (OUT_OF_MEMORY) {
		string error;
		error = "unable to allocate buffer in ThresholdImage_Triangle::do_thresholding()\r";
		throw OUT_OF_MEMORY(error);
	}
	
	maximum_bin = gsl_histogram_max_bin(hist);
	max_bin_double = (double)maximum_bin;
	max_val = gsl_histogram_max_val(hist);
	end_val = gsl_histogram_get(hist, number_of_bins - 1);	// the bin that contains the largest intensity is the last bin in the histogram
	end_bin_double = (double)(number_of_bins - 1);
	
	// catch an unlikely case where the maximum corresponds to the last bin
	if (maximum_bin == (number_of_bins - 1)) {
		gsl_histogram_free(hist);
		threshold_image = boost::shared_ptr<PALMMatrix <unsigned char> >(new PALMMatrix<unsigned char>(x_size, y_size));
		threshold_image->set_all(0);
		return threshold_image;
	}
	
	// calculate the line that connects the maximum and highest-intensity value
	slope = (end_val - max_val) / (end_bin_double - max_bin_double);
	intercept = max_val / (slope * max_bin_double);
	
	// catch an unlikely case where the connecting line is flat (the histogram is apparently uniform)
	if (slope == 0) {
		gsl_histogram_free(hist);
		threshold_image = boost::shared_ptr<PALMMatrix <unsigned char> >(new PALMMatrix<unsigned char>(x_size, y_size));
		threshold_image->set_all(0);
		return threshold_image;
	}
	
	
	// calculate the slope of a line perpendicular to the connecting line
	perpendicular_slope = - 1.0 / slope;
	
	for (size_t i = maximum_bin + 1; i < number_of_bins; i++) {	// determine the minimum distance in the triangle
		
		// what is the intercept for the perpendicular line if it has to go through the bin that we're currently looking at?
		current_bin_value = gsl_histogram_get(hist, i);
		double_i = (double)i;
		
		perpendicular_intercept = current_bin_value / (perpendicular_slope * double_i);
		
		// where does the perpendicular line intercept the connecting line?
		// x = (b1 - b2) / (a1 - a2) and y = a1 * x + b1
		intercept_x = (intercept - perpendicular_intercept) / (slope - perpendicular_slope);
		intercept_y = slope * intercept_x + intercept;
		
		// what is the distance to the connecting line?
		distance = sqrt((intercept_x - double_i) * (intercept_x - double_i) + (intercept_y - current_bin_value) * (intercept_y - current_bin_value));
		
		if (distance > max_distance) {
			max_index = i;
			max_distance = distance;
		}
	}
	
	// translate the maximal index to a threshold value
	gsl_histogram_get_range(hist, max_index, &lower_bin_limit, &upper_bin_limit);
	intensity_threshold = lower_bin_limit;
	
	// get another threshold class to do the work for us
	ThresholdImage_Direct thresholder(intensity_threshold);
	
	try {
		threshold_image = thresholder.do_thresholding(image);
	}
	catch (OUT_OF_MEMORY) {
		string error;
		error = "unable to do direct thresholding in ThresholdImage_Triangle::do_thresholding()\r";
		gsl_histogram_free(hist);
		throw OUT_OF_MEMORY(error);
	}
	
	gsl_histogram_free(hist);
	
	return threshold_image;
}


boost::shared_ptr<PALMMatrix <unsigned char> > ThresholdImage_GLRT::do_thresholding(boost::shared_ptr<PALMMatrix<double> > image) {
	// the code is based on a series of matlab files sent by Didier Marguet, corresponding author of the original paper
	boost::shared_ptr<PALMMatrix <unsigned char> > threshold_image;
	size_t x_size = image->getXSize();
	size_t y_size = image->getYSize();
	
	boost::shared_ptr<PALMMatrix<double> > averages(new PALMMatrix<double> (x_size, y_size));	// contains an estimation of the average value at every position, calculation over the number of pixels in the window
	boost::shared_ptr<PALMMatrix<double> > image_squared(new PALMMatrix<double> (x_size, y_size));
	boost::shared_ptr<PALMMatrix<double> > summed_squares(new PALMMatrix<double> (x_size, y_size));
	boost::shared_ptr<PALMMatrix<double> > null_hypothesis(new PALMMatrix<double> (x_size, y_size));
	boost::shared_ptr<PALMMatrix<double> > Gaussian_window(new PALMMatrix<double> (x_size, y_size));
	boost::shared_ptr<PALMMatrix<double> > image_Gaussian_convolved(new PALMMatrix<double> (x_size, y_size));	// this is 'alpha' in the original matlab code
	boost::shared_ptr<PALMMatrix<double> > hypothesis_test(new PALMMatrix<double> (x_size, y_size));	// this is 'test' in the original matlab code
	
	
	double average = 0;
	size_t window_size = 13;	// this is the size of the window over which we calculate the hypotheses
	size_t half_window_size = window_size / 2;	// integer division takes care of the floor() aspect
	size_t window_pixels = window_size * window_size;
	double double_window_pixels = (double)(window_pixels);
	double current_value;
	double distance_x, distance_y;
	double sum;
	double sum_squared_Gaussian;
	
	
	threshold_image = boost::shared_ptr<PALMMatrix <unsigned char> >(new PALMMatrix<unsigned char>(x_size, y_size));
	threshold_image->set_all(0);
	
	averages->set_all(0);
	summed_squares->set_all(0);
	hypothesis_test->set_all(0);
	
	// calculate the square of the pixel values
	// we'll use this later
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {
			current_value = image->get(i, j);
			current_value = current_value * current_value;
			image_squared->set(i, j, current_value);
		}
	}
	
	// NULL HYPOTHESIS: there is no emitter at a certain position
	
	// start by estimating the mean at every position
	// in the original code this done by a convolution of a unity matrix with the window size-> This is done using an FFT in the matlab code, but we'll do it directly
	
	for (size_t k = half_window_size; k < x_size - half_window_size; k++) {
		for (size_t l = half_window_size; l < y_size - half_window_size; l++) {
			// these two loops run over the entire image-> We exclude the edges where the window doesn't fit
			
			// now we loop over the size of the window, to determine the average at every pixel
			average = 0;
			
			for (size_t j = l - half_window_size; j <= l + half_window_size; j++) {
				for (size_t i = k - half_window_size; i <= k + half_window_size; i++) {
					average += image->get(i, j);
				}
			}
			average /= double_window_pixels;
			
			averages->set(k, l, average);
		}
	}
	
	// now we do the same, but for the square of the pixel values, and we don't divide (no averaging)
	for (size_t k = half_window_size; k < x_size - half_window_size; k++) {
		for (size_t l = half_window_size; l < y_size - half_window_size; l++) {
			// these two loops run over the entire image-> We exclude the edges where the window doesn't fit
			current_value = 0;
			for (size_t i = k - half_window_size; i <= k + half_window_size; i++) {
				for (size_t j = l - half_window_size; j <= l + half_window_size; j++) {
					current_value += image_squared->get(i, j);
				}
			}
			summed_squares->set(k, l, current_value);	// summed_squares is now equal to "Sim2" in the orignal matlab code
		}
	}
	
	// now calculate the null hypothesis image-> This is T_sig0_2 in the original matlab source
	for (size_t l = half_window_size; l < y_size - half_window_size; l++) {
		for (size_t k = half_window_size; k < x_size - half_window_size; k++) {
			current_value = summed_squares->get(k, l) - double_window_pixels * averages->get(k, l) * averages->get(k, l);
			null_hypothesis->set(k, l, current_value);
		}
	}
	
	// calculate the hypothesis H1 that there is an emitter
	
	// store the values of a Gaussian in a matrix with the size of a window
	// this is so we can cache the values, and do not to calculate it every time
	sum = 0;
	for (size_t j = 0; j < window_size; j++) {
		for (size_t i = 0; i < window_size; i++) {
			// the Gaussian is assumed to be in the center of the window
			distance_x = (double)half_window_size - (double)i;
			distance_y = (double)half_window_size - (double)j;
			current_value = 1.0 / (1.77245385 * gaussianWidth) * exp(- 1.0 / (2.0 * gaussianWidth * gaussianWidth) * (distance_x * distance_x + distance_y * distance_y));
			
			Gaussian_window->set(i, j, current_value);
			
			sum += current_value;	// we will use this below
		}
	}
	
	// now we re-normalize this Gaussian matrix
	// at this point Gaussian_window becomes equal to 'gc' in the original matlab code
	sum /= double_window_pixels;
	sum_squared_Gaussian = 0;
	for (size_t j = 0; j < window_size; j++) {
		for (size_t i = 0; i < window_size; i++) {
			current_value = Gaussian_window->get(i, j);
			current_value = current_value - sum;
			Gaussian_window->set(i, j, current_value);
			sum_squared_Gaussian += current_value * current_value;	// this is 'Sgc2' in the original code
		}
	}
	
	// now we need to again convolve this Gaussian_window ('gc') with the original image-> As before, we'll do it directly
	// these two loops run over the entire image-> We exclude the edges where the window doesn't fit
	for (size_t l = half_window_size; l < y_size - half_window_size; l++) {
		for (size_t k = half_window_size; k < x_size - half_window_size; k++) {
			
			sum = 0;
			for (size_t j = l - half_window_size; j <= l + half_window_size; j++) {
				for (size_t i = k - half_window_size; i <= k + half_window_size; i++) {
					current_value = Gaussian_window->get(i - k + half_window_size, j - l + half_window_size);
					current_value *= image->get(i, j);
					sum += current_value;
				}
			}
			sum /= sum_squared_Gaussian;
			image_Gaussian_convolved->set(k, l, sum);
		}
	}
	
	// "image_Gaussian_convolved" is now equal to 'alpha' in the original matlab code
	
	// calculate the image that will determine whether to accept or reject the null hypothesis
	for (size_t l = half_window_size; l < y_size - half_window_size; l++) {
		for (size_t k = half_window_size; k < x_size - half_window_size; k++) {
			current_value = 1 - (sum_squared_Gaussian * image_Gaussian_convolved->get(k, l) * image_Gaussian_convolved->get(k, l)) / null_hypothesis->get(k , l);
			current_value = (current_value > 0) * current_value + (current_value <= 0);	// the equivalent of test = (test > 0) ->* test + (test <= 0) in the original code
			current_value = - double_window_pixels * log(current_value);
			hypothesis_test->set(k, l, current_value);
		}
	}
	
	// at this point 'hypothesis_test' is equal to 'carte_MV' in the original image
	// check where we have to reject the hypothesis test
	for (size_t l = half_window_size; l < y_size - half_window_size; l++) {
		for (size_t k = half_window_size; k < x_size - half_window_size; k++) {
			if (hypothesis_test->get(k, l) > PFA) {
				threshold_image->set(k, l, 255);
			}
		}
	}
	
	return threshold_image;
}


boost::shared_ptr<PALMMatrix <unsigned char> > ThresholdImage_GLRT_FFT::do_thresholding(boost::shared_ptr<PALMMatrix<double> > image) {
	// the code is based on a series of matlab files sent by Didier Marguet, corresponding author of the original paper
	boost::shared_ptr<PALMMatrix <unsigned char> > threshold_image;
	size_t x_size = image->getXSize();
	size_t y_size = image->getYSize();
	
	boost::shared_ptr<PALMMatrix<double> > averages;
	boost::shared_ptr<PALMMatrix<double> > image_squared;
	boost::shared_ptr<PALMMatrix<double> > summed_squares;
	boost::shared_ptr<PALMMatrix<double> > null_hypothesis;
	boost::shared_ptr<PALMMatrix<double> > image_Gaussian_convolved;
	boost::shared_ptr<PALMMatrix<double> > hypothesis_test;
	boost::shared_ptr<PALMMatrix<double> > Gaussian_window;
	
	size_t window_size = 13;
	size_t half_window_size = window_size / 2;	// integer division takes care of the floor() aspect
	size_t center_x = x_size / 2;
	size_t center_y = y_size / 2;
	size_t window_pixels = window_size * window_size;
	double double_window_pixels = (double)(window_pixels);
	double current_value;
	double distance_x, distance_y;
	double sum;
	
	threshold_image = boost::shared_ptr<PALMMatrix <unsigned char> >(new PALMMatrix<unsigned char>(x_size, y_size));
	threshold_image->set_all(0);
	
	averages = boost::shared_ptr<PALMMatrix<double> >(new PALMMatrix<double>(x_size, y_size));
	averages->set_all(0);
	
	image_squared = boost::shared_ptr<PALMMatrix<double> >(new PALMMatrix<double>(x_size, y_size));
	
	summed_squares = boost::shared_ptr<PALMMatrix<double> >(new PALMMatrix<double>(x_size, y_size));
	summed_squares->set_all(0);
	
	null_hypothesis = boost::shared_ptr<PALMMatrix<double> >(new PALMMatrix<double>(x_size, y_size));
	
	Gaussian_window = boost::shared_ptr<PALMMatrix<double> >(new PALMMatrix<double>(x_size, y_size));
	
	image_Gaussian_convolved = boost::shared_ptr<PALMMatrix<double> >(new PALMMatrix<double>(x_size, y_size));	// this is 'alpha' in the original matlab code
	
	hypothesis_test = boost::shared_ptr<PALMMatrix<double> >(new PALMMatrix<double>(x_size, y_size));	// this is 'test' in the original matlab code
	hypothesis_test->set_all(0);
	
	// calculate the square of the pixel values
	// we'll use this later
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {
			current_value = image->get(i, j);
			current_value = current_value * current_value;
			image_squared->set(i, j, current_value);
		}
	}
	
	// NULL HYPOTHESIS: there is no emitter at a certain position
	
	// start by estimating the mean at every position
	// in the original code this done by a convolution of a unity matrix with the window size. We now do this using an FFT-based approach
	
	
	// convolve the image with a "box function", that will get us the average
	// if we have a lot of images then we only need to make this kernel once
	// but we need to make sure that only a single thread at a time is making the kernel, the others should block
	AverageKernelMutex.lock();
	
	if ((average_kernel.get() == NULL) || (averageKernelXSize != x_size) || (averageKernelYSize != y_size)) {
		averageCalculationMutex.lock();	// get a unique lock
		
		average_kernel = boost::shared_ptr<PALMMatrix<double> >(new PALMMatrix<double>(x_size, y_size));
		
		average_kernel->set_all(0);
		
		for (size_t j = center_y - half_window_size; j <= center_y + half_window_size; j++) {
			for (size_t i = center_x - half_window_size; i <= center_x + half_window_size; i++) {
				average_kernel->set(i, j, 1);
			}
		}
		
		averageKernelXSize = x_size;
		averageKernelYSize = y_size;
		
		averageCalculationMutex.unlock();
	}
	
	averageCalculationMutex.lock_shared();
	AverageKernelMutex.unlock();
	
	averages = matrixConvolver.ConvolveMatricesWithFFT(image, average_kernel);
	
	// do the same for the squared image
	summed_squares = matrixConvolver.ConvolveMatricesWithFFT(image_squared, average_kernel);
	
	averageCalculationMutex.unlock_shared();
	
	// normalize the result, so that we get averages
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {
			current_value = averages->get(i, j);
			current_value /= double_window_pixels;
			averages->set(i, j, current_value);
		}
	}
	
	// now calculate the null hypothesis image. This is T_sig0_2 in the original matlab source
	for (size_t l = half_window_size; l < y_size - half_window_size; l++) {
		for (size_t k = half_window_size; k < x_size - half_window_size; k++) {
			current_value = summed_squares->get(k, l) - double_window_pixels * averages->get(k, l) * averages->get(k, l);
			null_hypothesis->set(k, l, current_value);
		}
	}
	
	// calculate the hypothesis H1 that there is an emitter
	
	// create a Gaussian kernel for the convolution
	// we only need to do this once if we are looking at a series of images
	// but make sure that only a single thread is doing this
	GaussianKernelMutex.lock();
	
	if ((Gaussian_kernel.get() == NULL) || (GaussianKernelXSize != x_size) || (GaussianKernelYSize != y_size)) {
		gaussianCalculationMutex.lock();	// get a unique lock
		Gaussian_kernel = boost::shared_ptr<PALMMatrix<double> >(new PALMMatrix<double>(x_size, y_size));
		
		Gaussian_window = boost::shared_ptr<PALMMatrix<double> >(new PALMMatrix<double>(window_size, window_size));
		
		sum = 0;
		for (size_t j = 0; j < window_size; j++) {
			for (size_t i = 0; i < window_size; i++) {
				// the Gaussian is assumed to be in the center of the window
				distance_x = (double)half_window_size - (double)i;
				distance_y = (double)half_window_size - (double)j;
				current_value = 1.0 / (1.77245385 * gaussianWidth) * exp(- 1.0 / (2.0 * gaussianWidth * gaussianWidth) * (distance_x * distance_x + distance_y * distance_y));
				
				Gaussian_window->set(i, j, current_value);
				
				sum += current_value;	// we will use this below
			}
		}
		
		// now we re-normalize this Gaussian matrix
		// at this point Gaussian_window becomes equal to 'gc' in the original matlab code
		sum /= double_window_pixels;
		sum_squared_Gaussian = 0;
		for (size_t j = 0; j < window_size; j++) {
			for (size_t i = 0; i < window_size; i++) {
				current_value = Gaussian_window->get(i, j);
				current_value = current_value - sum;
				Gaussian_window->set(i, j, current_value);
				sum_squared_Gaussian += current_value * current_value;	// this is 'Sgc2' in the original code
			}
		}
		
		// now introduce this small kernel into a larger one that is the same size as the image
		Gaussian_kernel->set_all(0);
		
		for (size_t j = center_y - half_window_size; j <= center_y + half_window_size; j++) {
			for (size_t i = center_x - half_window_size; i <= center_x + half_window_size; i++) {
				Gaussian_kernel->set(i, j, Gaussian_window->get(i - center_x + half_window_size, j - center_y + half_window_size));
			}
		}
		
		GaussianKernelXSize = x_size;
		GaussianKernelYSize = y_size;
		
		gaussianCalculationMutex.unlock();
	}
	
	// now we need to again convolve this Gaussian_window ('gc') with the original image. 
	// we now do this using the FFT
	
	gaussianCalculationMutex.lock_shared();
	GaussianKernelMutex.unlock();
	
	image_Gaussian_convolved = matrixConvolver.ConvolveMatricesWithFFT(image, Gaussian_kernel);
	gaussianCalculationMutex.unlock_shared();
	
	// now normalize this convolved image so that it becomes equal to 'alpha' in the original matlab code
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {
			current_value = image_Gaussian_convolved->get(i, j);
			current_value /= sum_squared_Gaussian;
			image_Gaussian_convolved->set(i, j, current_value);
		}
	}
	
	// calculate the image that will determine whether to accept or reject the null hypothesis
	for (size_t l = half_window_size; l < y_size - half_window_size; l++) {
		for (size_t k = half_window_size; k < x_size - half_window_size; k++) {
			current_value = 1 - (sum_squared_Gaussian * image_Gaussian_convolved->get(k, l) * image_Gaussian_convolved->get(k, l)) / null_hypothesis->get(k , l);
			current_value = (current_value > 0) * current_value + (current_value <= 0);	// the equivalent of test = (test > 0) .* test + (test <= 0) in the original code
			current_value = - double_window_pixels * log(current_value);
			hypothesis_test->set(k, l, current_value);
		}
	}
	
	// at this point 'hypothesis_test' is equal to 'carte_MV' in the original image
	// check where we have to reject the hypothesis test
	for (size_t l = half_window_size; l < y_size - half_window_size; l++) {
		for (size_t k = half_window_size; k < x_size - half_window_size; k++) {
			if (hypothesis_test->get(k, l) > PFA) {
				threshold_image->set(k, l, 255);
			}
		}
	}
	
	return threshold_image;
}


boost::shared_ptr<PALMMatrix<double> > ThresholdImage_Preprocessor_MedianFilter::do_preprocessing(boost::shared_ptr<PALMMatrix<double> > image) {
	
	size_t kernel_size = kernel_x_size * kernel_y_size;
	size_t half_kernel_size_x = kernel_x_size / 2;
	size_t half_kernel_size_y = kernel_y_size / 2;
	size_t x_size = image->getXSize();
	size_t y_size = image->getYSize();
	size_t offset;
	double value, median;
	
	gsl_vector *median_environment;
	boost::shared_ptr<PALMMatrix<double> > filtered_image;
	
	// allocate a gsl_vector with the correct size
	median_environment = gsl_vector_alloc(kernel_size);
	size_t sorted_center = kernel_size / 2;
	
	// make a copy of the image
	// this copy will be median-filtered
	// close to the edges (where the kernel doesn't fit we will not modify the image)
	filtered_image = boost::shared_ptr<PALMMatrix<double> >(new PALMMatrix<double>(x_size, y_size));
	
	for (size_t j = 0; j < y_size; ++j) {
		for (size_t i = 0; i < x_size; ++i) {
			filtered_image->set(i, j, image->get(i, j));
		}
	}
	
	
	// the main loop
	// for now we handle this using a direct, slow calculation
	
	for (size_t j = half_kernel_size_y; j < y_size - half_kernel_size_y; j++) {
		for (size_t i = half_kernel_size_x; i < x_size - half_kernel_size_x; i++) {
			
			offset = 0;
			for (size_t l = j - half_kernel_size_y; l <= j + half_kernel_size_y; l++) {
				for (size_t k = i - half_kernel_size_x; k <= i + half_kernel_size_x; k++) {
					value = image->get(k, l);
					gsl_vector_set(median_environment, offset, value);
					offset++;
				}
			}
			gsl_sort_vector(median_environment);
			median = gsl_vector_get(median_environment, sorted_center);
			filtered_image->set(i, j, median);
		}
	}
	
	gsl_vector_free(median_environment);
	
	return filtered_image;
}


void ThresholdImage_Preprocessor_GaussianSmoothing::generate_Gaussian_kernel(size_t x_size, size_t y_size) {
	
	size_t window_size = 31;
	size_t half_window_size = window_size / 2;
	size_t center_x = x_size / 2;
	size_t center_y = y_size / 2;
	double current_value, distance_x, distance_y;
	
	boost::shared_ptr<PALMMatrix<double> > Gaussian_window(new PALMMatrix<double>(window_size, window_size));
	
	Gaussian_kernel = boost::shared_ptr<PALMMatrix<double> >(new PALMMatrix<double>(x_size, y_size));
	
	
	// calculate the values of a Gaussian with the correct width in a smaller window
	for (size_t j = 0; j < window_size; j++) {
		for (size_t i = 0; i < window_size; i++) {
			// the Gaussian is assumed to be in the center of the window
			distance_x = (double)half_window_size - (double)i;
			distance_y = (double)half_window_size - (double)j;
			current_value = 1.0 / (6.28318531 * width * width) * exp(- 1.0 / (2.0 * width * width) * (distance_x * distance_x + distance_y * distance_y));
			// normalized Gaussian in two dimensions
			
			Gaussian_window->set(i, j, current_value);
		}
	}
	
	// now introduce this small kernel into a larger one that is the same size as the image
	Gaussian_kernel->set_all(0);
	
	for (size_t j = center_y - half_window_size; j <= center_y + half_window_size; j++) {
		for (size_t i = center_x - half_window_size; i <= center_x + half_window_size; i++) {
			Gaussian_kernel->set(i, j, Gaussian_window->get(i - center_x + half_window_size, j - center_y + half_window_size));
		}
	}
}
	
	


boost::shared_ptr<PALMMatrix<double> > ThresholdImage_Preprocessor_GaussianSmoothing::do_preprocessing(boost::shared_ptr<PALMMatrix<double> > image) {
	
	boost::shared_ptr<PALMMatrix<double> > filtered_image;
	size_t x_size = image->getXSize();
	size_t y_size = image->getYSize();
	
	// do we already have a Gaussian kernel stored, or is this the first run?
	generateKernelMutex.lock();
	if (Gaussian_kernel.get() == NULL) {	// we don't have a kernel, we need to generate it
		
		generate_Gaussian_kernel(x_size, y_size);
		
	} else {	// we already have a kernel stored, is it the correct size?
				// if not we will calculate a new one
		if ((x_size != Gaussian_kernel->getXSize()) || (y_size != Gaussian_kernel->getYSize())) {
			generate_Gaussian_kernel(x_size, y_size);
		}
	}
	generateKernelMutex.unlock();
	
	filtered_image = matrixConvolver.ConvolveMatricesWithFFT(image, Gaussian_kernel);
	
	return filtered_image;
}


boost::shared_ptr<PALMMatrix<double> > ThresholdImage_Preprocessor_MeanFilter::do_preprocessing(boost::shared_ptr<PALMMatrix<double> > image) {
	
	size_t kernel_size = kernel_x_size * kernel_y_size;
	double double_kernel_pixels = (double)kernel_size;
	size_t half_kernel_size_x = kernel_x_size / 2;
	size_t half_kernel_size_y = kernel_y_size / 2;
	size_t x_size = image->getXSize();
	size_t y_size = image->getYSize();
	double mean;
	
	boost::shared_ptr<PALMMatrix<double> > filtered_image;
	
	// make a copy of the image
	// this copy will be mean-filtered
	// close to the edges, where the kernel doesn't fit we will not modify the image
	filtered_image = boost::shared_ptr<PALMMatrix<double> >(new PALMMatrix<double>(x_size, y_size));
	
	for (size_t i = 0; i < x_size; ++i) {
		for (size_t j = 0; j < y_size; ++j) {
			filtered_image->set(i, j, image->get(i, j));
		}
	}
	
	
	// the main loop
	// for now we handle this using a direct, slow calculation
	
	for (size_t j = half_kernel_size_y; j < y_size - half_kernel_size_y; j++) {
		for (size_t i = half_kernel_size_x; i < x_size - half_kernel_size_x; i++) {
			
			mean = 0;
			for (size_t l = j - half_kernel_size_y; l <= j + half_kernel_size_y; l++) {
				for (size_t k = i - half_kernel_size_x; k <= i + half_kernel_size_x; k++) {
					mean += image->get(k, l);
				}
			}
			mean /= double_kernel_pixels;
			filtered_image->set(i, j, mean);
		}
	}
	
	return filtered_image;
}
	
			
boost::shared_ptr<PALMMatrix <unsigned char> > ThresholdImage_Postprocessor_RemoveIsolatedPixels::do_postprocessing(boost::shared_ptr<PALMMatrix <unsigned char> > thresholded_image, boost::shared_ptr<PALMMatrix<double> > image) {
	// we don't care about the edges, they are ignored anyway in the fitting
	size_t x_size = thresholded_image->getXSize();
	size_t y_size = thresholded_image->getYSize();
	unsigned char value;
	double meanIntensity = 0;
	
	boost::shared_ptr<PALMMatrix <unsigned char> > processed_thresholded_image;
	
	processed_thresholded_image = boost::shared_ptr<PALMMatrix <unsigned char> >(new PALMMatrix<unsigned char>(x_size, y_size));
	
	processed_thresholded_image->set_all(0);
	
	// calculate the mean intensity
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {
			meanIntensity += image->get(i, j);
		}
	}
	meanIntensity /= x_size * y_size;
	
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {
			
			value = thresholded_image->get(i, j);
			if (value < 128) {	// this is an 'off' pixel
				continue;
			}
			
			// if we are here then the current pixel is active
			// we only mark the pixel as active in the processed image if it is above the mean
			if (image->get(i, j) > meanIntensity) {
				processed_thresholded_image->set(i, j, 255);
			}
		}
	}
	
	return processed_thresholded_image;
}


boost::shared_ptr<PALMMatrix <unsigned char> > ThresholdImage_Postprocessor_RemovePixelsBelowMean::do_postprocessing(boost::shared_ptr<PALMMatrix <unsigned char> > thresholded_image, boost::shared_ptr<PALMMatrix<double> > image) {
	// we don't care about the edges, they are ignored anyway in the fitting
	size_t x_size = thresholded_image->getXSize();
	size_t y_size = thresholded_image->getYSize();
	unsigned char value;
	bool neighbour_found;
	
	boost::shared_ptr<PALMMatrix <unsigned char> > processed_thresholded_image;
	
	processed_thresholded_image = boost::shared_ptr<PALMMatrix <unsigned char> >(new PALMMatrix<unsigned char>(x_size, y_size));
	
	for (size_t i = 0; i < x_size; ++i) {
		for (size_t j = 0; j < y_size; ++j) {
			processed_thresholded_image->set(i, j, thresholded_image->get(i, j));
		}
	}
	
	// we will return a copy
	
	for (size_t j = 1; j < y_size - 1; j++) {
		for (size_t i = 1; i < x_size - 1; i++) {
			
			value = thresholded_image->get(i, j);
			if (value < 128) {	// this is an 'off' pixel
				continue;
			}
			
			neighbour_found = 0;
			for (size_t l = j - 1; l <= j + 1; l++) {
				for (size_t k = i - 1; k <= i + 1; k++) {
					if ((k == i) && (l == j)) {
						continue;	// this is the pixel that we are considering itself, not the environment
					}
					
					if (thresholded_image->get(k, l) > 128) {
						neighbour_found = 1;
						break;
					}
				}
			}
			
			if (neighbour_found == 0) {
				// we didn't find an active point in the neighborhood, it was an isolated pixel
				processed_thresholded_image->set(i, j, 0);
			}
		}
	}
	
	return processed_thresholded_image;
}

ConvolveMatricesWithFFTClass::~ConvolveMatricesWithFFTClass() {
	if (forwardPlan != NULL) {
		fftw_destroy_plan(forwardPlan);
	}
	if (reversePlan != NULL) {
		fftw_destroy_plan(reversePlan);
	}
}

boost::shared_ptr<PALMMatrix<double> > ConvolveMatricesWithFFTClass::ConvolveMatricesWithFFT(boost::shared_ptr<PALMMatrix<double> > image1, boost::shared_ptr<PALMMatrix<double> > image2) {
	size_t x_size1, y_size1, x_size2, y_size2;
	
	x_size1 = image1->getXSize();
	y_size1 = image1->getYSize();
	x_size2 = image2->getXSize();
	y_size2 = image2->getYSize();
	
	size_t n_pixels, offset;
	size_t FFT_xSize, FFT_ySize, i, j;
	size_t n_FFT_values, nColumns;
	
	double *array1;
	double *array2;
	fftw_complex *array1_FFT;
	fftw_complex *array2_FFT;
	fftw_complex complex_value;
	boost::shared_ptr<PALMMatrix<double> > convolved_image(new PALMMatrix<double>(x_size1, y_size1));
	
	// are the dimensions equal?
	if ((x_size1 != x_size2) || (y_size1 != y_size2)) {
		XOPNotice("Error while convolving the images: the dimensions of the images are not equal.\r");
		throw INCOMPATIBLE_DIMENSIONING;
	}
	
	// does the image have dimension sizes that are odd? if so remove one column and/or row so that it becomes even
	// we will correct for this when returning the image by copying the values back in
	if ((x_size1 % 2) == 1) {	// odd x size
		FFT_xSize = x_size1 - 1;
	} else {	// even size, nothing needs to happen
		FFT_xSize = x_size1;
	}
	
	if ((y_size1 % 2) == 1) {	// odd y size
		FFT_ySize = y_size1 - 1;
	} else {	// even size, nothing needs to happen
		FFT_ySize = y_size1;
	}
	
	n_pixels = FFT_xSize * FFT_ySize;
	double normalization_factor = (double)(n_pixels);
	
	// convert both of the matrices to an array of doubles
	// we allocate the memory using fftw_malloc()
	// as this is recommended by the library
	array1 = (double *)fftw_malloc(sizeof(double) * n_pixels);
	if (array1 == NULL) {
		string error;
		error = "Error while convolving the images: unable to allocate array1\r";
		throw OUT_OF_MEMORY(error);
	}
	
	array2 = (double *)fftw_malloc(sizeof(double) * n_pixels);
	if (array2 == NULL) {
		fftw_free(array1);
		string error;
		error = "Error while convolving the images: unable to allocate array2\r";
		throw OUT_OF_MEMORY(error);
	}
	
	offset = 0;
	// now copy the matrices to the arrays
	
	for (size_t i = 0; i < FFT_xSize; i++) {
		for (size_t j = 0; j < FFT_ySize; j++) {
			// IMPORTANT: the data in the array is assumed to be in ROW-MAJOR order, so we loop over y first
			array1[offset] = image1->get(i, j);
			array2[offset] = image2->get(i, j);
			
			offset++;
		}
	}
	
	// now allocate the arrays that will hold the transformed result
	// the dimensions of these arrays are a bit unusual, and are x_size * (y_size / 2 + 1)
	n_FFT_values = FFT_xSize * (FFT_ySize / 2 + 1);
	nColumns = FFT_ySize / 2 + 1;
	
	array1_FFT = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * n_FFT_values);
	if (array1_FFT == NULL) {
		fftw_free(array1);
		fftw_free(array2);
		string error;
		error = "Error while convolving the images: unable to allocate array1_FFT\r";
		throw OUT_OF_MEMORY(error);
	}
	array2_FFT = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * n_FFT_values);
	if (array2_FFT == NULL) {
		fftw_free(array1);
		fftw_free(array2);
		fftw_free(array1_FFT);
		string error;
		error = "Error while convolving the images: unable to allocate array2_FFT\r";
		throw OUT_OF_MEMORY(error);
	}
	
	// prepare the transform and execute it on the first array
	// if there is no forward plan yet then create it
	forwardPlanMutex.lock();
	if ((forwardPlan == NULL) || (forwardPlanXSize != FFT_xSize) || (forwardPlanYSize != FFT_ySize)) {
		forwardCalculationMutex.lock();	// require exclusive ownership
		if (forwardPlan != NULL) {
			fftw_destroy_plan(forwardPlan);
		}
		
		forwardPlanXSize = FFT_xSize;
		forwardPlanYSize = FFT_ySize;
		
		forwardPlan = fftw_plan_dft_r2c_2d((int)(FFT_xSize), (int)(FFT_ySize), array1, array1_FFT, FFTW_ESTIMATE);
		forwardCalculationMutex.unlock();
	}
	
	forwardCalculationMutex.lock_shared();
	forwardPlanMutex.unlock();
	
	fftw_execute_dft_r2c(forwardPlan, array1, array1_FFT);
	
	// do the same on the second array
	fftw_execute_dft_r2c(forwardPlan, array2, array2_FFT);
	
	forwardCalculationMutex.unlock_shared();
	
	// now do the convolution
	for (size_t i = 0; i < n_FFT_values; i++) {
		complex_value[0] = array1_FFT[i][0] * array2_FFT[i][0] - array1_FFT[i][1] * array2_FFT[i][1];
		complex_value[1] = array1_FFT[i][0] * array2_FFT[i][1] + array1_FFT[i][1] * array2_FFT[i][0];
		
		// store the result in the first array
		// we add in a comb function so the origin is at the center of the image
		array1_FFT[i][0] = (2.0 * (double)((i / nColumns + i % nColumns) % 2) - 1.0) * -1.0 * complex_value[0];
		array1_FFT[i][1] = (2.0 * (double)((i / nColumns + i % nColumns) % 2) - 1.0) * -1.0 * complex_value[1];
	}
	
	// now do the reverse transform
	// we overwrite the original array
	// if there is no reverse plan yet then create it
	reversePlanMutex.lock();
	if ((reversePlan == NULL) || (reversePlanXSize != FFT_xSize) || (reversePlanYSize != FFT_ySize)) {
		reverseCalculationMutex.lock();	// require exclusive ownership
		if (reversePlan != NULL) {
			fftw_destroy_plan(reversePlan);
		}
		
		reversePlanXSize = FFT_xSize;
		reversePlanYSize = FFT_ySize;
		
		reversePlan = fftw_plan_dft_c2r_2d((int)(FFT_xSize), (int)(FFT_ySize), array1_FFT, array1, FFTW_ESTIMATE);
		reverseCalculationMutex.unlock();
	}
	
	reverseCalculationMutex.lock_shared();
	reversePlanMutex.unlock();
	
	fftw_execute_dft_c2r(reversePlan, array1_FFT, array1);
	
	reverseCalculationMutex.unlock_shared();
	
	// and store the result (we don't overwrite the input arguments)
	offset = 0;
	for (size_t i = 0; i < FFT_xSize; i++) {
		for (size_t j = 0; j < FFT_ySize; j++) {
			// the data in the array is assumed to be in ROW-MAJOR order, so we loop over x first
			// we also normalize the result
			convolved_image->set(i, j, (array1[offset] / normalization_factor));
			
			offset++;
		}
	}
	
	// if the number of rows was odd, make the last row (not included in the fft) a copy of that before it
	if ((x_size1 % 2) == 1) {
		i = x_size1 - 1;
		for (size_t j = 0; j < y_size1; ++j) {
			convolved_image->set(i, j, convolved_image->get(i - 1, j));
		}
	}
	
	// if the number of columns was odd, make the last column (not included in the fft) a copy of that before it
	if ((y_size1 % 2) == 1) {
		j = y_size1 - 1;
		for (size_t i = 0; i < x_size1; ++i) {
			convolved_image->set(i, j, convolved_image->get(i, j - 1));
		}
	}
	
	// if both the number of columns and the number of rows was odd, then the pixel at the top left (highest x, highest y) will be incorrect
	if (((x_size1 % 2) == 1) && ((y_size1 % 2) == 1)) {
		convolved_image->set(x_size1 - 1, y_size1 - 1, convolved_image->get(x_size1 - 2, y_size1 - 2));
	}
	
	// cleanup
	fftw_free(array1);
	fftw_free(array2);
	fftw_free(array1_FFT);
	fftw_free(array2_FFT);
	
	return convolved_image;
	
}

int Gauss_2D_fit_function(const gsl_vector *params, void *fitData_rhs, gsl_vector *deviations) {
	// params contains the current values of the parameters - amplitude, width, etc.
	// fitData is an object that contains the experimental data
	
	
	measured_data_Gauss_fits *fitDataLocal = (measured_data_Gauss_fits *)fitData_rhs;
	boost::shared_ptr<PALMMatrix<double> > imageSubset = fitDataLocal->imageSubset;
	
	size_t xSize = imageSubset->getXSize();
	size_t ySize = imageSubset->getXSize();
	size_t arrayOffset = 0;
	double xOffset = fitDataLocal->xOffset;
	double yOffset = fitDataLocal->yOffset;
	double sigma = fitDataLocal->sigma;
	
	double amplitude = gsl_vector_get(params, 0);
	double r = gsl_vector_get(params, 1);
	double x0 = gsl_vector_get(params, 2);
	double y0 = gsl_vector_get(params, 3);
	double offset = gsl_vector_get(params, 4);
	
	double x,y, function_value, square_deviation;
	
	if (r == 0) {
		return GSL_FAILURE;
	}
	
	for (size_t i = 0; i < xSize; ++i) {
		for (size_t j = 0; j < ySize; ++j) {
			x = xOffset + (double)i;
			y = yOffset + (double)j;
			function_value = offset + amplitude * exp(- (((x0 - x)/ r) * ((x0 - x)/ r) + ((y0 - y) / r) * ((y0 - y) / r)));
			square_deviation = (function_value - imageSubset->get(i, j)) / sigma;
			gsl_vector_set(deviations, arrayOffset, square_deviation);
			++arrayOffset;
		}
	}
	return GSL_SUCCESS;
}

int Gauss_2D_fit_function_Jacobian(const gsl_vector *params, void *fitData_rhs, gsl_matrix *jacobian) {
	measured_data_Gauss_fits *fitDataLocal = (measured_data_Gauss_fits *)fitData_rhs;
	boost::shared_ptr<PALMMatrix<double> > imageSubset = fitDataLocal->imageSubset;
	
	size_t xSize = imageSubset->getXSize();
	size_t ySize = imageSubset->getXSize();
	size_t arrayOffset = 0;
	double xOffset = fitDataLocal->xOffset;
	double yOffset = fitDataLocal->yOffset;
	double sigma = fitDataLocal->sigma;
	
	double amplitude = gsl_vector_get(params, 0);
	double r = gsl_vector_get(params, 1);
	double x0 = gsl_vector_get(params, 2);
	double y0 = gsl_vector_get(params, 3);
	// double offset = gsl_vector_get(params, 4);
	
	double x,y, exp_factor;
	double dfdA, dfdr, dfdx0, dfdy0,dfdoffset;
	
	if (r == 0) {
		return GSL_FAILURE;
	}
	
	for (size_t i = 0; i < xSize; ++i) {
		for (size_t j = 0; j < ySize; ++j) {
			x = xOffset + (double)i;
			y = yOffset + (double)j;
			
			exp_factor = exp(- ((x0 - x)/ r) * ((x0 - x)/ r) - ((y0 - y) / r) * ((y0 - y) / r));
			
			dfdA = exp_factor / sigma;
			dfdr = (2 * (y - y0) * (y - y0) / r / r / r + 2 * (x - x0) * (x - x0) / r / r / r) * exp_factor * amplitude / sigma;
			dfdx0 = (2 * (x - x0) * exp_factor * amplitude) / (r * r *sigma);
			dfdy0 = (2 * (y - y0) * exp_factor * amplitude) / (r * r *sigma);
			dfdoffset = 1/sigma;
			
			gsl_matrix_set(jacobian, arrayOffset, 0, dfdA);
			gsl_matrix_set(jacobian, arrayOffset, 1, dfdr);
			gsl_matrix_set(jacobian, arrayOffset, 2, dfdx0);
			gsl_matrix_set(jacobian, arrayOffset, 3, dfdy0);
			gsl_matrix_set(jacobian, arrayOffset, 4, dfdoffset);
			++arrayOffset;
		}
	}
	
	return GSL_SUCCESS;
	
}

int Gauss_2D_fit_function_and_Jacobian(const gsl_vector *params, void *measured_intensities_struct, gsl_vector *model_values, gsl_matrix *jacobian) {
	
	Gauss_2D_fit_function(params, measured_intensities_struct, model_values);
	Gauss_2D_fit_function_Jacobian(params, measured_intensities_struct, jacobian);
	
	return GSL_SUCCESS;
}

/*int Gauss_2D_Poissonian_fit_function(const gsl_vector *params, void *measured_intensities_struct, gsl_vector *model_values) {
	measured_data_Gauss_fits *intensities_local = (measured_data_Gauss_fits *)measured_intensities_struct;
	
	size_t number_of_intensities = intensities_local->get_number_of_intensities();
	double *measured_intensities = intensities_local->get_intensities();
	// double *sigma = intensities_local->get_sigma();
	size_t x_size = intensities_local->getXSize();
	size_t y_size = intensities_local->getYSize();
	double x_offset = intensities_local->get_x_offset();
	double y_offset = intensities_local->get_y_offset();
	
	double amplitude = gsl_vector_get(params, 0);
	double r = gsl_vector_get(params, 1);
	double x0 = gsl_vector_get(params, 2);
	double y0 = gsl_vector_get(params, 3);
	double offset = gsl_vector_get(params, 4);
	
	double x,y, function_value, square_deviation;
	
	if (r == 0) {
		return GSL_FAILURE;
	}
	
	for (size_t i = 0; i < number_of_intensities; i++) {
		return_xy_from_array(i, x_size, y_size, x_offset, y_offset, x, y);
		
		function_value = offset + amplitude * exp(- (((x0 - x)/ r) * ((x0 - x)/ r) + ((y0 - y) / r) * ((y0 - y) / r)));
		// square_deviation = (function_value - measured_intensities[i]) / sigma[i];
		square_deviation = sqrt(2 * measured_intensities[i] * log(measured_intensities[i] / function_value));
		gsl_vector_set(model_values, i, square_deviation);
	}
	return GSL_SUCCESS;
}

int Gauss_2D_Poissonian_fit_function_Jacobian(const gsl_vector *params, void *measured_intensities, gsl_matrix *jacobian) {
	
	measured_data_Gauss_fits *intensities_local = (measured_data_Gauss_fits *)measured_intensities;
	
	size_t number_of_intensities = intensities_local->get_number_of_intensities();
	// double *measured_intensities = intensities_local->get_intensities();
	// double *sigma = intensities_local->get_sigma();
	size_t x_size = intensities_local->getXSize();
	size_t y_size = intensities_local->getYSize();
	double x_offset = intensities_local->get_x_offset();
	double y_offset = intensities_local->get_y_offset();
	double *measured_intensities_array = intensities_local->get_intensities();
	
	double amplitude = gsl_vector_get(params, 0);
	double r = gsl_vector_get(params, 1);
	double x0 = gsl_vector_get(params, 2);
	double y0 = gsl_vector_get(params, 3);
	double offset = gsl_vector_get(params, 4);
	
	double x,y, exp_factor;
	double dfdA, dfdr, dfdx0, dfdy0,dfdoffset;
	
	if (r == 0) {
		return GSL_FAILURE;
	}
	
	// the maxima code to get the expression to derive from:
	// sqrt(2 * yi * log(yi / (offset + A * exp(-(((x0 - x) / r)^2 + ((y0 - y) / r)^2)))))
	
	double sqrt_2 = 1.414213562373095;
	double measured_intensity;
	double denominator;
	
	for (size_t i = 0; i < number_of_intensities; i++) {
		measured_intensity = measured_intensities_array[i];
		return_xy_from_array(i, x_size, y_size, x_offset, y_offset, x, y);
		
		exp_factor = exp(- ((x0 - x)/ r) * ((x0 - x)/ r) - ((y0 - y) / r) * ((y0 - y) / r));
		denominator = (exp_factor * amplitude  + offset) * sqrt(measured_intensity * log(measured_intensity / (exp_factor * amplitude + offset)));
		
		dfdA = - (sqrt_2 * exp_factor * measured_intensity) / (2.0 * denominator);
		dfdr = - sqrt_2 * (2.0 * (y - y0) * (y - y0) / r / r / r + 2 * (x - x0) * (x - x0) / r / r / r) * exp_factor * amplitude * measured_intensity / (2 * denominator);
		dfdx0 = (sqrt_2 * (x0 - x) * exp_factor * amplitude * measured_intensity) / (r * r * denominator);
		dfdy0 = (sqrt_2 * (y0 - y) * exp_factor * amplitude * measured_intensity) / (r * r * denominator);
		dfdoffset = (- sqrt_2 * measured_intensity) / (2.0 * denominator);
		
		gsl_matrix_set(jacobian, i, 0, dfdA);
		gsl_matrix_set(jacobian, i, 1, dfdr);
		gsl_matrix_set(jacobian, i, 2, dfdx0);
		gsl_matrix_set(jacobian, i, 3, dfdy0);
		gsl_matrix_set(jacobian, i, 4, dfdoffset);
	}
	
	return GSL_SUCCESS;
	
}

int Gauss_2D_Poissonian_fit_function_and_Jacobian(const gsl_vector *params, void *measured_intensities_struct, gsl_vector *model_values, gsl_matrix *jacobian) {
	
	Gauss_2D_Poissonian_fit_function(params, measured_intensities_struct, model_values);
	Gauss_2D_Poissonian_fit_function_Jacobian(params, measured_intensities_struct, jacobian);
	
	return GSL_SUCCESS;
} */


boost::shared_ptr<PALMMatrix<double> > FitPositionsGaussian::fit_positions(const boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<PALMMatrix<double> > positions) {
	
	size_t startPosition, endPosition;
	boost::shared_ptr<PALMMatrix<double> > fittedPositions;
	
	startPosition = 0;
	endPosition = positions->getXSize() - 1;
	
	fittedPositions = fit_positions(image, positions, startPosition, endPosition);
	
	return fittedPositions;
	
}

boost::shared_ptr<PALMMatrix<double> > FitPositionsGaussian::fit_positions(const boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<PALMMatrix<double> > positions, 
																		size_t startPos, size_t endPos) {
	
	// some safety checks
	if ((endPos >= positions->getXSize()) || (startPos >= positions->getXSize())) {
		string error;
		error = "Requested start and end positions are outside the range of positions supplied in FitPositionsGaussian::fit_positions";
		throw std::range_error(error);
	}
	
	if (startPos > endPos) {
		string error;
		error = "Start is beyond end in FitPositionsGaussian::fit_positions";
		throw std::range_error(error);
	}
	
	size_t number_of_positions = endPos - startPos + 1;
	size_t size_of_subset = 2 * cutoff_radius + 1;
	size_t x_offset, y_offset, x_max, y_max;
	size_t number_of_intensities = size_of_subset * size_of_subset;
	size_t xSize = image->getXSize();
	size_t ySize = image->getYSize();
	
	double x0_initial, y0_initial, initial_intensity, amplitude;
	double chi, degreesOfFreedom, c;
	long iterations = 0;
	int status;
	
	boost::shared_ptr<PALMMatrix<double> > image_subset;
	boost::shared_ptr<PALMMatrix<double> > fitted_positions;
	
	image_subset = boost::shared_ptr<PALMMatrix<double> > (new PALMMatrix<double>(size_of_subset, size_of_subset));
	fitted_positions = boost::shared_ptr<PALMMatrix<double> >(new PALMMatrix<double>(number_of_positions, N_OUTPUT_PARAMS_PER_FITTED_POSITION));
	
	fitted_positions->set_all(0);
	
	// initialize the solver
	const gsl_multifit_fdfsolver_type *solver;
	gsl_multifit_fdfsolver *fit_iterator;
	
	solver = gsl_multifit_fdfsolver_lmsder;
	measured_data_Gauss_fits fitData;
	gsl_multifit_function_fdf f;
	
	gsl_vector *fit_parameters = gsl_vector_alloc(5);
	if (fit_parameters == NULL) {
		string error;
		error = "unable to allocate fit_parameters in FitPositionsGaussian::fit_positions()\r";
		throw OUT_OF_MEMORY(error);
	}
	
	fit_iterator = gsl_multifit_fdfsolver_alloc(solver, number_of_intensities, 5);
	if (fit_iterator == NULL) {
		gsl_vector_free(fit_parameters);
		string error;
		error = "unable to allocate fit_iterator in FitPositionsGaussian::fit_positions()\r";
		throw OUT_OF_MEMORY(error);
	}
	
	gsl_matrix *covarianceMatrix = gsl_matrix_alloc(5, 5);
	if (covarianceMatrix == NULL) {
		gsl_vector_free(fit_parameters);
		gsl_multifit_fdfsolver_free(fit_iterator);
		string error;
		error = "unable to allocate covarianceMatrix in FitPositionsGaussian::fit_positions()\r";
		throw OUT_OF_MEMORY(error);
	}
	
	f.f = &Gauss_2D_fit_function;
	f.df = &Gauss_2D_fit_function_Jacobian;
	f.fdf = &Gauss_2D_fit_function_and_Jacobian;
	f.n = number_of_intensities;
	f.p = 5;
	f.params = (void *)&fitData;
	
	
	// iterate over all the determined positions
	for (size_t i = startPos; i <= endPos; i++) {
		iterations = 0;
		
		initial_intensity = positions->get(i, 0);
		x0_initial = positions->get(i, 1);
		y0_initial = positions->get(i, 2);
		
		amplitude = initial_intensity - background;
		
		x_offset = x0_initial - cutoff_radius;
		y_offset = y0_initial - cutoff_radius;
		x_max = x0_initial + cutoff_radius;
		y_max = y0_initial + cutoff_radius;
		
		if ((x_offset < 0) || (x_max > (xSize - 1)) || (y_offset < 0) || (y_max > (ySize - 1))) {	// this position is too close to the edge of the image
																									// we cannot include it
			continue;
		}
		
		for (size_t j = x_offset; j <= x_max; j++) {
			for (size_t k = y_offset; k <= y_max; k++) {
				image_subset->set(j - x_offset, k - y_offset, image->get(j, k));
			}
		}
		
		fitData.xOffset = (double)x_offset;
		fitData.yOffset = (double)y_offset;
		fitData.imageSubset = image_subset;
		fitData.sigma = sigma;
		
		// provide the initial parameters
		gsl_vector_set(fit_parameters, 0, amplitude);
		gsl_vector_set(fit_parameters, 1, r_initial * 1.414213562373095);	// because the fitting function is of the form 1/r^2, but standard deviation is 1/(2 r^2), we have to correct by sqrt(2)
		gsl_vector_set(fit_parameters, 2, x0_initial);
		gsl_vector_set(fit_parameters, 3, y0_initial);
		gsl_vector_set(fit_parameters, 4, background);
		
		// set the solver
		gsl_multifit_fdfsolver_set(fit_iterator, &f, fit_parameters);
		
		// run the iterations
		do {
			iterations++;
			status = gsl_multifit_fdfsolver_iterate(fit_iterator);
			if (status != 0)
				break;
			//			print_fit_state(iterations, fit_iterator);
			status = gsl_multifit_test_delta(fit_iterator->dx, fit_iterator->x, 10, 10);
		} while ((status = GSL_CONTINUE) && (iterations < 200));
		
		/*if ((status != GSL_SUCCESS) && (iterations == 200)) {
			// max number of iterations reached
		}
		if ((status != GSL_SUCCESS) && (iterations < 200)) {
			// some error occurred
			//			cout << gsl_strerror(status) << "\n";
		}*/
		
		// calculate the covariance matrix
		gsl_multifit_covar(fit_iterator->J, 0.0, covarianceMatrix);
		chi = gsl_blas_dnrm2(fit_iterator->f);
		degreesOfFreedom = (2 * cutoff_radius - 1) * (2 * cutoff_radius - 1) - 5;
		c = GSL_MAX_DBL(1, chi / sqrt(degreesOfFreedom));
		
		// store the data
		for (size_t j = 0; j < 5; ++j) {
			fitted_positions->set(i - startPos, j, gsl_vector_get(fit_iterator->x, j));
		}
		
		// store the errors
		for (size_t j = 5; j < 10; ++j) {
			fitted_positions->set(i - startPos, j, c * sqrt(gsl_matrix_get(covarianceMatrix, j - 5, j - 5)));
		}
		
		// store the number of iterations
		if ((status == GSL_SUCCESS) || (status == GSL_ETOLF) || (status == GSL_ETOLX)) {
			fitted_positions->set(i - startPos, 10, (double)iterations);
		} else {
			fitted_positions->set(i - startPos, 10, (double)(-1 * status));
		}
		
		// the width returned by the fit function is not equal to the standard deviation (a factor of sqrt 2 is missing)
		// so we correct for that
		
		fitted_positions->set(i - startPos, 1, fitted_positions->get(i - startPos, 1) / 1.414213562373095);
		fitted_positions->set(i - startPos, 6, fitted_positions->get(i - startPos, 6) / 1.414213562373095);	// the same for the error
	}
	
	gsl_multifit_fdfsolver_free(fit_iterator);
	gsl_vector_free(fit_parameters);
	
	return fitted_positions;
	
}


boost::shared_ptr<PALMMatrix<double> > FitPositionsMultiplication::fit_positions(const boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<PALMMatrix<double> > positions) {
	size_t startPosition, endPosition;
	boost::shared_ptr<PALMMatrix<double> > fittedPositions;
	
	startPosition = 0;
	endPosition = positions->getXSize() - 1;
	
	fittedPositions = fit_positions(image, positions, startPosition, endPosition);
	
	return fittedPositions;
}
		

boost::shared_ptr<PALMMatrix<double> > FitPositionsMultiplication::fit_positions(const boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<PALMMatrix<double> > positions, 
																			  size_t startPos, size_t endPos) {
	
	// some safety checks
	if ((endPos >= positions->getXSize()) || (startPos >= positions->getXSize())) {
		string error;
		error = "Requested start and end positions are outside the range of positions supplied in FitPositionsGaussian::fit_positions";
		throw std::range_error(error);
	}
	
	if (startPos > endPos) {
		string error;
		error = "Start is beyond end in FitPositionsGaussian::fit_positions";
		throw std::range_error(error);
	}
	
	size_t number_of_positions = endPos - startPos + 1;
	size_t size_of_subset = 2 * cutoff_radius + 1;
	size_t xSize = image->getXSize();
	size_t ySize = image->getYSize();
	size_t x_offset, y_offset, x_max, y_max;
	
	double x0_initial, y0_initial, initial_intensity, amplitude;
	size_t iterations = 0;
	
	double convergence_treshold_squared = convergence_threshold * convergence_threshold;
	double delta_squared = 10 * convergence_treshold_squared;	// this test the convergence of the position determined by the iteration
	// it is the distance between (xn-1, yn-1) and (xn, yn)
	// we initialize it to a value well over the treshold so that we will run at least two iterations
	double previous_position_x, previous_position_y;
	double current_x, current_y;
	
	boost::shared_ptr<PALMMatrix<double> > image_subset;
	boost::shared_ptr<PALMMatrix<double> > image_subset_mask;
	boost::shared_ptr<PALMMatrix<double> > fitted_positions;
	
	image_subset = boost::shared_ptr<PALMMatrix<double> >(new PALMMatrix<double>(size_of_subset, size_of_subset));
	image_subset_mask = boost::shared_ptr<PALMMatrix<double> > (new PALMMatrix<double>(size_of_subset, size_of_subset));
	fitted_positions = boost::shared_ptr<PALMMatrix<double> >(new PALMMatrix<double>(number_of_positions, N_OUTPUT_PARAMS_PER_FITTED_POSITION));
	
	fitted_positions->set_all(0);
	
	for (size_t i = startPos; i <= endPos; ++i) {
		initial_intensity = positions->get(i, 0);
		x0_initial = positions->get(i, 1);
		y0_initial = positions->get(i, 2);
		
		amplitude = initial_intensity - background;
		
		x_offset = (size_t)x0_initial - cutoff_radius;
		y_offset = (size_t)y0_initial - cutoff_radius;
		x_max = (size_t)x0_initial + cutoff_radius;
		y_max = (size_t)y0_initial + cutoff_radius;
		
		if ((x_offset < 0) || (x_max > (xSize - 1)) || (y_offset < 0) || (y_max > (ySize - 1))) {	// this position is too close to the edge of the image, we cannot include it
			continue;
		}
		
		for (size_t k = y_offset; k <= y_max; ++k) {
			for (size_t j = x_offset; j <= x_max; ++j) {
				image_subset->set(j - x_offset, k - y_offset, image->get(j, k));
			}
		}
		
		iterations = 0;
		
		current_x = x0_initial - (double)x_offset;	// correct the x- and y-values for the fact that we analyze in a subset of the image rather than the complete frame
		current_y = y0_initial - (double)y_offset;
		
		while (delta_squared > convergence_treshold_squared) {
			previous_position_x = current_x;
			previous_position_y = current_y;
			
			++iterations;
			
			if (iterations > 100) {	// the multiplication is not converging, we should stop
				SetNaN64(&current_x);
				SetNaN64(&current_y);
				break;
			}
			
			multiply_with_gaussian(image_subset, image_subset_mask, current_x, current_y, r_initial, background, amplitude);
			determine_x_y_position(image_subset_mask, current_x, current_y);
			
			if (iterations == 1)	// this is the first iteration, we should not check for termination
				continue;
			
			delta_squared = (current_x - previous_position_x) * (current_x - previous_position_x) + (current_y - previous_position_y) * (current_y - previous_position_y);
		}
		
		delta_squared = 10 * convergence_treshold_squared;
		
		fitted_positions->set(i - startPos, 2, current_x + (double)x_offset);
		fitted_positions->set(i - startPos, 3, current_y + (double)y_offset);
		fitted_positions->set(i - startPos, 10, (double)iterations);
	}
	
	return fitted_positions;
}
		

int FitPositionsMultiplication::multiply_with_gaussian(boost::shared_ptr<PALMMatrix<double> > original_image, boost::shared_ptr<PALMMatrix<double> > masked_image, double x, double y, 
													   double std_dev, double background, double amplitude) {
	// we will replace the contents of masked_image with the multiplication of original_image and a gaussian centered at position (x,y)
	
	size_t x_size = masked_image->getXSize();
	size_t y_size = masked_image->getYSize();
	
	double gaussian_value, distance_squared;
	
	if ((original_image->getXSize() != x_size) || (original_image->getYSize() != y_size)) {
		throw DIMENSIONS_SHOULD_BE_EQUAL();
	}
	
	for (size_t j = 0; j < y_size; j++) {
		for (size_t i = 0; i < x_size; i++) {
			distance_squared = (x - (double)i) * (x - (double)i) + (y - (double)j) * (y - (double)j);
			
			gaussian_value = amplitude * exp(- distance_squared / (2 * std_dev * std_dev)) + background;
			
			masked_image->set(i, j, gaussian_value * original_image->get(i, j));
		}
	}
	
	return 0;
}


int FitPositionsMultiplication::determine_x_y_position(boost::shared_ptr<PALMMatrix<double> > masked_image, double &x, double &y) {
	// based on eq (3) in Thompson Biophys J 2002
	
	size_t x_size = (size_t)masked_image->getXSize();
	size_t y_size = (size_t)masked_image->getYSize();
	
	double numerator_x = 0, denominator = 0;
	double numerator_y = 0;
	
	// start with determining the x-position
	for (size_t j = 0; j < y_size; j++) {
		for (size_t i = 0; i < x_size; i++) {
			numerator_x += (double)i * masked_image->get(i, j);
			numerator_y += (double)j * masked_image->get(i, j);
			denominator += masked_image->get(i, j);
		}
	}
	
	x = numerator_x / denominator;
	
	y = numerator_y / denominator;
	
	return 0;
}


boost::shared_ptr<PALMMatrix<double> > FitPositionsCentroid::fit_positions(const boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<PALMMatrix<double> > positions) {
	size_t startPosition, endPosition;
	boost::shared_ptr<PALMMatrix<double> > fittedPositions;
	
	startPosition = 0;
	endPosition = positions->getXSize() - 1;
	
	fittedPositions = fit_positions(image, positions, startPosition, endPosition);
	
	return fittedPositions;
}


boost::shared_ptr<PALMMatrix<double> > FitPositionsCentroid::fit_positions(const boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<PALMMatrix<double> > positions, 
																			  size_t startPos, size_t endPos) {
	
	// some safety checks
	if ((endPos >= positions->getXSize()) || (startPos >= positions->getXSize())) {
		string error;
		error = "Requested start and end positions are outside the range of positions supplied in FitPositionsCentroid::fit_positions";
		throw std::range_error(error);
	}
	
	if (startPos > endPos) {
		string error;
		error = "Start is beyond end in FitPositionsCentroid::fit_positions";
		throw std::range_error(error);
	}
	
	size_t number_of_positions = endPos - startPos + 1;
	size_t xSize = image->getXSize();
	size_t ySize = image->getYSize();
	size_t x_offset, y_offset, x_max, y_max;
	
	size_t x0_initial, y0_initial;
	double current_x, current_y;
	double denominator;
	
	boost::shared_ptr<PALMMatrix<double> > fitted_positions;
	
	fitted_positions = boost::shared_ptr<PALMMatrix<double> >(new PALMMatrix<double>(number_of_positions, N_OUTPUT_PARAMS_PER_FITTED_POSITION));
	
	fitted_positions->set_all(0);
	
	for (size_t i = startPos; i <= endPos; ++i) {
		x0_initial = positions->get(i, 1);
		y0_initial = positions->get(i, 2);
		current_x = 0;
		current_y = 0;
		denominator = 0;
		
		x_offset = x0_initial - cutoff_radius;
		y_offset = y0_initial - cutoff_radius;
		x_max = x0_initial + cutoff_radius;
		y_max = y0_initial + cutoff_radius;
		
		if ((x_offset > xSize) || (x_max > (xSize - 1)) || (y_offset > ySize) || (y_max > (ySize - 1))) {	// the point is too close to the edge
																											// because all the variables are unsigned, a negative value will
																											// actually end up being larger than xSize or ySize
			continue;
		}
		
		for (size_t j = x_offset; j <= x_max; ++j) {
			for (size_t k = y_offset; k <= y_max; ++k) {
				current_x += (double)j * image->get(j, k);
				current_y += (double)k * image->get(j, k);
				denominator += image->get(j, k);
			}
		}
		
		current_x /= denominator;
		current_y /= denominator;		
		
		
		fitted_positions->set(i - startPos, 2, current_x);
		fitted_positions->set(i - startPos, 3, current_y);
	}
	
	return fitted_positions;
}