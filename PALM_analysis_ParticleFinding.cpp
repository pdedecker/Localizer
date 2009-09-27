/*
 *  PALM_analysis_ParticleFinding.cpp
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 27/09/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "PALM_analysis_ParticleFinding.h"

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
	double backgroundIntensity = 0;
	size_t nBackgroundPixels = 0;
	
	// we run over all the points in the image to see if they are above the treshold
	for (size_t j = minDistanceFromEdge; j < y_size - minDistanceFromEdge; j++) {
		for (size_t i = minDistanceFromEdge; i < x_size - minDistanceFromEdge; i++) {
			
			if (threshold_image->get(i, j) < 128) { // we don't care about this point, it's not included in the thresholded image
				backgroundIntensity += (*image)(i, j);	// but use it to estimate the background intensity
				++nBackgroundPixels;
				continue;
			}
			
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
	
	output_positions = boost::shared_ptr<PALMMatrix<double> >(new PALMMatrix<double>(number_of_positions, 4));
	backgroundIntensity /= (double)nBackgroundPixels;
	
	for (size_t i = 0; i < number_of_positions; i++) {
		(*output_positions)(i, 0) = positions[i].get_intensity() - backgroundIntensity;
		(*output_positions)(i, 1) = positions[i].get_x();
		(*output_positions)(i, 2) = positions[i].get_y();
		(*output_positions)(i, 3) = backgroundIntensity;
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
	double backgroundIntensity = 0;
	size_t nBackgroundPixels = 0;
	
	mapped_image = boost::shared_ptr<PALMMatrix<long> >(new PALMMatrix<long>(x_size, y_size));
	
	mapped_image->set_all(-1);
	
	for (size_t j = minDistanceFromEdge; j < y_size - minDistanceFromEdge; ++j) {
		for (size_t i = minDistanceFromEdge; i < x_size - minDistanceFromEdge; ++i) {	// loop over the entire image
			
			if ((*threshold_image)(i, j) < 128) { // we don't care about this point, it's not included in the thresholded image
				backgroundIntensity += (*image)(i, j);	// but use it to estimate the background intensity
				++nBackgroundPixels;
				continue;
			}
			
			if ((*mapped_image)(i, j) != -1) {	// this point is already assigned to another particle
				continue;
			}
			
			// if we are still here then we found an active pixel that is not yet assigned to a new particle
			// we create a new particle at this position
			positionsInCurrentParticle.clear();
			positionsInCurrentParticleList.clear();
			
			(*mapped_image)(i, j) = particleIndex;
			
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
	
	output_positions = boost::shared_ptr<PALMMatrix<double> >(new PALMMatrix<double>(particles.size(), 4));
	backgroundIntensity /= (double)nBackgroundPixels;
	
	// now copy the output to a gsl matrix
	for (size_t k = 0; k < particles.size(); ++k) {
		(*output_positions)(k, 0) = particles[k].get_intensity() - backgroundIntensity;
		(*output_positions)(k, 1) = particles[k].get_x();
		(*output_positions)(k, 2) = particles[k].get_y();
		(*output_positions)(k, 3) = backgroundIntensity;
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
	
	long particleIndex = (*mapped_image)(x, y);
	
	if ((x < minDistanceFromEdge) || (x > x_size - minDistanceFromEdge - 1))
		return;
	if ((y < minDistanceFromEdge) || (y > y_size - minDistanceFromEdge - 1))
		return;
	
	// is the pixel to the left of the current one active?
	if (x > 0) {
		if ((*threshold_image)(x - 1, y) > 128) {
			// it's active
			// did we already include this pixel?
			if ((*mapped_image)(x - 1, y) == -1) {
				(*mapped_image)(x - 1, y) = particleIndex;
				// add the point to the vector
				currentPosition.set_x((double)x - 1);
				currentPosition.set_y((double)y);
				currentPosition.set_intensity(image->get(x - 1, y));
				positionsInCurrentParticle.push_back(currentPosition);
			}
		}
	}
	// is the pixel to the right of the current one active?
	if (x < x_size - 1) {
		if ((*threshold_image)(x + 1, y) > 128) {
			// it's active
			// did we already include this pixel?
			if ((*mapped_image)(x + 1, y) == -1) {
				(*mapped_image)(x + 1, y) = particleIndex;
				// add the point to the vector
				currentPosition.set_x((double)x + 1);
				currentPosition.set_y((double)y);
				currentPosition.set_intensity(image->get(x + 1, y));
				positionsInCurrentParticle.push_back(currentPosition);
			}
		}
	}
	// is the pixel above the current one active?
	if (y > 0) {
		if ((*threshold_image)(x, y - 1) > 128) {
			// it's active
			// did we already include this pixel?
			if ((*mapped_image)(x, y - 1) == -1) {
				(*mapped_image)(x, y - 1) = particleIndex;
				// add the point to the vector
				currentPosition.set_x((double)x);
				currentPosition.set_y((double)y - 1);
				currentPosition.set_intensity(image->get(x, y - 1));
				positionsInCurrentParticle.push_back(currentPosition);
			}
		}
	}
	// is the pixel below the current one active?
	if (y < y_size - 1) {
		if ((*threshold_image)(x, y + 1) > 128) {
			// it's active
			// did we already include this pixel?
			if ((*mapped_image)(x, y + 1) == -1) {
				(*mapped_image)(x, y + 1) = particleIndex;
				// add the point to the vector
				currentPosition.set_x((double)x);
				currentPosition.set_y((double)y + 1);
				currentPosition.set_intensity(image->get(x, y + 1));
				positionsInCurrentParticle.push_back(currentPosition);
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
	double backgroundIntensity = 0;
	size_t nBackgroundPixels = 0;
	
	mapped_image = boost::shared_ptr<PALMMatrix<long> >(new PALMMatrix<long>(x_size, y_size));
	
	mapped_image->set_all(-1);
	
	for (size_t j = minDistanceFromEdge; j < y_size - minDistanceFromEdge; ++j) {
		for (size_t i = minDistanceFromEdge; i < x_size - minDistanceFromEdge; ++i) {	// loop over the entire image
			
			if (threshold_image->get(i, j) < 128) { // we don't care about this point, it's not included in the thresholded image
				backgroundIntensity += (*image)(i, j);	// but use it to estimate the background intensity
				++nBackgroundPixels;
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
	
	output_positions = boost::shared_ptr<PALMMatrix<double> >(new PALMMatrix<double>(particles.size(), 4));
	backgroundIntensity /= (double)nBackgroundPixels;
	
	// now copy the output to a gsl matrix
	for (size_t k = 0; k < particles.size(); ++k) {
		output_positions->set(k, 0, particles[k].get_intensity() - backgroundIntensity);
		output_positions->set(k, 1, particles[k].get_x());
		output_positions->set(k, 2, particles[k].get_y());
		(*output_positions)(k, 3) = backgroundIntensity;
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