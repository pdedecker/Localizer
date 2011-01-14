/*
 *  PALM_analysis_ParticleFinding.cpp
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 27/09/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "PALM_analysis_ParticleFinding.h"

boost::shared_ptr<std::list<position> > ParticleFinder_radius::findPositions(boost::shared_ptr<Eigen::MatrixXd> image, boost::shared_ptr<Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> > threshold_image) {
	boost::shared_ptr<std::list<position> > positions (new std::list<position>());
	position currentPosition;
	// we store the pixels above the treshold as a vector containing x,y,intensity
	size_t x_size = image->rows(), y_size = image->cols();
	size_t number_of_positions;
	double current_intensity, previous_intensity, current_x, current_y;
	double distance_squared;
	double radius_squared = (double)(radius * radius);
	double x, y;
	boost::shared_ptr<Eigen::MatrixXd> output_positions;
	double backgroundIntensity = 0;
	size_t nBackgroundPixels = 0;
	
	// we run over all the points in the image to see if they are above the treshold
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {
			
			if ((*threshold_image)(i, j) < 128) { // we don't care about this point, it's not included in the thresholded image
				backgroundIntensity += (*image)(i, j);	// but use it to estimate the background intensity
				++nBackgroundPixels;
				continue;
			}
			
			current_intensity = (*image)(i, j);
			
			current_x = (double)i;
			current_y = (double)j;
			
			// if we are still here then we need to take a closer look at this point
			// check if the current point overlaps with a previous point
			number_of_positions = positions->size();
			for (std::list<position>::iterator it = positions->begin(); it != positions->end(); ++it) {
				x = (*it).get_x();
				y = (*it).get_y();
				distance_squared = (current_x - x) * (current_x - x) + (current_y - y) * (current_y - y);
				
				if (distance_squared < radius_squared) {
					// we need to skip one of the two pixels that we are comparing
					// we will keep the pixel with the largest intensity
					previous_intensity = (*it).get_intensity();
					if (current_intensity > previous_intensity) {
						(*it).set_intensity(current_intensity);
						(*it).set_x(current_x);
						(*it).set_y(current_y);
					}
					continue;
				}
			}
			
			// this position is a new emitter
			currentPosition.set_intensity(current_intensity);
			currentPosition.set_x(current_x);
			currentPosition.set_y(current_y);
			positions->push_back(currentPosition);
		}
	}
	
	backgroundIntensity /= (double)nBackgroundPixels;
	
	// update the amplitudes/intensities on all of the positions
	for (std::list<position>::iterator it = positions->begin(); it != positions->end(); ++it) {
		(*it).set_intensity((*it).get_intensity() - backgroundIntensity);
		(*it).set_background(backgroundIntensity);
	}
	
	
	return positions;
}


boost::shared_ptr<std::list<position> > ParticleFinder_adjacent4::findPositions(boost::shared_ptr<Eigen::MatrixXd> image, boost::shared_ptr<Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> > threshold_image) {
	
	boost::shared_ptr<std::list<position> > particles (new std::list<position>());
	std::list<position> positionsInCurrentParticleList;
	std::vector<position> positionsInCurrentParticle;
	position currentPosition;
	boost::shared_ptr<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > mapped_image;	// keeps track of which pixels have already been mapped to a particle
	size_t x_size = image->rows();
	size_t y_size = image->cols();
	size_t x, y;
	long particleIndex = 0;
	double average_x, average_y;
	double maxIntensity;
	double backgroundIntensity = 0;
	size_t nBackgroundPixels = 0;
	
	mapped_image = boost::shared_ptr<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >(new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>((int)x_size, (int)y_size));
	mapped_image->setConstant(-1);
	
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {	// loop over the entire image
			
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
			currentPosition.set_intensity((*image)(i, j));
			
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
			
			particles->push_back(currentPosition);
		}
	}
	
	backgroundIntensity /= (double)nBackgroundPixels;
	
	// update the amplitudes/intensities on all of the positions
	for (std::list<position>::iterator it = particles->begin(); it != particles->end(); ++it) {
		(*it).set_intensity((*it).get_intensity() - backgroundIntensity);
		(*it).set_background(backgroundIntensity);
	}
	
	return particles;
	
}

void ParticleFinder_adjacent4::growParticle(position centerPosition, std::list<position> &positionsInCurrentParticle, boost::shared_ptr<Eigen::MatrixXd> image, boost::shared_ptr<Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> > threshold_image, boost::shared_ptr<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > mapped_image) {
	// the pixel at position (x,y) belongs to a particle
	// do the surrounding pixels belong to the same particle?
	
	// the function checks which of the pixels surrounding pos are active
	// if one or more of these is active, then it checks if they are already assigned to the particle by checking mapped_image
	// if they are not known then they are added to to the list with positions of the current particle
	// and also added to mapped_image
	
	size_t x_size = image->rows();
	size_t y_size = image->cols();
	position currentPosition;
	
	size_t x = (size_t)(centerPosition.get_x() + 0.5);
	size_t y = (size_t)(centerPosition.get_y() + 0.5);
	
	long particleIndex = (*mapped_image)(x, y);
	
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
				currentPosition.set_intensity((*image)(x - 1, y));
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
				currentPosition.set_intensity((*image)(x + 1, y));
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
				currentPosition.set_intensity((*image)(x, y - 1));
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
				currentPosition.set_intensity((*image)(x, y + 1));
				positionsInCurrentParticle.push_back(currentPosition);
			}
		}
	}
}


boost::shared_ptr<std::list<position> > ParticleFinder_adjacent8::findPositions(boost::shared_ptr<Eigen::MatrixXd> image, boost::shared_ptr<Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> > threshold_image) {
	
	std::list<position> positionsInCurrentParticleList;
	std::vector<position> positionsInCurrentParticle;
	position currentPosition;
	boost::shared_ptr<std::list<position> > particles(new std::list<position>);
	boost::shared_ptr<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > mapped_image;	// keeps track of which pixels have already been mapped to a particle
	size_t x_size = image->rows();
	size_t y_size = image->cols();
	size_t x, y;
	long particleIndex = 0;
	double average_x, average_y;
	double maxIntensity;
	double backgroundIntensity = 0;
	size_t nBackgroundPixels = 0;
	
	mapped_image = boost::shared_ptr<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >(new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>((int)x_size, (int)y_size));
	mapped_image->setConstant(-1);
	
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j = 0; j < y_size; j++) {	// loop over the entire image
			
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
			currentPosition.set_intensity((*image)(i, j));
			
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
			
			particles->push_back(currentPosition);
		}
	}
	
	backgroundIntensity /= (double)nBackgroundPixels;
	
	// update the amplitudes/intensities on all of the positions
	for (std::list<position>::iterator it = particles->begin(); it != particles->end(); ++it) {
		(*it).set_intensity((*it).get_intensity() - backgroundIntensity);
		(*it).set_background(backgroundIntensity);
	}
	
	return particles;
	
}

void ParticleFinder_adjacent8::growParticle(position centerPosition, std::list<position> &positionsInCurrentParticle, boost::shared_ptr<Eigen::MatrixXd> image, boost::shared_ptr<Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> > threshold_image, boost::shared_ptr<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > mapped_image) {
	// the pixel at position (x,y) belongs to a particle
	// do the surrounding pixels belong to the same particle?
	
	// the function checks which of the pixels surrounding pos are active
	// if one or more of these is active, then it checks if they are already assigned to the particle by checking mapped_image
	// if they are not known then they are added to to the list with positions of the current particle
	// and also added to mapped_image
	
	size_t x_size = image->rows();
	size_t y_size = image->cols();
	position currentPosition;
	
	size_t x = (size_t)(centerPosition.get_x() + 0.5);
	size_t y = (size_t)(centerPosition.get_y() + 0.5);
	
	// make sure that we don't cross the boundaries of the image
	size_t lowerXBound = (x - 1) == (size_t)-1 ? 0 : x - 1;
	size_t upperXBound = (x + 1) >= x_size ? x_size - 1 : x + 1;
	size_t lowerYBound = (y - 1) == (size_t)-1 ? 0 : y - 1;
	size_t upperYBound = (y + 1) >= y_size ? y_size - 1 : y + 1;
	
	long particleIndex = (*mapped_image)(x, y);
	
	for (size_t i = lowerXBound; i <= upperXBound; ++i) {
		for (size_t j = lowerYBound; j <= upperYBound; ++j) {
			
			if ((i == x) && (j == y))
				continue;
			
			if ((*threshold_image)(i, j) < 128)
				continue;
			
			if ((*mapped_image)(i, j) != -1)
				continue;
			
			// add the current position
			(*mapped_image)(i, j) = particleIndex;
			currentPosition.set_x((double)i);
			currentPosition.set_y((double)j);
			currentPosition.set_intensity((*image)(i, j));
			positionsInCurrentParticle.push_back(currentPosition);
		}
	}
}

void ParticleVerifier_RemoveOverlappingParticles::VerifyParticles(boost::shared_ptr<Eigen::MatrixXd> image, boost::shared_ptr<std::list<position> > positions) {
	double distance;
	double minDistance = 4.0 * this->psfWidth;
	
	for (std::list<position>::iterator it1 = positions->begin(); it1 != positions->end(); ++it1) {
		// since the list::iterator does not support operator+() we have to be kludgy
		for (std::list<position>::iterator it2 = (++it1)--; it2 != positions->end(); ++it2) {
			distance = sqrt(((*it1).get_x() - (*it2).get_x()) * ((*it1).get_x() - (*it2).get_x()) + ((*it1).get_y() - (*it2).get_y()) * ((*it1).get_y() - (*it2).get_y()));
			if (distance < minDistance) {
				// these two points should be deleted
				// don't delete them now since there might be a third point overlapping with one of them
				// instead just mark them as needing deletion by changing their amplitude to some
				// unlikely value
				(*it1).set_intensity(-1.0e200);
				(*it2).set_intensity(-1.0e200);
			}
		}
	}
	
	// now delete the positions that failed the test
	for (std::list<position>::iterator it = positions->begin(); it != positions->end(); ++it) {
		if ((*it).get_intensity() == -1.0e200) {
			it = positions->erase(it);
			if (it != positions->begin())
				--it;
		}
	}
	
}
