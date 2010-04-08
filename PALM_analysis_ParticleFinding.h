/*
 *  PALM_analysis_ParticleFinding.h
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 27/09/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef PALM_ANALYSIS_PARTICLEFINDING_H
#define PALM_ANALYSIS_PARTICLEFINDING_H

#include "PALM_analysis_storage.h"
#include "boost/smart_ptr.hpp"
#include <list>


/**
 * @brief An abstract base class from which all other 'particle finder' classes must derive
 * 
 * This family of classes contains algorithms that accept a single binary segment images as input and reduce this to a list of (x, y) positions
 * where a localization fit will be attempted.
 *
 * This class is abstract, meaning that it can not be instantiated. However, it provides the virtual function 'findPositions()', 
 * which is where derived classes should do their work. Any necessary parameters should be passed in the constructors of derived classes.
 */
class ParticleFinder {
public:
	ParticleFinder() {;}
	virtual ~ParticleFinder() {;}
	
	virtual boost::shared_ptr<std::vector<position> > findPositions(boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<PALMMatrix <unsigned char> > threshold_image) = 0;
protected:
};

/**
 * @brief Assumes that all active pixels within a given radius belong to the same particle
 */
class ParticleFinder_radius : public ParticleFinder {
public:
	ParticleFinder_radius(double radius_rhs) {radius = radius_rhs;}
	~ParticleFinder_radius() {;}
	
	boost::shared_ptr<std::vector<position> > findPositions(boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<PALMMatrix <unsigned char> > threshold_image);
protected:
	double radius;
};

/**
 * @brief Assumes that all neighbouring active pixels (four-way adjacency) belong to the same particle
 */
class ParticleFinder_adjacent4 : public ParticleFinder {
public:
	ParticleFinder_adjacent4() {;}
	~ParticleFinder_adjacent4() {;}
	
	boost::shared_ptr<std::vector<position> > findPositions(boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<PALMMatrix <unsigned char> > threshold_image);
	
protected:
	void growParticle(position centerPosition, std::list<position> &positionsInCurrentParticle, boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<PALMMatrix <unsigned char> > threshold_image, boost::shared_ptr<PALMMatrix<long> > mapped_image);
	
};

/**
 * @brief Assumes that all neighbouring active pixels (eight-way adjacency) belong to the same particle
 */
class ParticleFinder_adjacent8 : public ParticleFinder {
public:
	ParticleFinder_adjacent8() {;}
	~ParticleFinder_adjacent8() {;}
	
	boost::shared_ptr<std::vector<position> > findPositions(boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<PALMMatrix <unsigned char> > threshold_image);
	
protected:
	void growParticle(position centerPosition, std::list<position> &positionsInCurrentParticle, boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<PALMMatrix <unsigned char> > threshold_image, boost::shared_ptr<PALMMatrix<long> > mapped_image);
	
};

#endif

