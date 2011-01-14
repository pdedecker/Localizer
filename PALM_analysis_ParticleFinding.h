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

#include <algorithm>
#include "PALM_analysis_storage.h"
#include "PALM_analysis_Localization.h"
#include "boost/smart_ptr.hpp"
#define EIGEN_DEFAULT_TO_ROW_MAJOR
#include <Eigen/Eigen>
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
	
	virtual boost::shared_ptr<std::list<position> > findPositions(boost::shared_ptr<Eigen::MatrixXd> image, boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > threshold_image) = 0;
protected:
};

/**
 * @brief Assumes that all active pixels within a given radius belong to the same particle
 */
class ParticleFinder_radius : public ParticleFinder {
public:
	ParticleFinder_radius(double radius_rhs) {radius = radius_rhs;}
	~ParticleFinder_radius() {;}
	
	boost::shared_ptr<std::list<position> > findPositions(boost::shared_ptr<Eigen::MatrixXd> image, boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > threshold_image);
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
	
	boost::shared_ptr<std::list<position> > findPositions(boost::shared_ptr<Eigen::MatrixXd> image, boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > threshold_image);
	
protected:
	void growParticle(position centerPosition, std::list<position> &positionsInCurrentParticle, boost::shared_ptr<Eigen::MatrixXd> image, boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > threshold_image, boost::shared_ptr<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > mapped_image);
	
};

/**
 * @brief Assumes that all neighbouring active pixels (eight-way adjacency) belong to the same particle
 */
class ParticleFinder_adjacent8 : public ParticleFinder {
public:
	ParticleFinder_adjacent8() {;}
	~ParticleFinder_adjacent8() {;}
	
	boost::shared_ptr<std::list<position> > findPositions(boost::shared_ptr<Eigen::MatrixXd> image, boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > threshold_image);
	
protected:
	void growParticle(position centerPosition, std::list<position> &positionsInCurrentParticle, boost::shared_ptr<Eigen::MatrixXd> image, boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > threshold_image, boost::shared_ptr<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > mapped_image);
	
};

/**
 * @brief An abstract base class from which all other 'particle verifer' classes must derive
 * 
 * This family of classes accepts lists of particles, derived from the segmented image, and the image itself and attempts to verify if these
 * positions are eligible as single emitters.
 */
class ParticleVerifier {
public:
	ParticleVerifier() {;}
	virtual ~ParticleVerifier() {;}
	
	/**
	 * Verify the positions. Does not return anything, but rather directly deletes elements from positions if not eligible
	 */
	virtual void VerifyParticles(boost::shared_ptr<Eigen::MatrixXd> image, boost::shared_ptr<std::list<position> > positions) = 0;
protected:
};

/**
 * Verify the particles by eliminating any two that are within 4 times the PSF deviation
 */
class ParticleVerifier_RemoveOverlappingParticles : public ParticleVerifier {
public:
	ParticleVerifier_RemoveOverlappingParticles(double psfWidth_rhs) : psfWidth(psfWidth_rhs) {}
	~ParticleVerifier_RemoveOverlappingParticles() {;}
	
	
	void VerifyParticles(boost::shared_ptr<Eigen::MatrixXd> image, boost::shared_ptr<std::list<position> > positions);
protected:
	double psfWidth;
};


/**
 * Verify the particles by seeing if they pass the fitting of a FitPositions_SymmetricGaussian class
 */
class ParticleVerifier_SymmetricGaussian : public ParticleVerifier {
public:
	ParticleVerifier_SymmetricGaussian(double initialPSFWidth_rhs, double sigma_rhs) : verifier(FitPositions_SymmetricGaussian(initialPSFWidth_rhs, sigma_rhs)) {}
	~ParticleVerifier_SymmetricGaussian() {;}
	
	
	void VerifyParticles(boost::shared_ptr<Eigen::MatrixXd> image, boost::shared_ptr<std::list<position> > positions) {
		verifier.fit_positions(image, positions);
	}
protected:
	FitPositions_SymmetricGaussian verifier;
};

/**
 * Verify the particles by seeing if they pass the fitting of a FitPositions_EllipsoidalGaussian class
 */
class ParticleVerifier_EllipsoidalGaussian : public ParticleVerifier {
public:
	ParticleVerifier_EllipsoidalGaussian(double initialPSFWidth_rhs, double sigma_rhs) : verifier(FitPositions_EllipsoidalGaussian(initialPSFWidth_rhs, sigma_rhs)) {}
	~ParticleVerifier_EllipsoidalGaussian() {;}
	
	
	void VerifyParticles(boost::shared_ptr<Eigen::MatrixXd> image, boost::shared_ptr<std::list<position> > positions) {
		verifier.fit_positions(image, positions);
	}
protected:
	FitPositions_EllipsoidalGaussian verifier;
};

#endif

