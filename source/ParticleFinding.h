/*
 Copyright 2008-2014 Peter Dedecker.
 
 This file is part of Localizer.
 
 Localizer is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 Localizer is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with Localizer.  If not, see <http://www.gnu.org/licenses/>.
 
 
 Additional permission under GNU GPL version 3 section 7
 
 If you modify this Program, or any covered work, by
 linking or combining it with libraries required for interaction
 with analysis programs such as Igor Pro or Matlab, or to acquire
 data from or control hardware related to an experimental measurement,
 the licensors of this Program grant you additional permission
 to convey the resulting work.
 */

#ifndef PALM_ANALYSIS_PARTICLEFINDING_H
#define PALM_ANALYSIS_PARTICLEFINDING_H

#include <algorithm>
#include "PALM_analysis_storage.h"
#include "PALM_analysis_Localization.h"
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
	
	virtual std::shared_ptr<std::list<Particle> > findPositions(ImagePtr image, std::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > threshold_image) = 0;
protected:
};

/**
 * @brief Assumes that all active pixels within a given radius belong to the same particle
 */
class ParticleFinder_radius : public ParticleFinder {
public:
	ParticleFinder_radius(double radius_rhs) {radius = radius_rhs;}
	~ParticleFinder_radius() {;}
	
	std::shared_ptr<std::list<Particle> > findPositions(ImagePtr image, std::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > threshold_image);
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
	
	std::shared_ptr<std::list<Particle> > findPositions(ImagePtr image, std::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > threshold_image);
	
protected:
	void growParticle(Particle centerPosition, std::list<Particle> &positionsInCurrentParticle, ImagePtr image, std::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > threshold_image, std::shared_ptr<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > mapped_image);
	
};

/**
 * @brief Assumes that all neighbouring active pixels (eight-way adjacency) belong to the same particle
 */
class ParticleFinder_adjacent8 : public ParticleFinder {
public:
	ParticleFinder_adjacent8() {;}
	~ParticleFinder_adjacent8() {;}
	
	std::shared_ptr<std::list<Particle> > findPositions(ImagePtr image, std::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > threshold_image);
	
protected:
	void growParticle(Particle centerPosition, std::list<Particle> &positionsInCurrentParticle, ImagePtr image, std::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > threshold_image, std::shared_ptr<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > mapped_image);
	
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
	virtual void VerifyParticles(ImagePtr image, std::shared_ptr<std::list<Particle> > positions) = 0;
protected:
};

/**
 * Verify the particles by eliminating any two that are within 4 times the PSF deviation
 */
class ParticleVerifier_RemoveOverlappingParticles : public ParticleVerifier {
public:
	ParticleVerifier_RemoveOverlappingParticles(double psfWidth_rhs) : psfWidth(psfWidth_rhs) {}
	~ParticleVerifier_RemoveOverlappingParticles() {;}
	
	
	void VerifyParticles(ImagePtr image, std::shared_ptr<std::list<Particle> > positions);
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
	
	
	void VerifyParticles(ImagePtr image, std::shared_ptr<std::list<Particle> > positions) {
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
	
	
	void VerifyParticles(ImagePtr image, std::shared_ptr<std::list<Particle> > positions) {
		verifier.fit_positions(image, positions);
	}
protected:
	FitPositions_EllipsoidalGaussian verifier;
};

/**
 * Verify the particles by seeing if they pass the fitting of a FitPositions_EllipsoidalGaussian_SymmetricPSF class
 */
class ParticleVerifier_EllipsoidalGaussian_SymmetricPSF : public ParticleVerifier {
public:
	ParticleVerifier_EllipsoidalGaussian_SymmetricPSF(double initialPSFWidth_rhs, double sigma_rhs) : verifier(FitPositions_EllipsoidalGaussian_SymmetricPSF(initialPSFWidth_rhs, sigma_rhs)) {}
	~ParticleVerifier_EllipsoidalGaussian_SymmetricPSF() {;}
	
	
	void VerifyParticles(ImagePtr image, std::shared_ptr<std::list<Particle> > positions) {
		verifier.fit_positions(image, positions);
	}
protected:
	FitPositions_EllipsoidalGaussian_SymmetricPSF verifier;
};

#endif

