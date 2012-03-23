/*
 Copyright 2008-2011 Peter Dedecker.
 
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

#ifndef PALM_ANALYSIS_STORAGE_H
#define PALM_ANALYSIS_STORAGE_H

#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include "boost/smart_ptr.hpp"
#include "PALM_analysis_defines.h"
#include "PALM_analysis_errors.h"

#ifdef WITH_IGOR
	#include "XOPStandardHeaders.h"
#endif

/**
 * @brief A simple container holding the coordinates of a single point in 3D space
 */
class Point3D {
public:
	double xPosition;
	double yPosition;
	double zPosition;
};

/**
 * @brief A more complex container containing the coordinates of a single point in 2D space
 * with its estimated intensity and background
 */
class Particle {
public:
	Particle() {x = 0; y = 0; intensity = 0;}
	Particle(double xLoc, double yLoc) {x = xLoc; y = yLoc;}
	Particle(double xLoc, double yLoc, double intensity_rhs) {x = xLoc; y = yLoc; intensity = intensity_rhs;}
	~Particle() {;}
	
	double x;
	double y;
	double intensity;
	double background;
};

/**
 *@brief A class designed to contain a single PALM position and its localization parameters
 */
class LocalizedPosition {
public:
	LocalizedPosition() {;}
	virtual ~LocalizedPosition() {;}
	
	virtual size_t getPositionType() const = 0;
};

/**
 * Stores a single position localized using fitting of a circularly symmetric Gaussian
 */
class LocalizedPosition_2DGauss : public LocalizedPosition {
public:
	LocalizedPosition_2DGauss() : nFramesPresent(1) {}
	~LocalizedPosition_2DGauss() {;}
	
	size_t frameNumber;
	size_t nFramesPresent;	// the number of frames this position was localized in
	
	double integral;
	double width;
	double xPosition;
	double yPosition;
	double background;
	
	double integralDeviation;
	double widthDeviation;
	double xPositionDeviation;
	double yPositionDeviation;
	double backgroundDeviation;
	
	size_t getPositionType() const {return LOCALIZED_POSITIONS_TYPE_2DGAUSS;}
};

class LocalizedPosition_2DGaussFixedWidth : public LocalizedPosition {
public:
	LocalizedPosition_2DGaussFixedWidth() : nFramesPresent(1) {}
	~LocalizedPosition_2DGaussFixedWidth() {;}
	
	size_t frameNumber;
	size_t nFramesPresent;	// the number of frames this position was localized in
	
	double integral;
	double xPosition;
	double yPosition;
	double background;
	
	double integralDeviation;
	double xPositionDeviation;
	double yPositionDeviation;
	double backgroundDeviation;
	
	size_t getPositionType() const {return LOCALIZED_POSITIONS_TYPE_2DGAUSS_FIXED_WIDTH;}
};

class LocalizedPosition_Ellipsoidal2DGauss : public LocalizedPosition {
public:
	LocalizedPosition_Ellipsoidal2DGauss() : nFramesPresent(1) {}
	~LocalizedPosition_Ellipsoidal2DGauss() {;}
	
	size_t frameNumber;
	size_t nFramesPresent;	// the number of frames this position was localized in
	
	double integral;
	double xWidth;
	double yWidth;
	double xPosition;
	double yPosition;
	double correlation;
	double background;
	
	double integralDeviation;
	double xWidthDeviation;
	double yWidthDeviation;
	double xPositionDeviation;
	double yPositionDeviation;
	double correlationDeviation;
	double backgroundDeviation;
	
	size_t getPositionType() const {return LOCALIZED_POSITIONS_TYPE_ELLIPSOIDAL2DGAUSS;}
};

class LocalizedPosition_Centroid : public LocalizedPosition {
public:
	LocalizedPosition_Centroid() : nFramesPresent(1) {}
	~LocalizedPosition_Centroid() {;}
	
	size_t frameNumber;
	size_t nFramesPresent;	// the number of frames this position was localized in
	
	double xPosition;
	double yPosition;
	
	size_t getPositionType() const {return LOCALIZED_POSITIONS_TYPE_CENTROID;}
};

class LocalizedPosition_Multiplication : public LocalizedPosition {
public:
	LocalizedPosition_Multiplication() : nFramesPresent(1) {}
	~LocalizedPosition_Multiplication() {;}
	
	size_t frameNumber;
	size_t nFramesPresent;	// the number of frames this position was localized in
	
	double xPosition;
	double yPosition;
	double width;
	
	size_t getPositionType() const {return LOCALIZED_POSITIONS_TYPE_MULTIPLICATION;}
};

class LocalizedPosition_ZeissPALM : public LocalizedPosition {
public:
	LocalizedPosition_ZeissPALM() : nFramesPresent(1) {}
	~LocalizedPosition_ZeissPALM() {;}
	
	size_t frameNumber;
	size_t nFramesPresent;	// the number of frames this position was localized in
	
	double integral;
	double xPosition;
	double yPosition;
	double positionDeviation;
	
	size_t getPositionType() const {return LOCALIZED_POSITIONS_TYPE_ZEISSPALM;}
};

class LocalizedPosition_MLEwG : public LocalizedPosition {
public:
	LocalizedPosition_MLEwG() : nFramesPresent(1) {}
	~LocalizedPosition_MLEwG() {;}
	
	size_t frameNumber;
	size_t nFramesPresent;	// the number of frames this position was localized in
	
	double integral;
	double xPosition;
	double yPosition;
	double positionDeviation;
	double width;
	double background;
	
	size_t getPositionType() const {return LOCALIZED_POSITIONS_TYPE_MLEWG;}
};

/**
 * @brief An abstract base class that holds localized positions. Has derived classes that handle specific types of positions.
 * Can import positions from waves or files, and write to them. The base class contains accessor methods for every data possible,
 * derived classes should return meaningful results for all of these, depending on the type of localization used
 */
class LocalizedPositionsContainer {
public:
	// these functions handle creating a LocalizedPositionsContainer object from an igor wave
	// or from a file containing positions written to disk
	// the functions will discern the type of positions and return a LocalizedPositionsContainer of the correct type
#ifdef WITH_IGOR
	static boost::shared_ptr<LocalizedPositionsContainer> GetPositionsFromWave(waveHndl positionsWave);
#endif
	static boost::shared_ptr<LocalizedPositionsContainer> GetPositionsFromFile(std::string filePath);
	
	// constructor and destructor
	LocalizedPositionsContainer() {;}
	virtual ~LocalizedPositionsContainer() {;}
	
	// accessor methods
	virtual size_t getNPositions() const = 0;
	virtual size_t getFrameNumber(size_t index) const = 0;
	virtual double getIntegral(size_t index) const {return 0;}
	virtual double getXWidth(size_t index) const {return 0;}
	virtual double getYWidth(size_t index) const {return 0;}
	virtual double getCorrelation(size_t index) const {return 0;}
	virtual double getXPosition(size_t index) const {return 0;}
	virtual double getYPosition(size_t index) const {return 0;}
	virtual double getZPosition(size_t index) const {return 0;}
	virtual double getBackground(size_t index) const {return 0;}
	
	
	virtual double getIntegralDeviation(size_t index) const {return 0;}
	virtual double getXWidthDeviation(size_t index) const {return 0;}
	virtual double getYWidthDeviation(size_t index) const {return 0;}
	virtual double getCorrelationDeviation(size_t index) const {return 0;}
	virtual double getXPositionDeviation(size_t index) const {return 0;}
	virtual double getYPositionDeviation(size_t index) const {return 0;}
	virtual double getZPositionDeviation(size_t index) const {return 0;}
	virtual double getBackgroundDeviation(size_t index) const {return 0;}
	
	// add positions
	virtual void addPosition(boost::shared_ptr<LocalizedPosition> newPosition) = 0;
	virtual void addPositions(boost::shared_ptr<LocalizedPositionsContainer> newPositionsContainer) = 0;
	
	// set the frame numbers for all positions
	virtual void setFrameNumbers(size_t frameNumber) = 0;
	
	// save the positions localized this far under different formats
	virtual ImagePtr getLocalizedPositionsAsMatrix() const = 0;
	virtual void writePositionsToFile(std::string filePath, std::string waveNote) const = 0;
	
	// sort the positions according to ascending frame number, with no guarantees
	// for the sorting of positions within the same frame
	virtual void sortPositionsByFrameNumber() = 0;
	
	virtual size_t getPositionsType() const = 0;
};

/**
 * @brief Contain localized positions created by fitting a circularly symmetric 2D Gaussian
 */
class LocalizedPositionsContainer_2DGauss : public LocalizedPositionsContainer {
public:
	LocalizedPositionsContainer_2DGauss() {;}
#ifdef WITH_IGOR
	LocalizedPositionsContainer_2DGauss(waveHndl wave);
#endif
	LocalizedPositionsContainer_2DGauss(const std::string filePath) {throw std::runtime_error("Loading positions from files is not yet supported");}
	
	~LocalizedPositionsContainer_2DGauss() {;}
	
	// accessor methods
	size_t getNPositions() const {return positionsVector.size();}
	size_t getFrameNumber(size_t index) const {return positionsVector.at(index).frameNumber;}
	double getIntegral(size_t index) const {return positionsVector.at(index).integral;}
	double getXWidth(size_t index) const {return positionsVector.at(index).width;}
	double getYWidth(size_t index) const {return positionsVector.at(index).width;}
	double getXPosition(size_t index) const {return positionsVector.at(index).xPosition;}
	double getYPosition(size_t index) const {return positionsVector.at(index).yPosition;}
	double getZPosition(size_t index) const {return 0;}
	double getBackground(size_t index) const {return positionsVector.at(index).background;}
	
	double getIntegralDeviation(size_t index) const {return positionsVector.at(index).integralDeviation;}
	double getXWidthDeviation(size_t index) const {return positionsVector.at(index).widthDeviation;}
	double getYWidthDeviation(size_t index) const {return positionsVector.at(index).widthDeviation;}
	double getXPositionDeviation(size_t index) const {return positionsVector.at(index).xPositionDeviation;}
	double getYPositionDeviation(size_t index) const {return positionsVector.at(index).yPositionDeviation;}
	double getBackgroundDeviation(size_t index) const {return positionsVector.at(index).backgroundDeviation;}
	
	// adding new positions
	void addPosition(boost::shared_ptr<LocalizedPosition> newPosition);
	void addPositions(boost::shared_ptr<LocalizedPositionsContainer> newPositionsContainer);
	
	// set the frame numbers for all positions
	void setFrameNumbers(size_t frameNumber) {
		for (std::vector<LocalizedPosition_2DGauss>::iterator it = this->positionsVector.begin(); it != this->positionsVector.end(); ++it) {
			(*it).frameNumber = frameNumber;
		}
	}
	
	ImagePtr getLocalizedPositionsAsMatrix() const;
	void writePositionsToFile(std::string filePath, std::string header) const;
	
	void sortPositionsByFrameNumber() {std::sort(positionsVector.begin(), positionsVector.end(), sortCompareFrameNumber);}
	static int sortCompareFrameNumber(LocalizedPosition_2DGauss left, LocalizedPosition_2DGauss right) {
		return ((left.frameNumber < right.frameNumber) ? 1 : 0);}
	
	size_t getPositionsType() const {return LOCALIZED_POSITIONS_TYPE_2DGAUSS;}
	
protected:
	std::vector<LocalizedPosition_2DGauss> positionsVector;
};

class LocalizedPositionsContainer_2DGaussFixedWidth : public LocalizedPositionsContainer {
public:
	LocalizedPositionsContainer_2DGaussFixedWidth() {;}
#ifdef WITH_IGOR
	LocalizedPositionsContainer_2DGaussFixedWidth(waveHndl wave);
#endif
	LocalizedPositionsContainer_2DGaussFixedWidth(const std::string filePath) {throw std::runtime_error("Loading positions from files is not yet supported");}
	
	~LocalizedPositionsContainer_2DGaussFixedWidth() {;}
	
	// accessor methods
	size_t getNPositions() const {return positionsVector.size();}
	size_t getFrameNumber(size_t index) const {return positionsVector.at(index).frameNumber;}
	double getIntegral(size_t index) const {return positionsVector.at(index).integral;}
	double getXPosition(size_t index) const {return positionsVector.at(index).xPosition;}
	double getYPosition(size_t index) const {return positionsVector.at(index).yPosition;}
	double getBackground(size_t index) const {return positionsVector.at(index).background;}
	double getIntegralDeviation(size_t index) const {return positionsVector.at(index).integralDeviation;}
	double getXPositionDeviation(size_t index) const {return positionsVector.at(index).xPositionDeviation;}
	double getYPositionDeviation(size_t index) const {return positionsVector.at(index).yPositionDeviation;}
	double getBackgroundDeviation(size_t index) const {return positionsVector.at(index).backgroundDeviation;}
	
	// adding new positions
	void addPosition(boost::shared_ptr<LocalizedPosition> newPosition);
	void addPositions(boost::shared_ptr<LocalizedPositionsContainer> newPositionsContainer);
	
	// set the frame numbers for all positions
	void setFrameNumbers(size_t frameNumber) {
		for (std::vector<LocalizedPosition_2DGaussFixedWidth>::iterator it = this->positionsVector.begin(); it != this->positionsVector.end(); ++it) {
			(*it).frameNumber = frameNumber;
		}
	}
	
	ImagePtr getLocalizedPositionsAsMatrix() const;
	void writePositionsToFile(std::string filePath, std::string header) const;
	
	void sortPositionsByFrameNumber() {std::sort(positionsVector.begin(), positionsVector.end(), sortCompareFrameNumber);}
	static int sortCompareFrameNumber(LocalizedPosition_2DGaussFixedWidth left, LocalizedPosition_2DGaussFixedWidth right) {
		return ((left.frameNumber < right.frameNumber) ? 1 : 0);}
	
	size_t getPositionsType() const {return LOCALIZED_POSITIONS_TYPE_2DGAUSS_FIXED_WIDTH;}
	
protected:
	std::vector<LocalizedPosition_2DGaussFixedWidth> positionsVector;
};

/**
 * @brief Contain localized positions created by fitting an ellipsoidal 2D Gaussian
 */
class LocalizedPositionsContainer_Ellipsoidal2DGaussian : public LocalizedPositionsContainer {
public:
	LocalizedPositionsContainer_Ellipsoidal2DGaussian() {;}
#ifdef WITH_IGOR
	LocalizedPositionsContainer_Ellipsoidal2DGaussian(waveHndl wave);
#endif
	LocalizedPositionsContainer_Ellipsoidal2DGaussian(const std::string filePath) {throw std::runtime_error("Loading positions from files is not yet supported");}
	
	~LocalizedPositionsContainer_Ellipsoidal2DGaussian() {;}
	
	// accessor methods
	size_t getNPositions() const {return positionsVector.size();}
	size_t getFrameNumber(size_t index) const {return positionsVector.at(index).frameNumber;}
	double getIntegral(size_t index) const {return positionsVector.at(index).integral;}
	double getXWidth(size_t index) const {return positionsVector.at(index).xWidth;}
	double getYWidth(size_t index) const {return positionsVector.at(index).yWidth;}
	double getCorrelation(size_t index) const {return positionsVector.at(index).correlation;}
	double getXPosition(size_t index) const {return positionsVector.at(index).xPosition;}
	double getYPosition(size_t index) const {return positionsVector.at(index).yPosition;}
	double getBackground(size_t index) const {return positionsVector.at(index).background;}
	
	double getIntegralDeviation(size_t index) const {return positionsVector.at(index).integralDeviation;}
	double getXWidthDeviation(size_t index) const {return positionsVector.at(index).xWidthDeviation;}
	double getYWidthDeviation(size_t index) const {return positionsVector.at(index).yWidthDeviation;}
	double getCorrelationDeviation(size_t index) const {return positionsVector.at(index).correlationDeviation;}
	double getXPositionDeviation(size_t index) const {return positionsVector.at(index).xPositionDeviation;}
	double getYPositionDeviation(size_t index) const {return positionsVector.at(index).yPositionDeviation;}
	double getBackgroundDeviation(size_t index) const {return positionsVector.at(index).backgroundDeviation;}
	
	// adding new positions
	void addPosition(boost::shared_ptr<LocalizedPosition> newPosition);
	void addPositions(boost::shared_ptr<LocalizedPositionsContainer> newPositionsContainer);
	
	// set the frame numbers for all positions
	void setFrameNumbers(size_t frameNumber) {
		for (std::vector<LocalizedPosition_Ellipsoidal2DGauss>::iterator it = this->positionsVector.begin(); it != this->positionsVector.end(); ++it) {
			(*it).frameNumber = frameNumber;
		}
	}
	
	ImagePtr getLocalizedPositionsAsMatrix() const;
	void writePositionsToFile(std::string filePath, std::string header) const;
	
	void sortPositionsByFrameNumber() {std::sort(positionsVector.begin(), positionsVector.end(), sortCompareFrameNumber);}
	static int sortCompareFrameNumber(LocalizedPosition_Ellipsoidal2DGauss left, LocalizedPosition_Ellipsoidal2DGauss right) {
		return ((left.frameNumber < right.frameNumber) ? 1 : 0);}
	
	size_t getPositionsType() const {return LOCALIZED_POSITIONS_TYPE_ELLIPSOIDAL2DGAUSS;}
	
protected:
	std::vector<LocalizedPosition_Ellipsoidal2DGauss> positionsVector;
};

class LocalizedPositionsContainer_Centroid : public LocalizedPositionsContainer {
public:
	LocalizedPositionsContainer_Centroid() {;}
#ifdef WITH_IGOR
	LocalizedPositionsContainer_Centroid(waveHndl wave);
#endif
	LocalizedPositionsContainer_Centroid(const std::string filePath) {throw std::runtime_error("Loading positions from files is not yet supported");}
	
	~LocalizedPositionsContainer_Centroid() {;}
	
	// accessor methods
	size_t getNPositions() const {return positionsVector.size();}
	size_t getFrameNumber(size_t index) const {return positionsVector.at(index).frameNumber;}
	double getXPosition(size_t index) const {return positionsVector.at(index).xPosition;}
	double getYPosition(size_t index) const {return positionsVector.at(index).yPosition;}
	
	// adding new positions
	void addPosition(boost::shared_ptr<LocalizedPosition> newPosition);
	void addPositions(boost::shared_ptr<LocalizedPositionsContainer> newPositionsContainer);
	
	// set the frame numbers for all positions
	void setFrameNumbers(size_t frameNumber) {
		for (std::vector<LocalizedPosition_Centroid>::iterator it = this->positionsVector.begin(); it != this->positionsVector.end(); ++it) {
			(*it).frameNumber = frameNumber;
		}
	}
	
	ImagePtr getLocalizedPositionsAsMatrix() const;
	void writePositionsToFile(std::string filePath, std::string header) const;
	
	void sortPositionsByFrameNumber() {std::sort(positionsVector.begin(), positionsVector.end(), sortCompareFrameNumber);}
	static int sortCompareFrameNumber(LocalizedPosition_Centroid left, LocalizedPosition_Centroid right) {
		return ((left.frameNumber < right.frameNumber) ? 1 : 0);}
	
	size_t getPositionsType() const {return LOCALIZED_POSITIONS_TYPE_CENTROID;}
	
protected:
	std::vector<LocalizedPosition_Centroid> positionsVector;
};

class LocalizedPositionsContainer_Multiplication : public LocalizedPositionsContainer {
public:
	LocalizedPositionsContainer_Multiplication() {;}
#ifdef WITH_IGOR
	LocalizedPositionsContainer_Multiplication(waveHndl wave);
#endif
	LocalizedPositionsContainer_Multiplication(const std::string filePath) {throw std::runtime_error("Loading positions from files is not yet supported");}
	
	~LocalizedPositionsContainer_Multiplication() {;}
	
	// accessor methods
	size_t getNPositions() const {return positionsVector.size();}
	size_t getFrameNumber(size_t index) const {return positionsVector.at(index).frameNumber;}
	double getXWidth(size_t index) const {return positionsVector.at(index).width;}
	double getYWidth(size_t index) const {return positionsVector.at(index).width;}
	double getXPosition(size_t index) const {return positionsVector.at(index).xPosition;}
	double getYPosition(size_t index) const {return positionsVector.at(index).yPosition;}
	
	// adding new positions
	void addPosition(boost::shared_ptr<LocalizedPosition> newPosition);
	void addPositions(boost::shared_ptr<LocalizedPositionsContainer> newPositionsContainer);
	
	// set the frame numbers for all positions
	void setFrameNumbers(size_t frameNumber) {
		for (std::vector<LocalizedPosition_Multiplication>::iterator it = this->positionsVector.begin(); it != this->positionsVector.end(); ++it) {
			(*it).frameNumber = frameNumber;
		}
	}
	
	ImagePtr getLocalizedPositionsAsMatrix() const;
	void writePositionsToFile(std::string filePath, std::string header) const;
	
	void sortPositionsByFrameNumber() {std::sort(positionsVector.begin(), positionsVector.end(), sortCompareFrameNumber);}
	static int sortCompareFrameNumber(LocalizedPosition_Multiplication left, LocalizedPosition_Multiplication right) {
		return ((left.frameNumber < right.frameNumber) ? 1 : 0);}
	
	size_t getPositionsType() const {return LOCALIZED_POSITIONS_TYPE_MULTIPLICATION;}
	
protected:
	std::vector<LocalizedPosition_Multiplication> positionsVector;
};

class LocalizedPositionsContainer_ZeissPALM : public LocalizedPositionsContainer {
public:
	LocalizedPositionsContainer_ZeissPALM() {;}
#ifdef WITH_IGOR
	LocalizedPositionsContainer_ZeissPALM(waveHndl wave);
#endif
	LocalizedPositionsContainer_ZeissPALM(const std::string filePath) {throw std::runtime_error("Loading positions from files is not yet supported");}
	
	~LocalizedPositionsContainer_ZeissPALM() {;}
	
	// accessor methods
	size_t getNPositions() const {return positionsVector.size();}
	size_t getFrameNumber(size_t index) const {return positionsVector.at(index).frameNumber;}
	double getIntegral(size_t index) const {return positionsVector.at(index).integral;}
	double getXPosition(size_t index) const {return positionsVector.at(index).xPosition;}
	double getYPosition(size_t index) const {return positionsVector.at(index).yPosition;}
	double getXPositionDeviation(size_t index) const {return positionsVector.at(index).positionDeviation;}
	double getYPositionDeviation(size_t index) const {return positionsVector.at(index).positionDeviation;}
	
	// adding new positions
	void addPosition(boost::shared_ptr<LocalizedPosition> newPosition);
	void addPositions(boost::shared_ptr<LocalizedPositionsContainer> newPositionsContainer);
	
	// set the frame numbers for all positions
	void setFrameNumbers(size_t frameNumber) {
		for (std::vector<LocalizedPosition_ZeissPALM>::iterator it = this->positionsVector.begin(); it != this->positionsVector.end(); ++it) {
			(*it).frameNumber = frameNumber;
		}
	}
	
	ImagePtr getLocalizedPositionsAsMatrix() const;
	void writePositionsToFile(std::string filePath, std::string header) const {throw std::runtime_error("There is no meaning to write Zeiss positions to a file");}
	
	void sortPositionsByFrameNumber() {std::sort(positionsVector.begin(), positionsVector.end(), sortCompareFrameNumber);}
	static int sortCompareFrameNumber(LocalizedPosition_ZeissPALM left, LocalizedPosition_ZeissPALM right) {
		return ((left.frameNumber < right.frameNumber) ? 1 : 0);}
	
	size_t getPositionsType() const {return LOCALIZED_POSITIONS_TYPE_ZEISSPALM;}
	
protected:
	std::vector<LocalizedPosition_ZeissPALM> positionsVector;
};

class LocalizedPositionsContainer_MLEwG : public LocalizedPositionsContainer {
public:
	LocalizedPositionsContainer_MLEwG() {;}
#ifdef WITH_IGOR
	LocalizedPositionsContainer_MLEwG(waveHndl wave);
#endif
	LocalizedPositionsContainer_MLEwG(const std::string filePath) {throw std::runtime_error("Loading positions from files is not yet supported");}
	
	~LocalizedPositionsContainer_MLEwG() {;}
	
	// accessor methods
	size_t getNPositions() const {return positionsVector.size();}
	size_t getFrameNumber(size_t index) const {return positionsVector.at(index).frameNumber;}
	double getIntegral(size_t index) const {return positionsVector.at(index).integral;}
	double getXWidth(size_t index) const {return positionsVector.at(index).width;}
	double getYWidth(size_t index) const {return positionsVector.at(index).width;}
	double getXPosition(size_t index) const {return positionsVector.at(index).xPosition;}
	double getYPosition(size_t index) const {return positionsVector.at(index).yPosition;}
	double getBackground(size_t index) const {return positionsVector.at(index).background;}
	
	double getXPositionDeviation(size_t index) const {return positionsVector.at(index).positionDeviation;}
	double getYPositionDeviation(size_t index) const {return positionsVector.at(index).positionDeviation;}
	
	// adding new positions
	void addPosition(boost::shared_ptr<LocalizedPosition> newPosition);
	void addPositions(boost::shared_ptr<LocalizedPositionsContainer> newPositionsContainer);
	
	// set the frame numbers for all positions
	void setFrameNumbers(size_t frameNumber) {
		for (std::vector<LocalizedPosition_MLEwG>::iterator it = this->positionsVector.begin(); it != this->positionsVector.end(); ++it) {
			(*it).frameNumber = frameNumber;
		}
	}
	
	ImagePtr getLocalizedPositionsAsMatrix() const;
	void writePositionsToFile(std::string filePath, std::string header) const;
	
	void sortPositionsByFrameNumber() {std::sort(positionsVector.begin(), positionsVector.end(), sortCompareFrameNumber);}
	static int sortCompareFrameNumber(LocalizedPosition_MLEwG left, LocalizedPosition_MLEwG right) {
		return ((left.frameNumber < right.frameNumber) ? 1 : 0);}
	
	size_t getPositionsType() const {return LOCALIZED_POSITIONS_TYPE_MLEWG;}
	
protected:
	std::vector<LocalizedPosition_MLEwG> positionsVector;
};

#endif
