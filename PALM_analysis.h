#ifndef PALM_ANALYSIS_H
#define PALM_ANALYSIS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "boost/thread.hpp"
#include "boost/bind.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include "PALM_analysis_errors.h"
#include "PALM_analysis_FileIO.h"
#include "PALM_analysis_defines.h"
#include "PALM_analysis_segmentation.h"
#include "PALM_analysis_ParticleFinding.h"
#include "PALM_analysis_Localization.h"
#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/filesystem.hpp"
#include "boost/filesystem/fstream.hpp"

#ifdef WITH_IGOR
#include "XOPStandardHeaders.h"
#endif

#define GSL_RANGE_CHECK_OFF	// this is not required since ublas::matrix<double> does range checks

namespace ublas = boost::numeric::ublas;

class ImageLoader;
class ThresholdImage;
class ThresholdImage_Preprocessor;
class ThresholdImage_Postprocessor;
class ParticleFinder;
class FitPositions;
class PALMAnalysisProgressReporter;


boost::shared_ptr<ublas::matrix <unsigned char> > do_processing_and_thresholding(boost::shared_ptr<ublas::matrix<double> > image, boost::shared_ptr<ThresholdImage_Preprocessor>preprocessor, 
																		 boost::shared_ptr<ThresholdImage> thresholder, boost::shared_ptr<ThresholdImage_Postprocessor> postprocessor);

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
	virtual double getIntegral(size_t index) const = 0;
	virtual double getXWidth(size_t index) const = 0;
	virtual double getYWidth(size_t index) const = 0;
	virtual double getCorrelation(size_t index) const = 0;
	virtual double getXPosition(size_t index) const = 0;
	virtual double getYPosition(size_t index) const = 0;
	virtual double getZPosition(size_t index) const = 0;
	virtual double getBackground(size_t index) const = 0;
	
	
	virtual double getIntegralDeviation(size_t index) const = 0;
	virtual double getXWidthDeviation(size_t index) const = 0;
	virtual double getYWidthDeviation(size_t index) const = 0;
	virtual double getCorrelationDeviation(size_t index) const = 0;
	virtual double getXPositionDeviation(size_t index) const = 0;
	virtual double getYPositionDeviation(size_t index) const = 0;
	virtual double getZPositionDeviation(size_t index) const = 0;
	virtual double getBackgroundDeviation(size_t index) const = 0;
	
	// add positions
	virtual void addPosition(boost::shared_ptr<LocalizedPosition> newPosition) = 0;
	virtual void addPositions(boost::shared_ptr<LocalizedPositionsContainer> newPositionsContainer) = 0;
	
	// set the frame numbers for all positions
	virtual void setFrameNumbers(size_t frameNumber) = 0;
	
	// save the positions localized this far under different formats
#ifdef WITH_IGOR
	virtual waveHndl writePositionsToWave(DataFolderAndName outputWaveParams, std::string waveNote) const = 0;
#endif
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
	double getCorrelation(size_t index) const {return 0;}
	double getXPosition(size_t index) const {return positionsVector.at(index).xPosition;}
	double getYPosition(size_t index) const {return positionsVector.at(index).yPosition;}
	double getZPosition(size_t index) const {return 0;}
	double getBackground(size_t index) const {return positionsVector.at(index).background;}
	
	double getIntegralDeviation(size_t index) const {return positionsVector.at(index).integralDeviation;}
	double getXWidthDeviation(size_t index) const {return positionsVector.at(index).widthDeviation;}
	double getYWidthDeviation(size_t index) const {return positionsVector.at(index).widthDeviation;}
	double getCorrelationDeviation(size_t index) const {return 0;}
	double getXPositionDeviation(size_t index) const {return positionsVector.at(index).xPositionDeviation;}
	double getYPositionDeviation(size_t index) const {return positionsVector.at(index).yPositionDeviation;}
	double getZPositionDeviation(size_t index) const {return 0;}
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
	
#ifdef WITH_IGOR
	waveHndl writePositionsToWave(DataFolderAndName outputWaveParams, std::string waveNote) const;
#endif
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
	double getXWidth(size_t index) const {return 0;}
	double getYWidth(size_t index) const {return 0;}
	double getCorrelation(size_t index) const {return 0;}
	double getXPosition(size_t index) const {return positionsVector.at(index).xPosition;}
	double getYPosition(size_t index) const {return positionsVector.at(index).yPosition;}
	double getZPosition(size_t index) const {return 0;}
	double getBackground(size_t index) const {return positionsVector.at(index).background;}
	
	double getIntegralDeviation(size_t index) const {return positionsVector.at(index).integralDeviation;}
	double getXWidthDeviation(size_t index) const {return 0;}
	double getYWidthDeviation(size_t index) const {return 0;}
	double getCorrelationDeviation(size_t index) const {return 0;}
	double getXPositionDeviation(size_t index) const {return positionsVector.at(index).xPositionDeviation;}
	double getYPositionDeviation(size_t index) const {return positionsVector.at(index).yPositionDeviation;}
	double getZPositionDeviation(size_t index) const {return 0;}
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
	
#ifdef WITH_IGOR
	waveHndl writePositionsToWave(DataFolderAndName outputWaveParams, std::string waveNote) const;
#endif
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
	double getZPosition(size_t index) const {return 0;}
	double getBackground(size_t index) const {return positionsVector.at(index).background;}
	
	double getIntegralDeviation(size_t index) const {return positionsVector.at(index).integralDeviation;}
	double getXWidthDeviation(size_t index) const {return positionsVector.at(index).xWidthDeviation;}
	double getYWidthDeviation(size_t index) const {return positionsVector.at(index).yWidthDeviation;}
	double getCorrelationDeviation(size_t index) const {return positionsVector.at(index).correlationDeviation;}
	double getXPositionDeviation(size_t index) const {return positionsVector.at(index).xPositionDeviation;}
	double getYPositionDeviation(size_t index) const {return positionsVector.at(index).yPositionDeviation;}
	double getZPositionDeviation(size_t index) const {return 0;}
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
	
#ifdef WITH_IGOR
	waveHndl writePositionsToWave(DataFolderAndName outputWaveParams, std::string waveNote) const;
#endif
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
	double getIntegral(size_t index) const {return 0;}
	double getXWidth(size_t index) const {return 0;}
	double getYWidth(size_t index) const {return 0;}
	double getCorrelation(size_t index) const {return 0;}
	double getXPosition(size_t index) const {return positionsVector.at(index).xPosition;}
	double getYPosition(size_t index) const {return positionsVector.at(index).yPosition;}
	double getZPosition(size_t index) const {return 0;}
	double getBackground(size_t index) const {return 0;}
	
	double getIntegralDeviation(size_t index) const {return 0;}
	double getXWidthDeviation(size_t index) const {return 0;}
	double getYWidthDeviation(size_t index) const {return 0;}
	double getCorrelationDeviation(size_t index) const {return 0;}
	double getXPositionDeviation(size_t index) const {return 0;}
	double getYPositionDeviation(size_t index) const {return 0;}
	double getZPositionDeviation(size_t index) const {return 0;}
	double getBackgroundDeviation(size_t index) const {return 0;}
	
	// adding new positions
	void addPosition(boost::shared_ptr<LocalizedPosition> newPosition);
	void addPositions(boost::shared_ptr<LocalizedPositionsContainer> newPositionsContainer);
	
	// set the frame numbers for all positions
	void setFrameNumbers(size_t frameNumber) {
		for (std::vector<LocalizedPosition_Centroid>::iterator it = this->positionsVector.begin(); it != this->positionsVector.end(); ++it) {
			(*it).frameNumber = frameNumber;
		}
	}
	
#ifdef WITH_IGOR
	waveHndl writePositionsToWave(DataFolderAndName outputWaveParams, std::string waveNote) const;
#endif
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
	double getIntegral(size_t index) const {return 0;}
	double getXWidth(size_t index) const {return positionsVector.at(index).width;}
	double getYWidth(size_t index) const {return positionsVector.at(index).width;}
	double getCorrelation(size_t index) const {return 0;}
	double getXPosition(size_t index) const {return positionsVector.at(index).xPosition;}
	double getYPosition(size_t index) const {return positionsVector.at(index).yPosition;}
	double getZPosition(size_t index) const {return 0;}
	double getBackground(size_t index) const {return 0;}
	
	double getIntegralDeviation(size_t index) const {return 0;}
	double getXWidthDeviation(size_t index) const {return 0;}
	double getYWidthDeviation(size_t index) const {return 0;}
	double getCorrelationDeviation(size_t index) const {return 0;}
	double getXPositionDeviation(size_t index) const {return 0;}
	double getYPositionDeviation(size_t index) const {return 0;}
	double getZPositionDeviation(size_t index) const {return 0;}
	double getBackgroundDeviation(size_t index) const {return 0;}
	
	// adding new positions
	void addPosition(boost::shared_ptr<LocalizedPosition> newPosition);
	void addPositions(boost::shared_ptr<LocalizedPositionsContainer> newPositionsContainer);
	
	// set the frame numbers for all positions
	void setFrameNumbers(size_t frameNumber) {
		for (std::vector<LocalizedPosition_Multiplication>::iterator it = this->positionsVector.begin(); it != this->positionsVector.end(); ++it) {
			(*it).frameNumber = frameNumber;
		}
	}
	
#ifdef WITH_IGOR
	waveHndl writePositionsToWave(DataFolderAndName outputWaveParams, std::string waveNote) const;
#endif
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
	double getXWidth(size_t index) const {return 0;}
	double getYWidth(size_t index) const {return 0;}
	double getCorrelation(size_t index) const {return 0;}
	double getXPosition(size_t index) const {return positionsVector.at(index).xPosition;}
	double getYPosition(size_t index) const {return positionsVector.at(index).yPosition;}
	double getZPosition(size_t index) const {return 0;}
	double getBackground(size_t index) const {return 0;}
	
	double getIntegralDeviation(size_t index) const {return 0;}
	double getXWidthDeviation(size_t index) const {return 0;}
	double getYWidthDeviation(size_t index) const {return 0;}
	double getCorrelationDeviation(size_t index) const {return 0;}
	double getXPositionDeviation(size_t index) const {return positionsVector.at(index).positionDeviation;}
	double getYPositionDeviation(size_t index) const {return positionsVector.at(index).positionDeviation;}
	double getZPositionDeviation(size_t index) const {return 0;}
	double getBackgroundDeviation(size_t index) const {return 0;}
	
	// adding new positions
	void addPosition(boost::shared_ptr<LocalizedPosition> newPosition);
	void addPositions(boost::shared_ptr<LocalizedPositionsContainer> newPositionsContainer);
	
	// set the frame numbers for all positions
	void setFrameNumbers(size_t frameNumber) {
		for (std::vector<LocalizedPosition_ZeissPALM>::iterator it = this->positionsVector.begin(); it != this->positionsVector.end(); ++it) {
			(*it).frameNumber = frameNumber;
		}
	}
	
#ifdef WITH_IGOR
	waveHndl writePositionsToWave(DataFolderAndName outputWaveParams, std::string waveNote) const;
#endif
	void writePositionsToFile(std::string filePath, std::string header) const {throw std::runtime_error("There is no meaning to write Zeiss positions to a file");}
	
	void sortPositionsByFrameNumber() {std::sort(positionsVector.begin(), positionsVector.end(), sortCompareFrameNumber);}
	static int sortCompareFrameNumber(LocalizedPosition_ZeissPALM left, LocalizedPosition_ZeissPALM right) {
		return ((left.frameNumber < right.frameNumber) ? 1 : 0);}
	
	size_t getPositionsType() const {return LOCALIZED_POSITIONS_TYPE_ZEISSPALM;}
	
protected:
	std::vector<LocalizedPosition_ZeissPALM> positionsVector;
};

/**
 *@brief A class that handles the actual PALM analysis
 *
 *
 * The 'PALMAnalysisController' class must be initialized with pointers to the thresholder, particle finder, fitter, and output writer objects.
 * It will coordinate the PALM analysis between those objects and take care of all the details
 */
class PALMAnalysisController {
public:
	PALMAnalysisController(boost::shared_ptr<ThresholdImage> thresholder_rhs,
						   boost::shared_ptr<ThresholdImage_Preprocessor> thresholdImagePreprocessor_rhs,
						   boost::shared_ptr<ThresholdImage_Postprocessor> thresholdImagePostprocessor_rhs,
						   boost::shared_ptr<ParticleFinder> particleFinder_rhs, boost::shared_ptr<FitPositions> fitPositions_rhs,
						   boost::shared_ptr<PALMAnalysisProgressReporter> progressReporter_rhs);
	~PALMAnalysisController() {;}
	
	boost::shared_ptr<LocalizedPositionsContainer> DoPALMAnalysis(boost::shared_ptr<ImageLoader> imageLoader_rhs);
	// runs a PALM analysis according to the parameters passed in as objects in the constructor
	// on the data file represented by the imageloader
	// returns a LocalizedPositionsContainer containing fitted PALM positions
	
protected:
	size_t nImages;
	
	std::queue<size_t> framesToBeProcessed;
	boost::shared_ptr<LocalizedPositionsContainer> localizedPositions;
	
	friend void ThreadPoolWorker(PALMAnalysisController* controller);
	
	boost::shared_ptr<ImageLoader> imageLoader;
	boost::shared_ptr<ThresholdImage> thresholder;
	boost::shared_ptr<ThresholdImage_Preprocessor> thresholdImagePreprocessor;
	boost::shared_ptr<ThresholdImage_Postprocessor> thresholdImagePostprocessor;
	boost::shared_ptr<ParticleFinder> particleFinder;
	boost::shared_ptr<FitPositions> fitPositions;
	boost::shared_ptr<PALMAnalysisProgressReporter> progressReporter;
	
	boost::mutex acquireFrameForProcessingMutex;
	boost::mutex addLocalizedPositionsMutex;
	boost::mutex errorReportingMutex;
	
	std::string errorMessage;
		// this string is an ugly hack to ensure that we can communicate errors encountered during the fitting back to the main thread
		// if one of the threads encounters an exception then it will set this message to some not-nil string
		// that is the sign for the main thread to kill the processing threads and throw an exception in the main thread
};

/**
 * @brief A worker function used by PALMAnalysisController and run in a separate thread. The function active requests images to process from PALMAnalysisController
 * and returns the results to the object
 */
void ThreadPoolWorker(PALMAnalysisController* controller);

/**
 * @briefWhen passed to PALMAnalysisController the methods of this object will be called to provide usage updates to the user
 */
class PALMAnalysisProgressReporter {
public:
	PALMAnalysisProgressReporter() {;}
	virtual ~PALMAnalysisProgressReporter() {;}
	
	virtual void CalculationStarted() = 0;
	virtual void UpdateCalculationProgress(double percentDone) = 0;
	virtual void CalculationDone() = 0;
	virtual void CalculationAborted() = 0;
};

class PALMAnalysisProgressReporter_Silent : public PALMAnalysisProgressReporter {
public:
	PALMAnalysisProgressReporter_Silent() {;}
	~PALMAnalysisProgressReporter_Silent() {;}
	
	void CalculationStarted() {;}
	void UpdateCalculationProgress(double percentDone) {;}
	void CalculationDone() {;}
	void CalculationAborted() {;}
	
protected:
	double previousPercentage;
};

#ifdef WITH_IGOR
/**
 * @brief Print updates on the calculation progress to the Igor command line
 */
class PALMAnalysisProgressReporter_IgorCommandLine : public PALMAnalysisProgressReporter {
public:
	PALMAnalysisProgressReporter_IgorCommandLine() {;}
	~PALMAnalysisProgressReporter_IgorCommandLine() {;}
	
	void CalculationStarted() {previousPercentage = 0; XOPNotice("Running calculation... ");}
	void UpdateCalculationProgress(double percentDone);
	void CalculationDone() {XOPNotice("Calculation finished!\r");}
	void CalculationAborted() {XOPNotice("Abort requested by user\r");}
	
protected:
	double previousPercentage;
};
#endif // WITH_IGOR

class PALMAnalysisProgressReporter_stdout : public PALMAnalysisProgressReporter {
public:
	PALMAnalysisProgressReporter_stdout() {;}
	~PALMAnalysisProgressReporter_stdout() {;}
	
	void CalculationStarted() {previousPercentage = 0; std::cout << "Running calculation... "; std::cout.flush();}
	void UpdateCalculationProgress(double percentDone);
	void CalculationDone() {std::cout << "Calculation finished!\n"; std::cout.flush();}
	void CalculationAborted() {std::cout << "Abort requested by user\n"; std::cout.flush();}
	
protected:
	double previousPercentage;
};

#endif
