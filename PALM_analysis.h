#ifndef PALM_ANALYSIS_H
#define PALM_ANALYSIS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "boost/thread.hpp"
#include "boost/bind.hpp"
#include "XOPStandardHeaders.h"
#include "PALM_analysis_errors.h"
#include "PALM_analysis_FileIO.h"
#include "PALM_analysis_defines.h"
#include "PALM_analysis_segmentation.h"
#include "PALM_analysis_ParticleFinding.h"
#include "PALM_analysis_Localization.h"
#include "boost/date_time/posix_time/posix_time.hpp"

#define GSL_RANGE_CHECK_OFF	// this is not required since PALMMatrix<double> does range checks

using namespace std;
class ImageLoader;
class ThresholdImage;
class ThresholdImage_Preprocessor;
class ThresholdImage_Postprocessor;
class ParticleFinder;
class PALMResultsWriter;
class FitPositions;
class PALMAnalysisProgressReporter;

int load_partial_ccd_image(ImageLoader *image_loader, size_t n_start, size_t n_end);

int parse_ccd_headers(ImageLoader *image_loader);


boost::shared_ptr<PALMMatrix <unsigned char> > do_processing_and_thresholding(boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<ThresholdImage_Preprocessor>preprocessor, 
																		 boost::shared_ptr<ThresholdImage> thresholder, boost::shared_ptr<ThresholdImage_Postprocessor> postprocessor);

/**
 *@brief A class designed to contain a single PALM position and its localization parameters
 */
class LocalizedPosition {
public:
	LocalizedPosition() {amplitude = 0; width = 0; xPos = 0; yPos = 0; offset = 0; amplitudeError = 0, widthError = 0,
		xPosError = 0; yPosError = 0; offsetError = 0; nIterations = 0;}
	~LocalizedPosition() {;}
	
	double amplitude;
	double width;
	double xPos;
	double yPos;
	double offset;
	
	double amplitudeError;
	double widthError;
	double xPosError;
	double yPosError;
	double offsetError;
	
	size_t nIterations;
};


/**
 *@brief Stores the result from fitting a particular frame, remembering both the frame index as well as the fitted positions
 */
class PALMResults {
public:
	PALMResults(size_t frameIndex_rhs, boost::shared_ptr<std::vector<LocalizedPosition> > fittedPositions_rhs) {
		frameIndex = frameIndex_rhs;
		fittedPositions = fittedPositions_rhs;
	}
	~PALMResults() {;}
	
	size_t getFrameIndex() {return frameIndex;}
	size_t getNumberOfFittedPositions() {return fittedPositions->size();}
	boost::shared_ptr<std::vector<LocalizedPosition> > getFittedPositions() {return fittedPositions;}
	
protected:
	size_t frameIndex;
	boost::shared_ptr<std::vector<LocalizedPosition> > fittedPositions;
};

/**
 * Stores a single position localized using fitting of a circularly symmetric Gaussian
 */
class LocalizedPosition_2DGauss {
public:
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
};

class LocalizedPosition_2DGaussFixedWidth {
public:
	size_t frameNumber;
	size_t nFramesPresent;	// the number of frames this position was localized in
	
	double integral;
	double width;
	double xPosition;
	double yPosition;
	double background;
	
	double integralDeviation;
	double xPositionDeviation;
	double yPositionDeviation;
	double backgroundDeviation;
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
	static boost::shared_ptr<LocalizedPositionsContainer> GetPositionsFromWave(waveHndl positionsWave);
	static boost::shared_ptr<LocalizedPositionsContainer> GetPositionsFromFile(std::string& filePath);
	
	// constructor and destructor
	LocalizedPositionsContainer() {;}
	virtual ~LocalizedPositionsContainer() {;}
	
	// accessor methods
	virtual size_t getNPositions() const = 0;
	virtual size_t getFrameNumber(size_t index) const = 0;
	virtual double getIntegral(size_t index) const = 0;
	virtual double getXWidth(size_t index) const = 0;
	virtual double getYWidth(size_t index) const = 0;
	virtual double getRotationAngle(size_t index) const = 0;
	virtual double getXPosition(size_t index) const = 0;
	virtual double getYPosition(size_t index) const = 0;
	virtual double getZPosition(size_t index) const = 0;
	virtual double getBackground(size_t index) const = 0;
	
	
	virtual double getIntegralDeviation(size_t index) const = 0;
	virtual double getXWidthDeviation(size_t index) const = 0;
	virtual double getYWidthDeviation(size_t index) const = 0;
	virtual double getRotationAngleDeviation(size_t index) const = 0;
	virtual double getXPositionDeviation(size_t index) const = 0;
	virtual double getYPositionDeviation(size_t index) const = 0;
	virtual double getZPositionDeviation(size_t index) const = 0;
	virtual double getBackgroundDeviation(size_t index) const = 0;
	
	// this abstract base class defines no methods to add positions
	// that is left up to derived class for type safety
	// this means that it is impossible to simultaneously read and write from the same positions
	
	// save the positions localized this far under different formats
	virtual waveHndl writePositionsToWave(std::string waveName) const = 0;
	virtual void writePositionsToFile(std::string filePath) const = 0;
	
	// sort the positions according to ascending frame number, with no guarantees
	// for the sorting of positions within the same frame
	virtual void sortPositionsByFrameNumber() = 0;
};

/**
 * @brief Contain localized positions created by fitting a circularly symmetric 2D Gaussian
 */
class LocalizedPositionsContainer_2DGauss : public LocalizedPositionsContainer {
public:
	LocalizedPositionsContainer_2DGauss() {;}
	LocalizedPositionsContainer_2DGauss(waveHndl wave);
	LocalizedPositionsContainer_2DGauss(const std::string& filePath) {throw std::runtime_error("Loading positions from files is not yet supported");}
	
	~LocalizedPositionsContainer_2DGauss() {;}
	
	// accessor methods
	size_t getNPositions() const {return positionsVector.size();}
	size_t getFrameNumber(size_t index) const {return positionsVector.at(index).frameNumber;}
	double getIntegral(size_t index) const {return positionsVector.at(index).integral;}
	double getXWidth(size_t index) const {return positionsVector.at(index).width;}
	double getYWidth(size_t index) const {return positionsVector.at(index).width;}
	double getRotationAngle(size_t index) const {return 0;}
	double getXPosition(size_t index) const {return positionsVector.at(index).xPosition;}
	double getYPosition(size_t index) const {return positionsVector.at(index).yPosition;}
	double getZPosition(size_t index) const {return 0;}
	double getBackground(size_t index) const {return positionsVector.at(index).background;}
	
	double getIntegralDeviation(size_t index) const {return positionsVector.at(index).integralDeviation;}
	double getXWidthDeviation(size_t index) const {return positionsVector.at(index).widthDeviation;}
	double getYWidthDeviation(size_t index) const {return positionsVector.at(index).widthDeviation;}
	double getRotationAngleDeviation(size_t index) const {return 0;}
	double getXPositionDeviation(size_t index) const {return positionsVector.at(index).xPositionDeviation;}
	double getYPositionDeviation(size_t index) const {return positionsVector.at(index).yPositionDeviation;}
	double getZPositionDeviation(size_t index) const {return 0;}
	double getBackgroundDeviation(size_t index) const {return positionsVector.at(index).backgroundDeviation;}
	
	// adding new positions
	void addPosition(LocalizedPosition_2DGauss& newPosition) {positionsVector.push_back(newPosition);}
	void addPositions(LocalizedPositionsContainer_2DGauss& newPositionsContainer);
	
	waveHndl writePositionsToWave(std::string waveName) const {;}
	void writePositionsToFile(std::string filePath) const {;}
	
	void sortPositionsByFrameNumber() {std::sort(positionsVector.begin(), positionsVector.end(), sortCompareFrameNumber);}
	static int sortCompareFrameNumber(LocalizedPosition_2DGauss left, LocalizedPosition_2DGauss right) {
		return ((left.frameNumber <= right.frameNumber) ? 1 : 0);}
protected:
	std::vector<LocalizedPosition_2DGauss> positionsVector;
};

class LocalizedPositionsContainer_2DGaussFixedWidth : public LocalizedPositionsContainer {
public:
	LocalizedPositionsContainer_2DGaussFixedWidth() {;}
	LocalizedPositionsContainer_2DGaussFixedWidth(waveHndl wave);
	LocalizedPositionsContainer_2DGaussFixedWidth(const std::string& filePath) {throw std::runtime_error("Loading positions from files is not yet supported");}
	
	~LocalizedPositionsContainer_2DGaussFixedWidth() {;}
	
	// accessor methods
	size_t getNPositions() const {return positionsVector.size();}
	size_t getFrameNumber(size_t index) const {return positionsVector.at(index).frameNumber;}
	double getIntegral(size_t index) const {return positionsVector.at(index).integral;}
	double getXWidth(size_t index) const {return positionsVector.at(index).width;}
	double getYWidth(size_t index) const {return positionsVector.at(index).width;}
	double getRotationAngle(size_t index) const {return 0;}
	double getXPosition(size_t index) const {return positionsVector.at(index).xPosition;}
	double getYPosition(size_t index) const {return positionsVector.at(index).yPosition;}
	double getZPosition(size_t index) const {return 0;}
	double getBackground(size_t index) const {return positionsVector.at(index).background;}
	
	double getIntegralDeviation(size_t index) const {return positionsVector.at(index).integralDeviation;}
	double getXWidthDeviation(size_t index) const {return 0;}
	double getYWidthDeviation(size_t index) const {return 0;}
	double getRotationAngleDeviation(size_t index) const {return 0;}
	double getXPositionDeviation(size_t index) const {return positionsVector.at(index).xPositionDeviation;}
	double getYPositionDeviation(size_t index) const {return positionsVector.at(index).yPositionDeviation;}
	double getZPositionDeviation(size_t index) const {return 0;}
	double getBackgroundDeviation(size_t index) const {return positionsVector.at(index).backgroundDeviation;}
	
	// adding new positions
	void addPosition(LocalizedPosition_2DGaussFixedWidth& newPosition) {positionsVector.push_back(newPosition);}
	void addPositions(LocalizedPositionsContainer_2DGaussFixedWidth& newPositionsContainer);
	
	waveHndl writePositionsToWave(std::string waveName) const {;}
	void writePositionsToFile(std::string filePath) const {;}
	
	void sortPositionsByFrameNumber() {std::sort(positionsVector.begin(), positionsVector.end(), sortCompareFrameNumber);}
	static int sortCompareFrameNumber(LocalizedPosition_2DGaussFixedWidth left, LocalizedPosition_2DGaussFixedWidth right) {
		return ((left.frameNumber <= right.frameNumber) ? 1 : 0);}
protected:
	std::vector<LocalizedPosition_2DGaussFixedWidth> positionsVector;
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
	PALMAnalysisController(boost::shared_ptr<ImageLoader> imageLoader_rhs, boost::shared_ptr<ThresholdImage> thresholder_rhs,
						   boost::shared_ptr<ThresholdImage_Preprocessor> thresholdImagePreprocessor_rhs,
						   boost::shared_ptr<ThresholdImage_Postprocessor> thresholdImagePostprocessor_rhs,
						   boost::shared_ptr<ParticleFinder> particleFinder_rhs, boost::shared_ptr<FitPositions> fitPositions_rhs,
						   boost::shared_ptr<PALMResultsWriter> resultsWriter_rhs, boost::shared_ptr<PALMAnalysisProgressReporter> progressReporter_rhs);
	~PALMAnalysisController() {;}
	
	void DoPALMAnalysis();
	
protected:
	size_t nImages;
	
	std::queue<size_t> framesToBeProcessed;
	std::list<boost::shared_ptr<PALMResults> > fittedPositionsList;
	
	friend void ThreadPoolWorker(PALMAnalysisController* controller);
	
	boost::shared_ptr<ImageLoader> imageLoader;
	boost::shared_ptr<ThresholdImage> thresholder;
	boost::shared_ptr<ThresholdImage_Preprocessor> thresholdImagePreprocessor;
	boost::shared_ptr<ThresholdImage_Postprocessor> thresholdImagePostprocessor;
	boost::shared_ptr<ParticleFinder> particleFinder;
	boost::shared_ptr<FitPositions> fitPositions;
	boost::shared_ptr<PALMResultsWriter> resultsWriter;
	boost::shared_ptr<PALMAnalysisProgressReporter> progressReporter;
	
	boost::mutex acquireFrameForProcessingMutex;
	boost::mutex addPALMResultsMutex;
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
 * Compare two PALMResults to see which one originated from an earlier frame
 */
int ComparePALMResults(boost::shared_ptr<PALMResults> result1, boost::shared_ptr<PALMResults> result2);

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
	PALMAnalysisProgressReporter_Silent() {previousPercentage = 0;}
	~PALMAnalysisProgressReporter_Silent() {;}
	
	void CalculationStarted() {;}
	void UpdateCalculationProgress(double percentDone) {;}
	void CalculationDone() {;}
	void CalculationAborted() {;}
	
protected:
	double previousPercentage;
};

/**
 * @brief Print updates on the calculation progress to the Igor command line
 */
class PALMAnalysisProgressReporter_IgorCommandLine : public PALMAnalysisProgressReporter {
public:
	PALMAnalysisProgressReporter_IgorCommandLine() {previousPercentage = 0;}
	~PALMAnalysisProgressReporter_IgorCommandLine() {;}
	
	void CalculationStarted() {XOPNotice("Running PALM analysis... ");}
	void UpdateCalculationProgress(double percentDone);
	void CalculationDone() {XOPNotice("Calculation finished!\r");}
	void CalculationAborted() {XOPNotice("Abort requested by user\r");}
	
protected:
	double previousPercentage;
};

class RipleysKFunctionCalculator {
public:
	RipleysKFunctionCalculator() {;}
	~RipleysKFunctionCalculator() {;}
	
	boost::shared_ptr<vector<double> > CalculateRipleysKFunction(boost::shared_ptr<PALMMatrix<double> > positions, double startBin, double endBin, double binWidth);
};

int construct_summed_intensity_trace(ImageLoader *image_loader, string output_wave_name, long startX, long startY, long endX, long endY);

int construct_average_image(ImageLoader *image_loader, string output_wave_name, long startX, long startY, long endX, long endY);

void calculateStandardDeviationImage(ImageLoader *image_loader, string output_wave_name, long startX, long startY, long endX, long endY);

#endif
