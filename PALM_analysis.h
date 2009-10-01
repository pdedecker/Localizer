#ifndef PALM_ANALYSIS_H
#define PALM_ANALYSIS_H

#include <iostream>
#include <fstream>
#include <vector>
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
class PALMLocalizationResult {
public:
	PALMLocalizationResult() {amplitude = 0; width = 0; xPos = 0; yPos = 0; offset = 0; amplitudeError = 0, widthError = 0,
		xPosError = 0; yPosError = 0; offsetError = 0; nIterations = 0;}
	~PALMLocalizationResult() {;}
	
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
	PALMResults(size_t frameIndex_rhs, boost::shared_ptr<std::vector<PALMLocalizationResult> > fittedPositions_rhs) {
		frameIndex = frameIndex_rhs;
		fittedPositions = fittedPositions_rhs;
	}
	~PALMResults() {;}
	
	size_t getFrameIndex() {return frameIndex;}
	size_t getNumberOfFittedPositions() {return fittedPositions->size();}
	boost::shared_ptr<std::vector<PALMLocalizationResult> > getFittedPositions() {return fittedPositions;}
	
protected:
	size_t frameIndex;
	boost::shared_ptr<std::vector<PALMLocalizationResult> > fittedPositions;
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
	boost::mutex loadImagesMutex;
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

int construct_summed_intensity_trace(ImageLoader *image_loader, string output_wave_name, long startX, long startY, long endX, long endY);

int construct_average_image(ImageLoader *image_loader, string output_wave_name, long startX, long startY, long endX, long endY);

void calculateStandardDeviationImage(ImageLoader *image_loader, string output_wave_name, long startX, long startY, long endX, long endY);

#endif
