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
 *@brief Stores the result from fitting a particular frame, remembering both the frame index as well as the fitted positions
 */
class PALMResults {
public:
	PALMResults(size_t frameIndex_rhs, boost::shared_ptr<PALMMatrix<double> > fittedPositions_rhs) {
		frameIndex = frameIndex_rhs;
		fittedPositions = fittedPositions_rhs;
	}
	~PALMResults() {;}
	
	size_t getFrameIndex() {return frameIndex;}
	size_t getNumberOfFittedPositions() {return fittedPositions->getXSize();}
	boost::shared_ptr<PALMMatrix<double> > getFittedPositions() {return fittedPositions;}
	
protected:
	size_t frameIndex;
	boost::shared_ptr<PALMMatrix<double> > fittedPositions;
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

class END_SHOULD_BE_LARGER_THAN_START{};
class GET_NTH_IMAGE_RETURNED_NULL{};


class PALMBitmapImageDeviationCalculator {
public:
	PALMBitmapImageDeviationCalculator() {;}
	virtual ~PALMBitmapImageDeviationCalculator() {;}
	
	virtual double getDeviation(boost::shared_ptr<PALMMatrix<double> > positions, size_t index) = 0;
};

class PALMBitmapImageDeviationCalculator_FitUncertainty : public PALMBitmapImageDeviationCalculator {
public:
	PALMBitmapImageDeviationCalculator_FitUncertainty(double scaleFactor_rhs, double upperLimit_rhs) {scaleFactor = scaleFactor_rhs; upperLimit = upperLimit_rhs;}
	~PALMBitmapImageDeviationCalculator_FitUncertainty() {;}
	
	double getDeviation(boost::shared_ptr<PALMMatrix<double> > positions, size_t index) {return ((positions->get(index, 7) < upperLimit) ? (positions->get(index, 7) * scaleFactor) : -1);}
	
private:
	double scaleFactor;
	double upperLimit;
};

class PALMBitmapImageDeviationCalculator_Constant : public PALMBitmapImageDeviationCalculator {
public:
	PALMBitmapImageDeviationCalculator_Constant(double deviation_rhs) {deviation = deviation_rhs;}
	PALMBitmapImageDeviationCalculator_Constant() {;}
	
	double getDeviation(boost::shared_ptr<PALMMatrix<double> > positions, size_t index) {return deviation;}
	
private:
	double deviation;
};

class PALMBitmapImageDeviationCalculator_AmplitudeSquareRoot : public PALMBitmapImageDeviationCalculator {
public:
	PALMBitmapImageDeviationCalculator_AmplitudeSquareRoot(double PSFWidth_rhs, double scaleFactor_rhs) {PSFWidth = PSFWidth_rhs; scaleFactor = scaleFactor_rhs;}
	PALMBitmapImageDeviationCalculator_AmplitudeSquareRoot() {;}
	
	double getDeviation(boost::shared_ptr<PALMMatrix<double> > positions, size_t index) {return PSFWidth / (scaleFactor * sqrt(positions->get(index, 1)));}
	
private:
	double PSFWidth;
	double scaleFactor;
};

class calculate_PALM_bitmap_image_ThreadStartParameters {
public:
	calculate_PALM_bitmap_image_ThreadStartParameters() {;}
	~calculate_PALM_bitmap_image_ThreadStartParameters() {;}
	
	boost::shared_ptr<PALMMatrix<double> > positions;
	boost::shared_ptr<PALMVolume <unsigned short> > image;
	boost::shared_ptr<PALMMatrix<double> > totalIntensities;
	
	boost::shared_ptr<PALMMatrix<double> > colors;
	boost::shared_ptr<PALMBitmapImageDeviationCalculator> deviationCalculator;
	
	int normalizeColors;
	size_t nFrames;
	size_t startIndex;
	size_t endIndex;
	size_t imageWidth;
	size_t imageHeight;
	size_t xSize;
	size_t ySize;
	double maxAmplitude;	// maximum amplitude of a fitted Gaussian over the entire fitted positions
	double scaleFactor;
};

boost::shared_ptr<PALMVolume <unsigned short> > calculate_PALM_bitmap_image(boost::shared_ptr<PALMMatrix<double> > positions, boost::shared_ptr<PALMMatrix<double> > colors, boost::shared_ptr<PALMBitmapImageDeviationCalculator> deviationCalculator,
																size_t xSize, size_t ySize, size_t imageWidth, size_t imageHeight, int normalizeColors);

boost::shared_ptr<PALMVolume <unsigned short> > calculate_PALM_bitmap_image_parallel(boost::shared_ptr<PALMMatrix<double> > positions, boost::shared_ptr<PALMMatrix<double> > colors, boost::shared_ptr<PALMBitmapImageDeviationCalculator> deviationCalculator,
																		 size_t xSize, size_t ySize, size_t imageWidth, size_t imageHeight, int normalizeColors);

void calculate_PALM_bitmap_image_ThreadStart(boost::shared_ptr<calculate_PALM_bitmap_image_ThreadStartParameters> startParameters);


int construct_summed_intensity_trace(ImageLoader *image_loader, string output_wave_name, long startX, long startY, long endX, long endY);

int construct_average_image(ImageLoader *image_loader, string output_wave_name, long startX, long startY, long endX, long endY);

void calculateStandardDeviationImage(ImageLoader *image_loader, string output_wave_name, long startX, long startY, long endX, long endY);

#endif
