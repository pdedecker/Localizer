#ifndef PALM_ANALYSIS_H
#define PALM_ANALYSIS_H

#include <vector>
#include <algorithm>
#include "boost/thread.hpp"
#include "boost/bind.hpp"
#include <Eigen/Eigen>
#include "PALM_analysis_errors.h"
#include "PALM_analysis_FileIO.h"
#include "PALM_analysis_defines.h"
#include "PALM_analysis_segmentation.h"
#include "PALM_analysis_ParticleFinding.h"
#include "PALM_analysis_Localization.h"
#include "PALM_analysis_storage.h"
#include "PALM_analysis_ProgressReporting.h"
#include "boost/date_time/posix_time/posix_time.hpp"

#ifdef WITH_IGOR
#include "XOPStandardHeaders.h"
#endif

#define GSL_RANGE_CHECK_OFF	// this is not required since Eigen::MatrixXddoes range checks

class ImageLoader;
class ThresholdImage;
class ThresholdImage_Preprocessor;
class ThresholdImage_Postprocessor;
class ParticleFinder;
class FitPositions;


boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > do_processing_and_thresholding(boost::shared_ptr<Eigen::MatrixXd> image, boost::shared_ptr<ThresholdImage_Preprocessor>preprocessor, 
																		 boost::shared_ptr<ThresholdImage> thresholder, boost::shared_ptr<ThresholdImage_Postprocessor> postprocessor);



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
						   boost::shared_ptr<ParticleFinder> particleFinder_rhs, 
						   std::vector<boost::shared_ptr<ParticleVerifier> > particleVerifiers_rhs,
						   boost::shared_ptr<FitPositions> fitPositions_rhs,
						   boost::shared_ptr<ProgressReporter> progressReporter_rhs,
						   size_t firstFrame = (size_t)-1, size_t lastFrame = (size_t)-1);
	~PALMAnalysisController() {;}
	
	boost::shared_ptr<LocalizedPositionsContainer> DoPALMAnalysis(boost::shared_ptr<ImageLoader> imageLoader_rhs);
	// runs a PALM analysis according to the parameters passed in as objects in the constructor
	// on the data file represented by the imageloader
	// returns a LocalizedPositionsContainer containing fitted PALM positions
	
protected:
	size_t nImages;
	size_t firstFrameToAnalyze;
	size_t lastFrameToAnalyze;
	
	std::queue<size_t> framesToBeProcessed;
	boost::shared_ptr<LocalizedPositionsContainer> localizedPositions;
	
	friend void ThreadPoolWorker(PALMAnalysisController* controller);
	
	boost::shared_ptr<ImageLoader> imageLoader;
	boost::shared_ptr<ThresholdImage> thresholder;
	boost::shared_ptr<ThresholdImage_Preprocessor> thresholdImagePreprocessor;
	boost::shared_ptr<ThresholdImage_Postprocessor> thresholdImagePostprocessor;
	boost::shared_ptr<ParticleFinder> particleFinder;
	std::vector<boost::shared_ptr<ParticleVerifier> > particleVerifiers;
	boost::shared_ptr<FitPositions> fitPositions;
	boost::shared_ptr<ProgressReporter> progressReporter;
	
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
 * @brief A higher-level positions fitter that implements deflation, but needs another FitPositions to do the actual fitting work
 *
 * Deflation is the idea that fitted positions can be effectively subtracted from the image to reveal other underlying emitters that were obscured
 * at first. The validity of this approach is open to debate, however, some papers in the literature make use of this approach.
 * This class does no fitting by itself, it merely takes another FitPositions as well as all the classes required for segmentation.
 * The process is iterative: the class will keep running localize-deflate circles until no more positions are recovered.
 * The returned positions will contain all positions localized, both originally and after deflation.
 */
class FitPositionsDeflate : public FitPositions {
public:
	FitPositionsDeflate(boost::shared_ptr <ThresholdImage_Preprocessor> preprocessor_rhs, boost::shared_ptr <ThresholdImage_Postprocessor> postprocessor_rhs,
						boost::shared_ptr<ThresholdImage> thresholder_rhs, boost::shared_ptr<ParticleFinder> particleFinder_rhs, 
						boost::shared_ptr<FitPositions> positionsFitter_rhs) {
		preprocessor = preprocessor_rhs; postprocessor = postprocessor_rhs; thresholder = thresholder_rhs; positionsFitter = positionsFitter_rhs; positionsFitter = positionsFitter_rhs;}
	
	~FitPositionsDeflate() {;}
	
	boost::shared_ptr<LocalizedPositionsContainer> fit_positions(const boost::shared_ptr<Eigen::MatrixXd> image, boost::shared_ptr<std::list<position> > positions);
	
protected:
	boost::shared_ptr<Eigen::MatrixXd> subtractLocalizedPositions(boost::shared_ptr<Eigen::MatrixXd> image, boost::shared_ptr<LocalizedPositionsContainer> positions);
	
	boost::shared_ptr <ThresholdImage_Preprocessor> preprocessor;
	boost::shared_ptr <ThresholdImage_Postprocessor> postprocessor;
	boost::shared_ptr<ThresholdImage> thresholder;
	boost::shared_ptr<ParticleFinder> particleFinder;
	boost::shared_ptr<FitPositions> positionsFitter;
};

#endif
