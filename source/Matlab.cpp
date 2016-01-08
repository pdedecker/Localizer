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
 with analysis programs such as Igor Pro or Matlab, 
 the licensors of this Program grant you additional permission 
 to convey the resulting work.
 */

#include "Matlab.h"

#include "Defines.h"
#include "PALMAnalysis.h"
#include "ParticleFinding.h"
#include "FileIO.h"
#include "SOFI.h"
#include "NewSOFI.h"

/**
 The Matlab interface exports just a single 'gateway' function, mexFunction. This is the only
 function that can be called from Matlab. By convention the first argument to this function should
 be a string that identifies the operation to be performed. Depending on the value of this string
 the rhs and lhs parameters passed here are applied to the more specific functions.
*/
void mexFunction(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs) {
    gsl_set_error_handler_off();	// we will handle GSL errors ourselves
    TIFFSetErrorHandler(NULL);      // we will handle libtiff errors ourselves
	
	if (nrhs < 1)
		mexErrMsgTxt("First argument must be a string describing the operation to perform (\"readccdimages\", \"writeccdimages\", \"localize\", \"testsegmentation\", \"sofi\", or \"newsofi\")");
	
	const mxArray* inputArray = prhs[0];
	if ((mxGetClassID(inputArray) != mxCHAR_CLASS) || (mxGetM(inputArray) > 1))
		mexErrMsgTxt("First argument must be a string describing the operation to perform (\"readccdimages\", \"writeccdimages\", \"localize\", \"testsegmentation\", \"sofi\", or \"newsofi\")");

	std::string selector = GetMatlabString(inputArray);
	if (boost::iequals(selector, "localize")) {
		MatlabLocalization(nlhs, plhs, nrhs, prhs);
	} else if (boost::iequals(selector, "testsegmentation")) {
		MatlabTestSegmentation(nlhs, plhs, nrhs, prhs);
	} else if (boost::iequals(selector, "sofi")) {
		MatlabSOFI(nlhs, plhs, nrhs, prhs);
	} else if (boost::iequals(selector, "newsofi")) {
		MatlabNewSOFI(nlhs, plhs, nrhs, prhs);
	} else if (boost::iequals(selector, "readccdimages")) {
		MatlabReadCCDImages(nlhs, plhs, nrhs, prhs);
    } else if (boost::iequals(selector, "writeccdimages")) {
        MatlabWriteCCDImages(nlhs, plhs, nrhs, prhs);
	} else {
		mexErrMsgTxt("Unknown selector (should be one of \"readccdimages\", \"writeccdimages\", \"localize\", \"testsegmentation\", \"sofi\", or \"newsofi\")");
	}
}

/**
 Function that will handle localization. The input arguments must be of the form
 0 - the string "localize"
 1 - the expect standard deviation of the psf
 2 - segmentation algorithm, 'glrt' or 'smoothsigma'
 3 - the PFA if GLRT, otherwise smoothsigma factor
 4 - localization algorithm. '2DGauss' only for now
 5 - file path or a matrix containing numeric data
*/
void MatlabLocalization(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs) {
	// input validation
	if (nrhs != 6)
		mexErrMsgTxt("Must have exactly 6 input arguments for localize\n\
					 1 - the string \"localize\"\n\
					 2. the estimated standard deviation of the PSF (in pixels)\n\
					 3. the segmentation algorithm (\"glrt\" or \"smoothsigma\"\n\
					 4. the PFA if GLRT, otherwise smoothsigma factor\n\
					 5. the localization algorithm. '2DGauss', '2DGaussFixedWidth', 'IterativeMultiplication', 'Centroid', 'Ellipsoidal2DGauss', 'Ellipsoidal2DGaussAstigmatism', or 'MLEwG'\n\
					 6. string containing the path to the file, or a 2D or 3D matrix containing image data\n\
					 \n\n\
					 The format of the returned positions depends on the localization algorithm. See localize.m for the details.");

	if (nlhs != 1)
		mexErrMsgTxt("Must have exactly one left hand side argument for localize");
	
	const mxArray* array;
	// the mxArray at index 0 will have been checked already by mexFunction
	
	// index 1 - must be a single number
	array = prhs[1];
	if ((mxGetN(array) != 1) || (mxGetM(array) != 1) || (mxGetClassID(array) != mxDOUBLE_CLASS))
		mexErrMsgTxt("2nd argument must be a scalar of type double (the PSF standard deviation)");
	double psfWidth = *(mxGetPr(array));
	if (psfWidth <= 0.0)
		mexErrMsgTxt("Expected positive non-zero number for the PSF width");
	
	// index 2 - must be a string, either 'GLRT' or 'SmoothSigma' (not case sensitive)
	array = prhs[2];
	if (mxGetClassID(array) != mxCHAR_CLASS)
		mexErrMsgTxt("3rd argument must be a string (the segmentation algorithm)");
	std::string segmentationStr = GetMatlabString(array);
	
	int segmentationAlgorithm;
	if (boost::iequals(segmentationStr, "GLRT")) {
		segmentationAlgorithm = THRESHOLD_METHOD_GLRT;
	} else if (boost::iequals(segmentationStr, "SmoothSigma")) {
		segmentationAlgorithm = THRESHOLD_METHOD_SMOOTHSIGMA;
	} else {
		mexErrMsgTxt("Unknown segmentation algorithm");
	}
	
	// index 3 - must be a single number
	array = prhs[3];
	if ((mxGetN(array) != 1) || (mxGetM(array) != 1) || (mxGetClassID(array) != mxDOUBLE_CLASS))
		mexErrMsgTxt("4th argument must be a scalar of type double (the segmentation parameter)");
	double thresholdParam = *(mxGetPr(array));
	if (thresholdParam <= 0.0)
		mexErrMsgTxt("Expected positive non-zero number for the segmentation parameter");

	// index 4 - must be a string
	array = prhs[4];
	if (mxGetClassID(array) != mxCHAR_CLASS)
		mexErrMsgTxt("5th argument must be a string (the localization algorithm)");
	std::string localizationStr = GetMatlabString(array);

	int localizationAlgorithm;
	if (boost::iequals(localizationStr, "2DGauss")) {
		localizationAlgorithm = LOCALIZATION_METHOD_2DGAUSS;
	} else if (boost::iequals(localizationStr, "2DGaussFixedWidth")) {
		localizationAlgorithm = LOCALIZATION_METHOD_2DGAUSS_FIXEDWIDTH;
	} else if (boost::iequals(localizationStr, "IterativeMultiplication")) {
		localizationAlgorithm = LOCALIZATION_METHOD_MULTIPLICATION;
	} else if (boost::iequals(localizationStr, "Centroid")) {
		localizationAlgorithm = LOCALIZATION_METHOD_CENTROID;
	} else if (boost::iequals(localizationStr, "Ellipsoidal2DGauss")) {
		localizationAlgorithm = LOCALIZATION_METHOD_2DGAUSS_ELLIPSOIDAL;
	} else if (boost::iequals(localizationStr, "Ellipsoidal2DGaussAstigmatism")) {
		localizationAlgorithm = LOCALIZATION_METHOD_2DGAUSS_ELLIPSOIDAL_ASTIGMATISM;
	} else if (boost::iequals(localizationStr, "MLEwG")) {
		localizationAlgorithm = LOCALIZATION_METHOD_MLEwG;
	} else {
		mexErrMsgTxt("Unknown localization algorithm");
	}
	
	std::string filePath;
	mxArray* dataArray = NULL;
	// index 5 - must be the data
	array = prhs[5];
	if (mxGetClassID(array) == mxCHAR_CLASS) {
		filePath = GetMatlabString(array);
		if (filePath.empty())
			mexErrMsgTxt("Need a non-empty filepath to the data");
	} else {
		// assume that this is a valid data array
		// if it isn't then the image loader will throw an error
		dataArray = const_cast<mxArray*>(array);
	}
	
	// now run the actual localization
	try {
		std::shared_ptr<ImageLoader> imageLoader;
		if (!filePath.empty()) {
			imageLoader = GetImageLoader(filePath);
		} else {
			imageLoader = std::shared_ptr<ImageLoader>(new ImageLoaderMatlab(dataArray));
		}
		
		std::shared_ptr<ThresholdImage> thresholder;
		switch(segmentationAlgorithm) {
			case THRESHOLD_METHOD_GLRT:
				thresholder = std::shared_ptr<ThresholdImage>(new ThresholdImage_GLRT_FFT(thresholdParam, psfWidth));
				break;
			case THRESHOLD_METHOD_SMOOTHSIGMA:
				thresholder = std::shared_ptr<ThresholdImage>(new ThresholdImage_SmoothSigma(psfWidth, thresholdParam));
				break;
			default:
				throw std::runtime_error("Unknown segmentation method");
				break;
        }
		
		std::shared_ptr<ThresholdImage_Preprocessor> preprocessor(new ThresholdImage_Preprocessor_DoNothing());
		std::shared_ptr<ThresholdImage_Postprocessor> postprocessor(new ThresholdImage_Postprocessor_DoNothing());
		
		std::shared_ptr<ParticleFinder> particleFinder(new ParticleFinder_adjacent8());
		
		std::vector<std::shared_ptr<ParticleVerifier> > particleVerifiers;
		particleVerifiers.push_back(std::shared_ptr<ParticleVerifier> (new ParticleVerifier_RemoveOverlappingParticles(psfWidth)));
		
		std::shared_ptr<FitPositions> positionsFitter;
		double sigma = 1.0;
		switch (localizationAlgorithm) {
            case LOCALIZATION_METHOD_2DGAUSS:
                positionsFitter = std::shared_ptr<FitPositions>(new FitPositions_SymmetricGaussian(psfWidth, sigma));
                break;
            case LOCALIZATION_METHOD_2DGAUSS_FIXEDWIDTH:
                positionsFitter = std::shared_ptr<FitPositions>(new FitPositions_FixedWidthGaussian(psfWidth, sigma));
                break;
            case LOCALIZATION_METHOD_MULTIPLICATION:
                positionsFitter = std::shared_ptr<FitPositions>(new FitPositionsMultiplication(psfWidth, sigma));
                break;
            case LOCALIZATION_METHOD_CENTROID:
                positionsFitter = std::shared_ptr<FitPositions>(new FitPositionsCentroid(psfWidth));
                break;
            case LOCALIZATION_METHOD_2DGAUSS_ELLIPSOIDAL:
                positionsFitter = std::shared_ptr<FitPositions>(new FitPositions_EllipsoidalGaussian_SymmetricPSF(psfWidth, sigma));
                break;
			case LOCALIZATION_METHOD_2DGAUSS_ELLIPSOIDAL_ASTIGMATISM:
				positionsFitter = std::shared_ptr<FitPositions>(new FitPositions_EllipsoidalGaussian(psfWidth, sigma));
				break;
            case LOCALIZATION_METHOD_MLEwG:
                positionsFitter = std::shared_ptr<FitPositions>(new FitPositions_MLEwG(psfWidth));
                break;
            default:
                throw std::runtime_error("Unknown localization method");
                break;
		}
		
		std::shared_ptr<ProgressReporter> progressReporter(new ProgressReporter_MatlabWaitMex());
		
		std::shared_ptr<PALMAnalysisController> analysisController(new PALMAnalysisController(thresholder, preprocessor,
																								postprocessor, particleFinder, particleVerifiers,
																								positionsFitter,
																								progressReporter, -1,
																								-1));
		
		// do the actual calculation
        ImagePtr localizedPositions = analysisController->DoPALMAnalysis(imageLoader)->getLocalizedPositionsAsMatrix();
		mxArray* outputArray = mxCreateDoubleMatrix(localizedPositions->rows(), localizedPositions->cols(), mxREAL);
		if (outputArray == NULL)
			throw std::bad_alloc();
		double* positionsPtr = localizedPositions->data();
		double* outputPositionsPtr = mxGetPr(outputArray);
		memcpy(outputPositionsPtr, positionsPtr, localizedPositions->rows() * localizedPositions->cols() * sizeof(double));
		plhs[0] = outputArray;
	}
	catch (std::bad_alloc) {
        mexErrMsgTxt("Insufficient memory");
    }
    catch (int e) {
        mexErrMsgTxt("Int error");
    }
    catch (USER_ABORTED e) {
        mexErrMsgTxt("User abort");
    }
    catch (std::runtime_error e) {
        mexErrMsgTxt(e.what());
    }
    catch (...) {
        mexErrMsgTxt("Unknown error");
    }
}

/**
 Function that will handle segmentation testing. The input arguments must be of the form
 0 - the string "testsegmentation"
 1 - the expect standard deviation of the psf
 2 - segmentation algorithm, 'glrt' or 'smoothsigma'
 3 - the PFA if GLRT, otherwise smoothsigma factor
 4 - mxArray containing the image to segment
 */
void MatlabTestSegmentation(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs) {
	// input validation
	if (nrhs != 5)
		mexErrMsgTxt("Must have exactly 5 input arguments for testsegmentation\n\
					 0 - the string \"testsegmentation\"\n\
					 1 - the expect standard deviation of the psf\n\
					 2 - segmentation algorithm, \"glrt\" or \"smoothsigma\"\n\
					 3 - the PFA if GLRT, otherwise smoothsigma factor\n\
					 4 - 2D matrix containing the image to segment");

	if (nlhs != 2)
		mexErrMsgTxt("Must have exactly two left hand side arguments for testsegmentation: the\
					 binary segmented image and list of estimated particles (not localized!)");
	
	const mxArray* array;
	// the mxArray at index 0 will have been checked already by mexFunction
	
	// index 1 - must be a single number
	array = prhs[1];
	if ((mxGetN(array) != 1) || (mxGetM(array) != 1) || (mxGetClassID(array) != mxDOUBLE_CLASS))
		mexErrMsgTxt("2nd argument must be a scalar of type double (the PSF standard deviation)");
	double psfWidth = *(mxGetPr(array));
	if (psfWidth <= 0.0)
		mexErrMsgTxt("Expected positive non-zero number for the PSF width");
	
	// index 2 - must be a string, either 'GLRT' or 'SmoothSigma' (not case sensitive)
	array = prhs[2];
	if (mxGetClassID(array) != mxCHAR_CLASS)
		mexErrMsgTxt("3rd argument must be a string (the segmentation algorithm, 'GLRT' or 'SmoothSigma')");
	std::string segmentationStr = GetMatlabString(array);
	
	int segmentationAlgorithm;
	if (boost::iequals(segmentationStr, "GLRT")) {
		segmentationAlgorithm = THRESHOLD_METHOD_GLRT;
	} else if (boost::iequals(segmentationStr, "SmoothSigma")) {
		segmentationAlgorithm = THRESHOLD_METHOD_SMOOTHSIGMA;
	} else {
		mexErrMsgTxt("Unknown segmentation algorithm");
	}
	
	// index 3 - must be a single number
	array = prhs[3];
	if ((mxGetN(array) != 1) || (mxGetM(array) != 1) || (mxGetClassID(array) != mxDOUBLE_CLASS))
		mexErrMsgTxt("4th argument must be a scalar of type double (the segmentation parameter)");
	double segmentationParameter = *(mxGetPr(array));
	if (segmentationParameter <= 0.0)
		mexErrMsgTxt("Expected positive non-zero number for the segmentation parameter");
	
	// index 4 - must contain the image to segment
	array = prhs[4];
	if ((mxGetN(array) <= 1) || (mxGetM(array) <= 1))
		mexErrMsgTxt("5th argument must be a 2D matrix (the image to segment)");
	
	try {
		boost::scoped_ptr<ImageLoader> imageLoader(new ImageLoaderMatlab(const_cast<mxArray*>(array)));
		ImagePtr imageToSegment = imageLoader->readImage(0);
		
		std::shared_ptr<ThresholdImage> thresholder;
		switch(segmentationAlgorithm) {
			case THRESHOLD_METHOD_GLRT:	// the GLRT test proposed by Arnauld et al in Nat Methods 5:687 2008
				thresholder = std::shared_ptr<ThresholdImage>(new ThresholdImage_GLRT_FFT(segmentationParameter, psfWidth));
				break;
			case THRESHOLD_METHOD_SMOOTHSIGMA:
				thresholder = std::shared_ptr<ThresholdImage>(new ThresholdImage_SmoothSigma(psfWidth, segmentationParameter));
				break;
			default:
				throw std::runtime_error("Unknown segmentation method");
				break;
        }
		
		std::shared_ptr<ThresholdImage_Preprocessor> preprocessor(new ThresholdImage_Preprocessor_DoNothing());
		std::shared_ptr<ThresholdImage_Postprocessor> postprocessor(new ThresholdImage_Postprocessor_DoNothing());
		std::shared_ptr<ParticleFinder> particlefinder(new ParticleFinder_adjacent8());
		
		// calculate the threshold
        std::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > thresholded_image = do_processing_and_thresholding(imageToSegment, preprocessor, thresholder, postprocessor);
		std::shared_ptr<std::list<Particle> > located_particles = particlefinder->findPositions(imageToSegment, thresholded_image);
		
		// and copy the results back out
		size_t segmentedRows = thresholded_image->rows();
		size_t segmentedCols = thresholded_image->cols();
		size_t nSegmentedPoints = segmentedRows * segmentedCols;
		
		// store the segmented image
		mwSize ndims = 2;
		mwSize dims[2];
		dims[0] = segmentedRows;
		dims[1] = segmentedCols;
		mxArray* segmentedOutputMatrix = mxCreateNumericArray(ndims, dims, mxUINT8_CLASS, mxREAL);
		if (segmentedOutputMatrix == NULL)
			throw std::bad_alloc();
		
		uint8_t* segmentedOutputPtr = reinterpret_cast<uint8_t*>(mxGetPr(segmentedOutputMatrix));
		int* segmentedImagePtr = thresholded_image->data();
		for (int i = 0; i < nSegmentedPoints; ++i) {
			segmentedOutputPtr[i] = segmentedImagePtr[i];
		}
		
		// store the particles
		size_t nParticles = located_particles->size();
		ndims = 2;
		dims[0] = nParticles;
		dims[1] = 4;
		mxArray* particleOutputMatrix = mxCreateNumericArray(ndims, dims, mxDOUBLE_CLASS, mxREAL);
		double* particleOutputPtr = reinterpret_cast<double*>(mxGetPr(particleOutputMatrix));
		if (particleOutputMatrix == NULL)
			throw std::bad_alloc();
		
		int offset = 0;
		for (std::list<Particle>::iterator it = located_particles->begin(); it != located_particles->end(); ++it, ++offset) {
			particleOutputPtr[offset] = (*it).intensity;
			particleOutputPtr[offset + nParticles] = (*it).x;
			particleOutputPtr[offset + 2 * nParticles] = (*it).y;
			particleOutputPtr[offset + 3 * nParticles] = (*it).background;
		}
		
		plhs[0] = segmentedOutputMatrix;
		plhs[1] = particleOutputMatrix;
	}
	catch (std::bad_alloc) {
        mexErrMsgTxt("Insufficient memory");
    }
    catch (int e) {
        mexErrMsgTxt("Int error");
    }
    catch (USER_ABORTED e) {
        mexErrMsgTxt("User abort");
    }
    catch (std::runtime_error e) {
        mexErrMsgTxt(e.what());
    }
    catch (...) {
        mexErrMsgTxt("Unknown error");
    }
}

/**
 Function that will handle SOFI calculations. The input arguments must be of the form
 0 - the string "sofi"
 1 - file path to a data file, or a 2D or 3D matrix containing numeric data
 2 - 1D array containing the orders to calculate
 all other arguments must be keyword-value pairs, meaning that each consists of
 two arguments: a string containing the keyword, followed by the keyword's value.
 */
void MatlabSOFI(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs) {
	// input validation
	if (nrhs != 7)
		mexErrMsgTxt("Must have exactly 7 input arguments for sofi\n\
					 1. the string \"sofi\"\n\
					 2. the desired order (2 or 3)\n\
					 3. single number: 0 for autocumulants, 1 for crosscumulants\n\
					 4. 1D matrix containing time lags to use in the calculation. An empty matrix\n\
						is equivalent to specifying time lags of zero. For a calculation of order N\n\
						at least (N - 1) values must be specified. Excess values will be ignored.\n\
					 5. number of frames to skip at the beginning\n\
					 6. number of frames to include in the calculation (<=0 means all frames up to the end)\n\
					 7. file path to a data file, or a 2D or 3D matrix containing numeric data");
	
	if (nlhs != 2)
		mexErrMsgTxt("Must have exactly two left hand side arguments for sofi (sofi output image and average output image");
	
	mxArray* array;
	// the mxArray at index 0 (the string "sofi") will have been checked already by mexFunction
	
	// index 1 - must be a single number (order)
	array = const_cast<mxArray*>(prhs[1]);
	if ((mxGetN(array) != 1) || (mxGetM(array) != 1) || (mxGetClassID(array) != mxDOUBLE_CLASS))
		mexErrMsgTxt("2nd argument must be a double scalar (the order of the calculation)");
	int correlationOrder = *(mxGetPr(array));
	if ((correlationOrder < 2) || (correlationOrder > 3))
		mexErrMsgTxt("Expected 2 or 3 for the correlation order");
	
	// index 2 - must be a boolean (auto or cross)
	bool doCrossCorrelation;
	array = const_cast<mxArray*>(prhs[2]);
	if ((mxGetN(array) != 1) || (mxGetM(array) != 1) || (mxGetClassID(array) != mxDOUBLE_CLASS))
		mexErrMsgTxt("3nd argument must be a double scalar (zero for autocorrelation, non-zero for crosscorrelation)");
	doCrossCorrelation = (*(mxGetPr(array)) != 0.0);

	// index 3 - must be an empty or 1D matrix containing time lags
	std::vector<int> lagTimes;
	array = const_cast<mxArray*>(prhs[3]);
	if ((mxGetM(array) != 1) || (mxGetClassID(array) != mxDOUBLE_CLASS))
		mexErrMsgTxt("3nd argument must be a 1D or empty matrix containing lag times");
	int nLagTimes = mxGetN(array);
	if ((nLagTimes > 0) && (nLagTimes < correlationOrder - 1))
		mexErrMsgTxt("a calculation of order N requires at least (N-1) lag times");
	for (int i = 0; i < nLagTimes; ++i) {
		int thisLagTime = static_cast<int>(*(mxGetPr(array) + i));
		if (thisLagTime < 0)
			mexErrMsgTxt("all time lags must be >= 0");
		lagTimes.push_back(thisLagTime);
	}
	
	// index 4 - must be a single number (number of frames to skip)
	int nFramesToSkip = 0;
	array = const_cast<mxArray*>(prhs[4]);
	if ((mxGetN(array) != 1) || (mxGetM(array) != 1) || (mxGetClassID(array) != mxDOUBLE_CLASS))
		mexErrMsgTxt("4th argument must be a double scalar (number of frames to skip)");
	if (*mxGetPr(array) <= 0) {
		nFramesToSkip = 0;
	} else {
		nFramesToSkip = static_cast<int>(*mxGetPr(array) + 0.5);
	}
	
	// index 5 - must be a single number (number of frames to include)
	int nFramesToInclude = -1;
	array = const_cast<mxArray*>(prhs[5]);
	if ((mxGetN(array) != 1) || (mxGetM(array) != 1) || (mxGetClassID(array) != mxDOUBLE_CLASS))
		mexErrMsgTxt("5th argument must be a double scalar (number of frames to include)");
	if (*mxGetPr(array) <= 0) {
		nFramesToInclude = -1;
	} else {
		nFramesToInclude = static_cast<int>(*mxGetPr(array) + 0.5);
	}
	
	// index 6 - must be the data
	std::string filePath;
	mxArray* dataArray = NULL;
	array = const_cast<mxArray*>(prhs[6]);
	if (mxGetClassID(array) == mxCHAR_CLASS) {
		filePath = GetMatlabString(array);
		if (filePath.empty())
			mexErrMsgTxt("Need a non-empty filepath to the data");
	} else {
		// assume that this is a valid data array
		// if it isn't then the image loader will throw an error
		dataArray = array;
	}
	
	try {
		std::shared_ptr<ImageLoader> imageLoader;
		if (!filePath.empty()) {
			imageLoader = GetImageLoader(filePath);
		} else {
			imageLoader = std::shared_ptr<ImageLoader>(new ImageLoaderMatlab(dataArray));
		}
		
		std::shared_ptr<ProgressReporter> progressReporter(new ProgressReporter_MatlabWaitMex());
		
		// no frame verifiers for now
		std::vector<std::shared_ptr<SOFIFrameVerifier> > frameVerifiers;
		std::vector<ImagePtr> sofiOutputImages, averageOutputImages;
		DoSOFIAnalysis(imageLoader, frameVerifiers, progressReporter,nFramesToSkip, nFramesToInclude, lagTimes, correlationOrder, doCrossCorrelation, 0, sofiOutputImages, averageOutputImages);

		plhs[0] = ConvertImagesToArray(sofiOutputImages);   // was saving only the first image. Noticed by Matthieu Dumont, who also created a fix.
		plhs[1] = ConvertImagesToArray(averageOutputImages);
	}
	catch (std::bad_alloc) {
        mexErrMsgTxt("Insufficient memory");
    }
    catch (int e) {
        mexErrMsgTxt("Int error");
    }
    catch (USER_ABORTED e) {
        mexErrMsgTxt("User abort");
    }
    catch (std::runtime_error e) {
        mexErrMsgTxt(e.what());
    }
    catch (...) {
        mexErrMsgTxt("Unknown error");
    }
}

void ParseSOFIKeywordArguments(const mxArray** prhs, int nrhs, int lhs, SOFIOptions& sofiOptions, std::vector<int>& lagTimes) {
    int firstKeywordIndex = 3;
    
    // check that we have an even number of arguments remaining
    if ((nrhs < firstKeywordIndex) || (((nrhs - firstKeywordIndex) % 2) != 0)) {
        char buf[128];
        sprintf(buf, "all arguments after the first %d must be keyword-value pairs", firstKeywordIndex);
        mexErrMsgTxt(buf);
    }
    
    // check that all remaining arguments are string - value pairs
    for (int i = firstKeywordIndex; i < nrhs; i += 2) {
        const mxArray* array = prhs[i];
        if ((mxGetClassID(array) != mxCHAR_CLASS) || (mxGetM(array) > 1)) {
            char buf[128];
            sprintf(buf, "all arguments after the first %d must be keyword (string) - value (double) pairs", firstKeywordIndex);
            mexErrMsgTxt(buf);
        }
        if (!mxIsDouble(prhs[i + 1])) {
            mexErrMsgTxt("all keyword arguments (not the keywords themselves) must be of numeric type double");
        }
    }
    
    // most keyword arguments must be 1x1 matrices (scalars) of type double
    for (int i = firstKeywordIndex; i < nrhs; i += 2) {
        const std::string keyword = GetMatlabString(prhs[i]);
        const mxArray* argument = prhs[i + 1];
        
        if (!mxIsDouble(argument)) {
            mexErrMsgTxt("all keyword arguments must be of type 'double'");
        }
        
        if (boost::iequals(keyword, "lagtimes")) {
            continue;   // not 1x1 in size
        }
        
        if ((mxGetNumberOfDimensions(argument) > 2) || (mxGetM(argument) != 1) || (mxGetN(argument) != 1)) {
            mexErrMsgTxt("all keyword arguments must be 1x1 matrices (scalar values)");
        }
    }
    
    // extract the values of the keyword arguments and perform more detailed error checking
    for (int i = firstKeywordIndex; i < nrhs; i += 2) {
        const std::string keyword = GetMatlabString(prhs[i]);
        const mxArray* argument = prhs[i + 1];
        
        if (boost::iequals(keyword, "pixelationCorrection")) {
            sofiOptions.doPixelationCorrection = (*mxGetPr(argument) != 0.0);
        } else if (boost::iequals(keyword, "alsoCorrectVariance")) {
            sofiOptions.alsoCorrectVariance = (*mxGetPr(argument) != 0.0);
        } else if (boost::iequals(keyword, "batchSize")) {
            if (*mxGetPr(argument) <= 0.0) {
                mexErrMsgTxt("'batchSize' argument requires a positive non-zero value");
            }
            sofiOptions.batchSize = *mxGetPr(argument);
        } else if (boost::iequals(keyword, "pixelcombinations")) {
            if (*mxGetPr(argument) <= 0.0) {
                mexErrMsgTxt("'pixelcombinations' argument requires a postive non-zero value");
            }
        } else if (boost::iequals(keyword, "jackknife")) {
            sofiOptions.wantJackKnife = (*mxGetPr(argument) != 0.0);
        } else if (boost::iequals(keyword, "lagtimes")) {
            if ((mxGetNumberOfDimensions(argument) > 2) || (mxGetM(argument) != 1)) {
                mexErrMsgTxt("lag times must be provided as a row vector");
            }
            double* lagPtr = mxGetPr(argument);
            for (size_t i = 0; i < mxGetN(argument); ++i) {
                lagTimes.push_back(*lagPtr);
                lagPtr += 1;
            }
        }
        
    }
}

/**
 Function that will handle SOFI calculations using the "new" algorithm. The input arguments must be of the form
 0 - the string "newsofi"
 1 - the desired order
 2 - single number: 0 for autocorrelation, 1 for crosscorrelation
 3 - file path or a matrix containing numeric data
 */
void MatlabNewSOFI(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs) {
	// input validation
	if (nrhs != 4)
		mexErrMsgTxt("Must have exactly 4 input arguments for new sofi\n\
					 1. the string \"newsofi\"\n\
					 2. the desired order (up to 6)\n\
					 3. single number: 0 for no pixelation correction, 1 for pixelation correction\n\
					 4. file path to a data file, or a 2D or 3D matrix containing numeric data");
	
	if (nlhs != 2)
		mexErrMsgTxt("Must have exactly two left hand side arguments for sofi (sofi output image and average output image");
	
	mxArray* array;
	// the mxArray at index 0 (the string "newsofi") will have been checked already by mexFunction
	
	// index 1 - must be a vector containing the orders to process
	std::vector<int> orders;
	array = const_cast<mxArray*>(prhs[1]);
	if ((mxGetM(array) != 1) || (mxGetClassID(array) != mxDOUBLE_CLASS))
		mexErrMsgTxt("1st argument must be a vector of double scalars (the orders to be calculated)");
	int nOrders = mxGetN(array);
	if (nOrders < 1)
		mexErrMsgTxt("at least one order should be specified");
	for (int i = 0; i < nOrders; ++i) {
		int thisOrder = static_cast<int>(*(mxGetPr(array) + i));
		if ((thisOrder < 1) || (thisOrder > 6))
			mexErrMsgTxt("all orders must be >=2 and <=6");
		orders.push_back(thisOrder);
	}
	
	// index 2 - must be a boolean (do pixelation correction)
	bool doPixelationCorrection;
	array = const_cast<mxArray*>(prhs[2]);
	if ((mxGetN(array) != 1) || (mxGetM(array) != 1) || (mxGetClassID(array) != mxDOUBLE_CLASS))
		mexErrMsgTxt("2nd argument must be a double scalar (zero for no pixelation correction, non-zero for pixelation correction)");
	doPixelationCorrection = (*(mxGetPr(array)) != 0.0);

	// index 3 - must be the data
	std::string filePath;
	mxArray* dataArray = NULL;
	array = const_cast<mxArray*>(prhs[3]);
	if (mxGetClassID(array) == mxCHAR_CLASS) {
		filePath = GetMatlabString(array);
		if (filePath.empty())
			mexErrMsgTxt("Need a non-empty filepath to the data");
	} else {
		// assume that this is a valid data array
		// if it isn't then the image loader will throw an error
		dataArray = array;
	}
	
	try {
		std::shared_ptr<ImageLoader> imageLoader;
		if (!filePath.empty()) {
			imageLoader = GetImageLoader(filePath);
		} else {
			imageLoader = std::shared_ptr<ImageLoader>(new ImageLoaderMatlab(dataArray));
		}
		
		std::shared_ptr<ProgressReporter> progressReporter(new ProgressReporter_MatlabWaitMex());
		
		std::vector<std::shared_ptr<SOFIFrameVerifier> > frameVerifiers; // no frame verifiers for now
		std::vector<ImagePtr> sofiOutputImages;
		// set SOFI calculation options
		SOFIOptions sofiOptions;
        sofiOptions.orders = orders;
        sofiOptions.doPixelationCorrection = doPixelationCorrection;
        sofiOptions.frameVerifiers = frameVerifiers;
        sofiOptions.wantAverageImage = true;
        DoNewSOFI(imageLoader, sofiOptions, progressReporter, sofiOutputImages);

		mwSize dims[2] = {static_cast<mwSize>(1), static_cast<mwSize>(nOrders)};
		mxArray * sofiOutputCellArray = mxCreateCellArray(2, dims);
		for (int i = 0; i < nOrders; ++i) {
			mxSetCell(sofiOutputCellArray, i, ConvertImageToArray(sofiOutputImages.at(i)));
		}
		plhs[0] = sofiOutputCellArray;
		plhs[1] = ConvertImageToArray(sofiOptions.averageImage);
	}
	catch (std::bad_alloc) {
        mexErrMsgTxt("Insufficient memory");
    }
    catch (int e) {
        mexErrMsgTxt("Int error");
    }
    catch (USER_ABORTED e) {
        mexErrMsgTxt("User abort");
    }
    catch (std::runtime_error e) {
        mexErrMsgTxt(e.what());
    }
    catch (...) {
        mexErrMsgTxt("Unknown error");
    }
}

/**
 Function that will load images from disk. The input arguments must be of the form
 0 - the string "readccdimages"
 1 - the index of the first image to read (starting from 0)
 2 - the number of images to read
 3 - file path or a matrix containing numeric data
 */
void MatlabReadCCDImages(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs) {
	// input validation
	if (nrhs != 4)
		mexErrMsgTxt("Must have exactly 4 input arguments for readccdimages\n\
					 1. the string \"readccdimages\"\n\
					 2. the index of the first image to read (starting from 0)\n\
					 3. the number of images to read (or -1 for all images)\n\
					 4. file path to a data file, or a 2D or 3D matrix containing numeric data");
	
	if (nlhs != 1)
		mexErrMsgTxt("Must have exactly one left hand side argument for readccdimages");

	mxArray* array;
	// the mxArray at index 0 will have been checked already by mexFunction
	
	// index 1 - must be a single number - the index of the first image to read
	array = const_cast<mxArray*>(prhs[1]);
	if ((mxGetN(array) != 1) || (mxGetM(array) != 1) || (mxGetClassID(array) != mxDOUBLE_CLASS))
		mexErrMsgTxt("2nd argument must be a double scalar: the index of the first image to read (starting from 0)");
	int firstImageToRead = *(mxGetPr(array));
	if (firstImageToRead < 0)
		firstImageToRead = 0;
	
	// index 2 - must be a single number - the number of images to read
	array = const_cast<mxArray*>(prhs[2]);
	if ((mxGetN(array) != 1) || (mxGetM(array) != 1) || (mxGetClassID(array) != mxDOUBLE_CLASS))
		mexErrMsgTxt("3nd argument must be a double scalar: the number of images to read (or -1 for all images)");
	int nImagesToRead = *(mxGetPr(array));
	if (nImagesToRead < 0)
		nImagesToRead = -1;

	// index 3 - must be the data
	std::string filePath;
	mxArray* dataArray = NULL;
	array = const_cast<mxArray*>(prhs[3]);
	if (mxGetClassID(array) == mxCHAR_CLASS) {
		filePath = GetMatlabString(array);
		if (filePath.empty())
			mexErrMsgTxt("Need a non-empty filepath to the data");
	} else {
		// assume that this is a valid data array
		// if it isn't then the image loader will throw an error
		dataArray = array;
	}

	try {
		std::shared_ptr<ImageLoader> imageLoader;
		if (!filePath.empty()) {
			imageLoader = GetImageLoader(filePath);
		} else {
			imageLoader = std::shared_ptr<ImageLoader>(new ImageLoaderMatlab(dataArray));
		}

		plhs[0] = LoadImagesIntoArray(imageLoader, firstImageToRead, nImagesToRead);
	}
	catch (std::bad_alloc) {
        mexErrMsgTxt("Insufficient memory");
    }
    catch (int e) {
        mexErrMsgTxt("Int error");
    }
    catch (USER_ABORTED e) {
        mexErrMsgTxt("User abort");
    }
    catch (std::runtime_error e) {
        mexErrMsgTxt(e.what());
    }
    catch (...) {
        mexErrMsgTxt("Unknown error");
    }
}

/**
 Function that will write images to disk as a TIFF file. The numerical precision will be preserved.
 Any existing file will be overwritten.
 The input arguments must be of the form
 0 - the string "writeccdimages"
 1 - 3D matrix or string containing the full path to the file to read the data from.
 2 - the index of the first image to write (starting from 0)
 3 - the number of images to write (-1 for all images)
 4 - string containing the full path to the TIFF file to write
 */
void MatlabWriteCCDImages(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs) {
	// input validation
	if (nrhs != 5)
		mexErrMsgTxt("Must have exactly 5 input arguments for writeccdimages\n\
					 1. the string \"writeccdimages\"\n\
                     2. 2D or 3D matrix or string containing the full path to the file to read the data from.\n\
					 3. the index of the first image to write (starting from 0)\n\
					 4. the number of images to write (or -1 for all images)\n\
					 5. path to the TIFF file to write. Existing files will be overwritten");
	
	if (nlhs != 0)
		mexErrMsgTxt("no return argument for writeccdimages");
    
	mxArray* array;
	// the mxArray at index 0 will have been checked already by mexFunction
	
    // index 1 - must be a 2D or 3D matrix, or a string with a filepath.
    std::string inputFilePath;
	mxArray* dataArray = NULL;
	array = const_cast<mxArray*>(prhs[1]);
	if (mxGetClassID(array) == mxCHAR_CLASS) {
		inputFilePath = GetMatlabString(array);
		if (inputFilePath.empty())
			mexErrMsgTxt("Need a non-empty filepath to the data");
	} else {
		// assume that this is a valid data array
		// if it isn't then the image loader will throw an error
		dataArray = array;
	}
    
	// index 2 - must be a single number - the index of the first image to write
	array = const_cast<mxArray*>(prhs[2]);
	if ((mxGetN(array) != 1) || (mxGetM(array) != 1) || (mxGetClassID(array) != mxDOUBLE_CLASS))
		mexErrMsgTxt("3rd argument must be a double scalar: the index of the first image to write (starting from 0)");
	int firstImageToWrite = *(mxGetPr(array));
	if (firstImageToWrite < 0)
		firstImageToWrite = 0;
	
	// index 3 - must be a single number - the number of images to write
	array = const_cast<mxArray*>(prhs[3]);
	if ((mxGetN(array) != 1) || (mxGetM(array) != 1) || (mxGetClassID(array) != mxDOUBLE_CLASS))
		mexErrMsgTxt("4th argument must be a double scalar: the number of images to write (or -1 for all images)");
	int nImagesToWrite = *(mxGetPr(array));
	if (nImagesToWrite < 0)
		nImagesToWrite = -1;
    
	// index 4 - must be the output path
	std::string outputFilePath;
	array = const_cast<mxArray*>(prhs[4]);
	if (mxGetClassID(array) == mxCHAR_CLASS) {
		outputFilePath = GetMatlabString(array);
		if (outputFilePath.empty())
			mexErrMsgTxt("Need a non-empty filepath to the data");
	} else {
		mexErrMsgTxt("5th argument must be a string: the path to the output TIFF file)");
	}
    
	try {
		std::shared_ptr<ImageLoader> imageLoader;
		if (!inputFilePath.empty()) {
			imageLoader = GetImageLoader(inputFilePath);
		} else {
			imageLoader = std::shared_ptr<ImageLoader>(new ImageLoaderMatlab(dataArray));
		}
        LocalizerTIFFImageOutputWriter imageWriter(outputFilePath, true, false, imageLoader->getStorageType());
        Clip(firstImageToWrite, 0, imageLoader->getNImages() - 1);
        if (nImagesToWrite <= 0)
            nImagesToWrite = imageLoader->getNImages();
        Clip(nImagesToWrite, 1, imageLoader->getNImages() - firstImageToWrite);
        imageLoader->spoolTo(firstImageToWrite);
        for (int i = 0; i < nImagesToWrite; ++i) {
            imageWriter.write_image(imageLoader->readNextImage());
        }
	}
	catch (std::bad_alloc) {
        mexErrMsgTxt("Insufficient memory");
    }
    catch (int e) {
        mexErrMsgTxt("Int error");
    }
    catch (USER_ABORTED e) {
        mexErrMsgTxt("User abort");
    }
    catch (std::runtime_error e) {
        mexErrMsgTxt(e.what());
    }
    catch (...) {
        mexErrMsgTxt("Unknown error");
    }
}

std::string GetMatlabString(const mxArray* array) {
	char *str = NULL;
	str = mxArrayToString(array);
	if (str == NULL)
		mexErrMsgTxt("Error: cannot allocate string (non-string argument supplied?)");
	std::string stdStr(str);
	mxFree(str);
	return stdStr;
}

mxArray* ConvertImageToArray(ImagePtr image) {

	mxArray* outputMatrix = mxCreateDoubleMatrix(image->rows(), image->cols(), mxREAL);
	if (outputMatrix == NULL)
		throw std::bad_alloc();

	int nElements = image->rows() * image->cols();

	if (nElements > 0) {
		double *firstOutputElement = mxGetPr(outputMatrix);
		memcpy(firstOutputElement, image->data(), nElements * sizeof(double));
	}

	return outputMatrix;
}

mxArray* ConvertImagesToArray(const std::vector<ImagePtr>& images) {
    int nImages = images.size();
    int nRows, nCols;
    if (nImages == 0) {
        nRows = 0;
        nCols = 0;
    } else {
        nRows = images[0]->rows();
        nCols = images[0]->cols();
    }
    int nPixels = nRows * nCols;
    
    mwSize dims[3];
    dims[0] = nRows;
    dims[1] = nCols;
    dims[2] = nImages;
    mxArray* outputArray = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
	if (outputArray == NULL)
		throw std::bad_alloc();
    
    double* arrayPtr = mxGetPr(outputArray);
    for (int i = 0; i < nImages; ++i) {
        memcpy(arrayPtr + i * nPixels, images[i]->data(), nPixels * sizeof(double));
    }
    
    return outputArray;
}

mxArray* LoadImagesIntoArray(std::shared_ptr<ImageLoader> imageLoader, int firstImage, int nImagesToRead) {
    if ((firstImage < 0) || (firstImage > imageLoader->getNImages() - 1))
        throw std::runtime_error("invalid first image to read");
    if (nImagesToRead < 0)
        nImagesToRead = imageLoader->getNImages();
    nImagesToRead = Clip(nImagesToRead, 0, imageLoader->getNImages() - firstImage);
    MatlabImageOutputWriter outputWriter(nImagesToRead, imageLoader->getStorageType());
    imageLoader->spoolTo(firstImage);
    for (int i = 0; i < nImagesToRead; ++i) {
        outputWriter.write_image(imageLoader->readNextImage());
    }
    return outputWriter.getArray();
}
