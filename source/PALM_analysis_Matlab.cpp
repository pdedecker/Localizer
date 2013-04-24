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

#include "PALM_analysis_Matlab.h"

#include "PALM_analysis_defines.h"
#include "PALM_analysis.h"
#include "PALM_analysis_ParticleFinding.h"
#include "PALM_analysis_FileIO.h"
#include "PALM_analysis_SOFI.h"

/**
 The Matlab interface exports just a single 'gateway' function, mexFunction. This is the only
 function that can be called from Matlab. By convention the first argument to this function should
 be a string that identifies the operation to be performed. Depending on the value of this string
 the rhs and lhs parameters passed here are applied to the more specific functions.
*/
void mexFunction(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs) {
	
	if (nrhs < 1)
		mexErrMsgTxt("First argument must be a string describing the operation to perform (\"localize\", \"testsegmentation\", or \"sofi\")");
	
	const mxArray* inputArray = prhs[0];
	if ((mxGetClassID(inputArray) != mxCHAR_CLASS) || (mxGetM(inputArray) > 1))
		mexErrMsgTxt("First argument must be a string describing the operation to perform (\"localize\", \"testsegmentation\", or \"sofi\")");

	std::string selector = GetMatlabString(inputArray);
	if (boost::iequals(selector, "localize")) {
		MatlabLocalization(nlhs, plhs, nrhs, prhs);
	} else if (boost::iequals(selector, "testsegmentation")) {
		MatlabTestSegmentation(nlhs, plhs, nrhs, prhs);
	} else if (boost::iequals(selector, "sofi")) {
		MatlabSOFI(nlhs, plhs, nrhs, prhs);
	} else if (boost::iequals(selector, "readccdimages")) {
		MatlabReadCCDImages(nlhs, plhs, nrhs, prhs);
	} else {
		mexErrMsgTxt("Unknown selector (should be one of \"readccdimages\", \"localize\", \"testsegmentation\", or \"sofi\")");
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
		boost::shared_ptr<ImageLoader> imageLoader;
		if (!filePath.empty()) {
			imageLoader = GetImageLoader(filePath);
		} else {
			imageLoader = boost::shared_ptr<ImageLoader>(new ImageLoaderMatlab(dataArray));
		}
		
		boost::shared_ptr<ThresholdImage> thresholder;
		switch(segmentationAlgorithm) {
			case THRESHOLD_METHOD_GLRT:
				thresholder = boost::shared_ptr<ThresholdImage>(new ThresholdImage_GLRT_FFT(thresholdParam, psfWidth));
				break;
			case THRESHOLD_METHOD_SMOOTHSIGMA:
				thresholder = boost::shared_ptr<ThresholdImage>(new ThresholdImage_SmoothSigma(psfWidth, thresholdParam));
				break;
			default:
				throw std::runtime_error("Unknown segmentation method");
				break;
        }
		
		boost::shared_ptr<ThresholdImage_Preprocessor> preprocessor(new ThresholdImage_Preprocessor_DoNothing());
		boost::shared_ptr<ThresholdImage_Postprocessor> postprocessor(new ThresholdImage_Postprocessor_DoNothing());
		
		boost::shared_ptr<ParticleFinder> particleFinder(new ParticleFinder_adjacent8());
		
		std::vector<boost::shared_ptr<ParticleVerifier> > particleVerifiers;
		particleVerifiers.push_back(boost::shared_ptr<ParticleVerifier> (new ParticleVerifier_RemoveOverlappingParticles(psfWidth)));
		
		boost::shared_ptr<FitPositions> positionsFitter;
		double sigma = 1.0;
		switch (localizationAlgorithm) {
            case LOCALIZATION_METHOD_2DGAUSS:
                positionsFitter = boost::shared_ptr<FitPositions>(new FitPositions_SymmetricGaussian(psfWidth, sigma));
                break;
            case LOCALIZATION_METHOD_2DGAUSS_FIXEDWIDTH:
                positionsFitter = boost::shared_ptr<FitPositions>(new FitPositions_FixedWidthGaussian(psfWidth, sigma));
                break;
            case LOCALIZATION_METHOD_MULTIPLICATION:
                positionsFitter = boost::shared_ptr<FitPositions>(new FitPositionsMultiplication(psfWidth, sigma));
                break;
            case LOCALIZATION_METHOD_CENTROID:
                positionsFitter = boost::shared_ptr<FitPositions>(new FitPositionsCentroid(psfWidth));
                break;
            case LOCALIZATION_METHOD_2DGAUSS_ELLIPSOIDAL:
                positionsFitter = boost::shared_ptr<FitPositions>(new FitPositions_EllipsoidalGaussian_SymmetricPSF(psfWidth, sigma));
                break;
			case LOCALIZATION_METHOD_2DGAUSS_ELLIPSOIDAL_ASTIGMATISM:
				positionsFitter = boost::shared_ptr<FitPositions>(new FitPositions_EllipsoidalGaussian(psfWidth, sigma));
				break;
            case LOCALIZATION_METHOD_MLEwG:
                positionsFitter = boost::shared_ptr<FitPositions>(new FitPositions_MLEwG(psfWidth));
                break;
            default:
                throw std::runtime_error("Unknown localization method");
                break;
		}
		
		boost::shared_ptr<ProgressReporter> progressReporter(new ProgressReporter_MatlabWaitMex());
		
		boost::shared_ptr<PALMAnalysisController> analysisController(new PALMAnalysisController(thresholder, preprocessor,
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
		
		boost::shared_ptr<ThresholdImage> thresholder;
		switch(segmentationAlgorithm) {
			case THRESHOLD_METHOD_GLRT:	// the GLRT test proposed by Arnauld et al in Nat Methods 5:687 2008
				thresholder = boost::shared_ptr<ThresholdImage>(new ThresholdImage_GLRT_FFT(segmentationParameter, psfWidth));
				break;
			case THRESHOLD_METHOD_SMOOTHSIGMA:
				thresholder = boost::shared_ptr<ThresholdImage>(new ThresholdImage_SmoothSigma(psfWidth, segmentationParameter));
				break;
			default:
				throw std::runtime_error("Unknown segmentation method");
				break;
        }
		
		boost::shared_ptr<ThresholdImage_Preprocessor> preprocessor(new ThresholdImage_Preprocessor_DoNothing());
		boost::shared_ptr<ThresholdImage_Postprocessor> postprocessor(new ThresholdImage_Postprocessor_DoNothing());
		boost::shared_ptr<ParticleFinder> particlefinder(new ParticleFinder_adjacent8());
		
		// calculate the threshold
        boost::shared_ptr<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > thresholded_image = do_processing_and_thresholding(imageToSegment, preprocessor, thresholder, postprocessor);
		boost::shared_ptr<std::list<Particle> > located_particles = particlefinder->findPositions(imageToSegment, thresholded_image);
		
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
 1 - the desired order
 2 - single number: 0 for autocorrelation, 1 for crosscorrelation
 3 - file path or a matrix containing numeric data
 */
void MatlabSOFI(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs) {
	// input validation
	if (nrhs != 6)
		mexErrMsgTxt("Must have exactly 6 input arguments for sofi\n\
					 1. the string \"sofi\"\n\
					 2. the desired order (2 or 3)\n\
					 3. single number: 0 for autocorrelation, 1 for crosscorrelation\n\
					 4. number of frames to skip at the beginning\n\
					 5. number of frames to include in the calculation (<=0 means all frames up to the end)\n\
					 6. file path to a data file, or a 2D or 3D matrix containing numeric data");
	
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
	
	// index 3 - must be a single number (number of frames to skip)
	size_t nFramesToSkip;
	array = const_cast<mxArray*>(prhs[3]);
	if ((mxGetN(array) != 1) || (mxGetM(array) != 1) || (mxGetClassID(array) != mxDOUBLE_CLASS))
		mexErrMsgTxt("4th argument must be a double scalar (number of frames to skip)");
	nFramesToSkip = static_cast<size_t>(*mxGetPr(array) + 0.5);
	
	// index 4 - must be a single number (number of frames to include)
	size_t nFramesToInclude;
	array = const_cast<mxArray*>(prhs[4]);
	if ((mxGetN(array) != 1) || (mxGetM(array) != 1) || (mxGetClassID(array) != mxDOUBLE_CLASS))
		mexErrMsgTxt("5th argument must be a double scalar (number of frames to include)");
	if (*mxGetPr(array) <= 0) {
		nFramesToInclude = static_cast<size_t>(-1);
	} else {
		nFramesToSkip = static_cast<size_t>(*mxGetPr(array) + 0.5);
	}
	
	// index 5 - must be the data
	std::string filePath;
	mxArray* dataArray = NULL;
	array = const_cast<mxArray*>(prhs[5]);
	if (mxGetClassID(array) == mxCHAR_CLASS) {
		filePath = GetMatlabString(array);
		if (filePath.empty())
			mexErrMsgTxt("Need a non-empty filepath to the data");
	} else {
		// assume that this is a valid data array
		// if it isn't then the image loader will throw an error
		dataArray = array;
	}
	
	size_t nRows = mxGetN(array);
	size_t nCols = mxGetM(array);
	
	try {
		boost::shared_ptr<ImageLoader> imageLoader;
		if (!filePath.empty()) {
			imageLoader = GetImageLoader(filePath);
		} else {
			imageLoader = boost::shared_ptr<ImageLoader>(new ImageLoaderMatlab(dataArray));
		}
		
		boost::shared_ptr<ProgressReporter> progressReporter(new ProgressReporter_MatlabWaitMex());
		
		// no frame verifiers for now
		std::vector<boost::shared_ptr<SOFIFrameVerifier> > frameVerifiers;
		int lagTime = 0;
		std::vector<ImagePtr> sofiOutputImages, averageOutputImages;
		DoSOFIAnalysis(imageLoader, frameVerifiers, progressReporter,nFramesToSkip, nFramesToInclude, lagTime, correlationOrder, doCrossCorrelation, 50, sofiOutputImages, averageOutputImages);

		plhs[0] = ConvertImageToArray(sofiOutputImages.at(0));
		plhs[1] = ConvertImageToArray(averageOutputImages.at(0));
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
					 3. the number of images to read\n\
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
		firstImageToRead = -1;
	
	// index 2 - must be a single number - the number of images to read
	bool doCrossCorrelation;
	array = const_cast<mxArray*>(prhs[2]);
	if ((mxGetN(array) != 1) || (mxGetM(array) != 1) || (mxGetClassID(array) != mxDOUBLE_CLASS))
		mexErrMsgTxt("3nd argument must be a double scalar: the number of images to read (or -1)");
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
		boost::shared_ptr<ImageLoader> imageLoader;
		if (!filePath.empty()) {
			imageLoader = GetImageLoader(filePath);
		} else {
			imageLoader = boost::shared_ptr<ImageLoader>(new ImageLoaderMatlab(dataArray));
		}

		plhs[0] = ConvertImagesToArray(imageLoader, firstImageToRead, nImagesToRead);
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

boost::shared_ptr<ImageLoader> GetImageLoader(std::string& data_file_path) {
    boost::shared_ptr<ImageLoader> image_loader;
	
    int estimatedCameraType = GetFileStorageType(data_file_path);
	
    switch (estimatedCameraType) {
		case CAMERA_TYPE_WINSPEC:	// spe files
			image_loader = boost::shared_ptr<ImageLoader>(new ImageLoaderSPE(data_file_path));
			break;
		case CAMERA_TYPE_ANDOR:
			image_loader = boost::shared_ptr<ImageLoader>(new ImageLoaderAndor(data_file_path));
			break;
		case CAMERA_TYPE_HAMAMATSU:
			image_loader = boost::shared_ptr<ImageLoader>(new ImageLoaderHamamatsu(data_file_path));
			break;
		case CAMERA_TYPE_TIFF:	// 3 is reserved for TIFF files
			image_loader = boost::shared_ptr<ImageLoader>(new ImageLoaderTIFF(data_file_path));
			break;
		case CAMERA_TYPE_PDE:
			image_loader = boost::shared_ptr<ImageLoader>(new ImageLoaderPDE(data_file_path));
			break;
		case CAMERA_TYPE_ZEISS:	// Zeiss lsm files
			image_loader = boost::shared_ptr<ImageLoader>(new ImageLoaderTIFF(data_file_path));
			break;
		default:
			throw std::runtime_error("Unsupported CCD file type");
			break;
    }
	
    return image_loader;
	
}

int GetFileStorageType(std::string &filePath) {
    size_t startOfExtension = filePath.rfind('.');
    if (startOfExtension == size_t(-1)) {
        // the filepath does not appear to contain an extension
        throw std::runtime_error("Unable to deduce the file type");
    }
	
    std::string extension = filePath.substr(startOfExtension + 1);
    if ((extension.length() < 3) || (extension.length() > 4)) {
        throw std::runtime_error("Unable to deduce the file type");
    }
	
    if (boost::algorithm::iequals(extension, "spe"))
        return CAMERA_TYPE_WINSPEC;
    if (boost::algorithm::iequals(extension, "sif"))
        return CAMERA_TYPE_ANDOR;
    if (boost::algorithm::iequals(extension, "his"))
        return CAMERA_TYPE_HAMAMATSU;
    if (boost::algorithm::iequals(extension, "pde"))
        return CAMERA_TYPE_PDE;
    if (boost::algorithm::iequals(extension, "tif"))
        return CAMERA_TYPE_TIFF;
    if (boost::algorithm::iequals(extension, "tiff"))
        return CAMERA_TYPE_TIFF;
	if (boost::algorithm::iequals(extension, "btf"))
        return CAMERA_TYPE_TIFF;
	if (boost::algorithm::iequals(extension, "tf8"))
        return CAMERA_TYPE_TIFF;
    if (boost::algorithm::iequals(extension, "lsm"))
        return CAMERA_TYPE_TIFF;
	
    // if we're still here then the extension was not recognized
    throw std::runtime_error("Unable to deduce the file type");
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

mxArray* ConvertImagesToArray(boost::shared_ptr<ImageLoader> imageLoader, int firstImage, int nImagesToRead) {
	size_t nImages = imageLoader->getNImages();
	size_t xSize = imageLoader->getXSize();
	size_t ySize = imageLoader->getYSize();

	// allocate the output matrix
	// try to match the required storage type
	int storageType = imageLoader->getStorageType();
	mxClassID classID;
	size_t bytesPixel;
	switch (storageType) {
		case STORAGE_TYPE_INT8:
			classID = mxINT8_CLASS;
			bytesPixel = 1;
			break;
		case STORAGE_TYPE_UINT8:
			classID = mxUINT8_CLASS;
			bytesPixel = 1;
			break;
		case STORAGE_TYPE_INT16:
			classID = mxINT16_CLASS;
			bytesPixel = 2;
			break;
		case STORAGE_TYPE_UINT16:
			classID = mxUINT16_CLASS;
			bytesPixel = 2;
			break;
		case STORAGE_TYPE_INT32:
			classID = mxINT32_CLASS;
			bytesPixel = 4;
			break;
		case STORAGE_TYPE_UINT32:
			classID = mxUINT32_CLASS;
			bytesPixel = 4;
			break;
		case STORAGE_TYPE_INT64:
			classID = mxINT8_CLASS;
			bytesPixel = 8;
			break;
		case STORAGE_TYPE_UINT64:
			classID = mxUINT8_CLASS;
			bytesPixel = 8;
			break;
			case STORAGE_TYPE_FP32:
			classID = mxSINGLE_CLASS;
			bytesPixel = 4;
			break;
		case STORAGE_TYPE_FP64:
			classID = mxDOUBLE_CLASS;
			bytesPixel = 8;
			break;
		default:
			throw std::runtime_error("unsupported storage type in ConvertImagesToArray()");
			break;
	}

	if (firstImage < 0)
		firstImage = 0;
	if ((nImagesToRead < 0) || (nImagesToRead > nImages - firstImage))
		nImagesToRead = nImages - firstImage;

	mwSize ndim = 3;
	mwSize dims[3] = {xSize, ySize, nImagesToRead};
	mxArray* outputMatrix = mxCreateNumericArray(ndim, dims, classID, mxREAL);
	if (outputMatrix == NULL)
		throw std::bad_alloc();

	size_t bytesPerImage = xSize * ySize * bytesPixel;
	char* arrayPtr = reinterpret_cast<char*>(mxGetPr(outputMatrix));

	imageLoader->spoolTo(firstImage);
	for (size_t i = 0; i < nImagesToRead; ++i) {
		ImagePtr thisImage = imageLoader->readNextImage();
		char* bufferPtr = arrayPtr + i * bytesPerImage;
		switch (storageType) {
			case STORAGE_TYPE_INT8:
				CopyImageToBuffer<int8_t>(thisImage, bufferPtr);
				break;
			case STORAGE_TYPE_UINT8:
				CopyImageToBuffer<uint8_t>(thisImage, bufferPtr);
				break;
			case STORAGE_TYPE_INT16:
				CopyImageToBuffer<int16_t>(thisImage, bufferPtr);
				break;
			case STORAGE_TYPE_UINT16:
				CopyImageToBuffer<uint16_t>(thisImage, bufferPtr);
				break;
			case STORAGE_TYPE_INT32:
				CopyImageToBuffer<int32_t>(thisImage, bufferPtr);
				break;
			case STORAGE_TYPE_UINT32:
				CopyImageToBuffer<uint32_t>(thisImage, bufferPtr);
				break;
			case STORAGE_TYPE_INT64:
				CopyImageToBuffer<int64_t>(thisImage, bufferPtr);
				break;
			case STORAGE_TYPE_UINT64:
				CopyImageToBuffer<uint64_t>(thisImage, bufferPtr);
				break;
			case STORAGE_TYPE_FP32:
				CopyImageToBuffer<float>(thisImage, bufferPtr);
				break;
			case STORAGE_TYPE_FP64:
				CopyImageToBuffer<double>(thisImage, bufferPtr);
			break;
		default:
			throw std::runtime_error("unsupported storage type in ConvertImagesToArray()");
			break;
		}
	}

	return outputMatrix;
}
