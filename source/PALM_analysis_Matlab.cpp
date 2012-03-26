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

/**
 The Matlab interface exports just a single 'gateway' function, mexFunction. This is the only
 function that can be called from Matlab. By convention the first argument to this function should
 be a string that identifies the operation to be performed. Depending on the value of this string
 the rhs and lhs parameters passed here are applied to the more specific functions.
*/
void mexFunction(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs) {
	
	if (nrhs < 1)
		mexErrMsgTxt("Must have more than 1 input argument");
	
	const mxArray* inputArray = prhs[0];
	if ((mxGetClassID(inputArray) != mxCHAR_CLASS) || (mxGetM(inputArray) > 1))
		mexErrMsgTxt("First input argument must be a selector string");
	
	std::string selector = GetMatlabString(inputArray);
	if (boost::iequals(selector, "localize")) {
		MatlabLocalization(nlhs, plhs, nrhs, prhs);
	} else {
		mexErrMsgTxt("Unknown selector");
	}

}

/**
 Function that will handle localization. The input arguments must be of the form
 0 - the string "localize"
 1 - the expect standard deviation of the psf
 2 - segmentation algorithm, 'glrt' or 'smoothsigma'
 3 - the PFA if GLRT, otherwise smoothsigma factor
 4 - localization algorithm. '2DGauss' only for now
 5 - file path
*/
void MatlabLocalization(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs) {
	// input validation
	if (nrhs != 6)
		mexErrMsgTxt("Must have exactly 6 input arguments for localize");

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
	
	// index 4 - must be a string, currently only '2DGauss' is allowed
	array = prhs[4];
	if (mxGetClassID(array) != mxCHAR_CLASS)
		mexErrMsgTxt("5th argument must be a string (the localization algorithm)");
	std::string localizationStr = GetMatlabString(array);
	
	int localizationAlgorithm;
	if (boost::iequals(localizationStr, "2DGauss")) {
		localizationAlgorithm = LOCALIZATION_METHOD_2DGAUSS;
	} else {
		mexErrMsgTxt("Unknown localization algorithm");
	}
	
	// index 5 - must be the data file path
	array = prhs[5];
	if (mxGetClassID(array) != mxCHAR_CLASS)
		mexErrMsgTxt("6th argument must be a string (the path to the file containing the data)");
	
	std::string filePath = GetMatlabString(array);
	
	// now run the actual localization
	try {
		boost::shared_ptr<ImageLoader> imageLoader = GetImageLoader(filePath);
		
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
		switch (localizationAlgorithm) {
            case LOCALIZATION_METHOD_2DGAUSS:
                positionsFitter = boost::shared_ptr<FitPositions>(new FitPositions_SymmetricGaussian(psfWidth, 1.0));
                break;
            default:
                throw std::runtime_error("Unknown localization method");
                break;
		}
		
		boost::shared_ptr<ProgressReporter> progressReporter(new ProgressReporter_Silent);
		
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
    if (boost::algorithm::iequals(extension, "lsm"))
        return CAMERA_TYPE_TIFF;
	
    // if we're still here then the extension was not recognized
    throw std::runtime_error("Unable to deduce the file type");
}
