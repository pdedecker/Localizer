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

#include "LocalizerDLL.h"

EXPORT int DoLocalizationAnalysis(char *filePath,           // path to the file on disk
                                  double pfa,               // pfa for GLRT
                                  double psfWidth,          // width of the psf
                                  double** positionsArray,  // pointer to pointer to double that will be set to 2D array with results
                                  int* nPositions,          // number of positions in result
                                  int* nColumns) {          // number of columns in the array
    try {
        boost::shared_ptr<ImageLoader> imageLoader = GetImageLoader(std::string(filePath));
        boost::shared_ptr<ThresholdImage_Preprocessor> preprocessor(new ThresholdImage_Preprocessor_DoNothing());
        boost::shared_ptr<ThresholdImage> thresholder(new ThresholdImage_GLRT_FFT(pfa, psfWidth));
        boost::shared_ptr<ThresholdImage_Postprocessor> postprocessor(new ThresholdImage_Postprocessor_DoNothing());
        boost::shared_ptr<ParticleFinder> particleFinder(new ParticleFinder_adjacent4());
        boost::shared_ptr<FitPositions> fitPositions(new FitPositions_SymmetricGaussian(psfWidth, 1.0));
        std::vector<boost::shared_ptr<ParticleVerifier> > particleVerifiers;
        
        boost::shared_ptr<ProgressReporter> progressReporter(new ProgressReporter_Silent);
        
        boost::shared_ptr<PALMAnalysisController> analysisController(new PALMAnalysisController(thresholder, preprocessor,
                                                                                                postprocessor, particleFinder, particleVerifiers,
                                                                                                fitPositions,
                                                                                                progressReporter, -1, -1));
        boost::shared_ptr<LocalizedPositionsContainer> localizedPositions = analysisController->doPALMAnalysis(imageLoader);
    }
    catch (...) {
        return -1;
    }
    
    return 0;
}

boost::shared_ptr<ImageLoader> GetImageLoader(std::string filePath) {
    std::string fileExtension = filePath.substr(filePath.length() - 3, 3);
	
	// convert the extension to lowercase for easy comparison
	std::transform(fileExtension.begin(), fileExtension.end(), fileExtension.begin(), ::tolower);
	
	if (fileExtension == std::string("spe"))
		return boost::shared_ptr<ImageLoader>(new ImageLoaderSPE(filePath));
	
	if (fileExtension == std::string("sif"))
		return boost::shared_ptr<ImageLoader>(new ImageLoaderAndor(filePath));
	
	if (fileExtension == std::string("his"))
		return boost::shared_ptr<ImageLoader>(new ImageLoaderHamamatsu(filePath));
	
	if (fileExtension == std::string("tif"))
		return boost::shared_ptr<ImageLoader>(new ImageLoaderTIFF(filePath));
	
	if (fileExtension == std::string("pde"))
		return boost::shared_ptr<ImageLoader>(new ImageLoaderPDE(filePath));
	
	// if we get here then we don't recognize the file type
	throw (std::runtime_error("Unknown data file type with extension " + fileExtension));
}

EXPORT void LocalizerFreeArray(double *arrayPtr) {
    delete[] arrayPtr;
}
