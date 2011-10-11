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

#include "PALM_analysis_CommandLine.h"
namespace po = boost::program_options;



// the main function
int main(int argc, char *argv[]) {
	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "produce help message")
		("preprocessing", po::value<std::string>()->default_value("none"), "segmentation preprocessing. Typically not required.")
		("postprocessing", po::value<std::string>()->default_value("none"), "segmentation postprocessing. Typically not required.")
		("segmentation", po::value<std::string>()->default_value("glrt"), "segmentation algorithm to use. The only recommended option is \"glrt\".")
		("particlefinding", po::value<std::string>()->default_value("4way"), "particle finding algorithm to use. Options are \"4way\", and \"8way\"")
		("particleverifier", po::value<std::vector<std::string> >()->composing(), "particle verification to use. Multiple options are allowed. Options are \"none\", \"removeoverlapping\", \"symmetric2dgauss\" and \"ellipsoidal2dgauss\".")
		("localization", po::value<std::string>()->default_value("symmetric2dgauss"), "localization algorithm to use. Options are \"symmetric2dgauss\", \"symmetric2dgaussfixedwidth\", \"ellipsoidal2dgauss\", \"multiplication\", \"centroid\", and \"mlewg\".")
		("pfa", po::value<double>(), "Threshold parameter for GLRT localization.")
		("threshold", po::value<double>(), "Threshold parameter for direct thresholding.")
		("psf-width", po::value<double>()->default_value(2.0), "Estimated standard deviation of the PSF (in pixels).")
		("input-file", po::value< std::vector<std::string> >(), "Input file containing CCD images")
	;
	
	// tell program options that arguments without a flag are input CCD files
	po::positional_options_description p;
	p.add("input-file", -1);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).
			  options(desc).positional(p).run(), vm);
	po::notify(vm);
	
	// if requested print a help message
	if (vm.count("help")) {
		std::cout << desc << "\n";
		return 1;
	}
	
	// did the user specify any input files?
	if (vm.count("input-file") == 0) {
		std::cout << "No input files specified. Aborting... (try adding \"--help\" to get a help message)\n";
		return 1;
	}
	
	std::string preprocessingName = vm["preprocessing"].as<std::string>();
	std::string postprocessingName = vm["postprocessing"].as<std::string>();
	std::string segmentationName = vm["segmentation"].as<std::string>();
	std::string particleFinderName = vm["particlefinding"].as<std::string>();
	std::string localizationName = vm["localization"].as<std::string>();
	std::vector<std::string> inputFiles = vm["input-file"].as<std::vector <std::string> >();
	
	std::vector<std::string> particleVerifierNames;
	if (vm.count("particleverifier") != 0) {
		particleVerifierNames = vm["particleverifier"].as<std::vector<std::string> >();
	}
	
	
	// glrt and direct thresholding require extra parameters
	double pfa, directThreshold;
	if (segmentationName == std::string("glrt")) {
		if (vm.count("pfa") == 0) {
			throw std::runtime_error("The 'glrt' segmentation algorithm was chosen, but a pfa value was not specified (--pfa flag)");
		} else {
			pfa = vm["pfa"].as<double>();
		}
	}
	
	if (segmentationName == std::string("direct")) {
		if (vm.count("pfa") == 0) {
			throw std::runtime_error("The 'direct' segmentation algorithm was chosen, but a threshold value was not specified (--threshold flag)");
		} else {
			directThreshold = vm["threshold"].as<double>();
		}
	}
	
	double psfWidth = vm["psf-width"].as<double>();
	
	// get the pre- and postprocessors
	boost::shared_ptr<ThresholdImage_Preprocessor> preprocessor = GetPreProcessorType(vm["preprocessing"].as<std::string>());
	boost::shared_ptr<ThresholdImage_Postprocessor> postprocessor = GetPostProcessorType(vm["postprocessing"].as<std::string>());
	
	// get the thresholder
	boost::shared_ptr<ThresholdImage> thresholder = GetSegmentationType(segmentationName, pfa, directThreshold, psfWidth);
	
	// get the particle finder
	boost::shared_ptr<ParticleFinder> particleFinder = GetParticleFinderType(particleFinderName);
	
	// get the positions verifier
	std::vector<boost::shared_ptr<ParticleVerifier> > particleVerifiers;
	for (std::vector<std::string>::iterator it = particleVerifierNames.begin(); it != particleVerifierNames.end(); ++it) {
		particleVerifiers.push_back(GetParticleVerifierType(*it, psfWidth, 1.0));
	}
	
	// get the positions fitter
	boost::shared_ptr<FitPositions> positionsFitter = GetPositionsFitter(localizationName, psfWidth);
	
	// get a progress reporter
	boost::shared_ptr<ProgressReporter> progressReporter(new ProgressReporter_stdout());
	
	// get an analysis controller
	boost::shared_ptr<PALMAnalysisController> analysisController (new PALMAnalysisController(thresholder, preprocessor, 
																							postprocessor, particleFinder, particleVerifiers,
																							positionsFitter,
																							progressReporter));
	// run the PALM analysis for every input file
	size_t nInputFiles = inputFiles.size();
	if (nInputFiles == 0)
		std::cout << "No input files specified!\n";
	
	boost::shared_ptr<ImageLoader> imageLoader;
	std::string outputFilePath;
	std::ostringstream header;
	boost::shared_ptr<LocalizedPositionsContainer> fittedPositions;
	
	for (size_t i = 0; i < nInputFiles; ++i) {
		try {
			// get an imageloader and output file path
			imageLoader = GetImageLoader(inputFiles.at(i));
			outputFilePath = GetOutputPositionsFilePath(inputFiles.at(i));
			
			// if there is more than one input file then tell the user which inputfile is being processed
			if (nInputFiles > 1) {
				std::cout << "Now processing " << inputFiles.at(i) << " (" << i + 1 << " of " << nInputFiles << ")\n";
				std::cout.flush();
			}
			
			// do the analysis
			fittedPositions = analysisController->DoPALMAnalysis(imageLoader);
			
			// write a header
			header.str("");
			header << "X SIZE:" << imageLoader->getXSize() << "\n";
			header << "Y SIZE:" << imageLoader->getYSize() << "\n";
			header << "PFA:" << pfa << "\n";
			header << "PSF WIDTH:" << psfWidth << "\n";
			header << "THRESHOLD METHOD:" << segmentationName << "\n";
			header << "PARTICLE FINDING:" << particleFinderName << "\n";
			header << "PARTICLE VERIFIER:";
			for (std::vector<std::string>::iterator it = particleVerifierNames.begin(); it != particleVerifierNames.end(); ++it) {
				header << *it << ',';
			}
			if (particleVerifierNames.size() == 0)
				header << "none";
			header << "\n";
			header << "LOCALIZATION METHOD:" << localizationName << "\n";
			
			// write the output
			fittedPositions->writePositionsToFile(outputFilePath, header.str());
		}
		catch (std::runtime_error e) {
			std::cerr << "the file at " << inputFiles.at(i) << " failed due to " << e.what() << std::endl;
		}
		catch (...) {
			std::cerr << "the file at " << inputFiles.at(i) << " failed for an unknown reason" << std::endl;
		}
			
	}
	
	return 0;
}

boost::shared_ptr<ThresholdImage_Preprocessor> GetPreProcessorType(std::string name) {
	
	if (name == std::string("none"))
		return boost::shared_ptr<ThresholdImage_Preprocessor>(new ThresholdImage_Preprocessor_DoNothing());
	
	if (name == std::string("3x3median"))
		return boost::shared_ptr<ThresholdImage_Preprocessor>(new ThresholdImage_Preprocessor_MedianFilter(3,3));
	
	if (name == std::string("5x5median"))
		return boost::shared_ptr<ThresholdImage_Preprocessor>(new ThresholdImage_Preprocessor_MedianFilter(5,5));
	
	if (name == std::string("1x1gaussian"))
		return boost::shared_ptr<ThresholdImage_Preprocessor>(new ThresholdImage_Preprocessor_GaussianSmoothing(1));
	
	if (name == std::string("2x2gaussian"))
		return boost::shared_ptr<ThresholdImage_Preprocessor>(new ThresholdImage_Preprocessor_GaussianSmoothing(2));
	
	if (name == std::string("3x3mean"))
		return boost::shared_ptr<ThresholdImage_Preprocessor>(new ThresholdImage_Preprocessor_MeanFilter(3,3));
	
	if (name == std::string("5x5mean"))
		return boost::shared_ptr<ThresholdImage_Preprocessor>(new ThresholdImage_Preprocessor_MeanFilter(5,5));
	
	// if we get here then we didn't recognize the preprocessor
	throw std::runtime_error("Unknown preprocessor option");
}

boost::shared_ptr<ThresholdImage_Postprocessor> GetPostProcessorType(std::string name) {
	if (name == std::string("none"))
		return boost::shared_ptr<ThresholdImage_Postprocessor>(new ThresholdImage_Postprocessor_DoNothing());
	
	if (name == std::string("removeisolated"))
		return boost::shared_ptr<ThresholdImage_Postprocessor>(new ThresholdImage_Postprocessor_RemoveIsolatedPixels());
	
	// if we get here then we didn't recognize the postprocessor
	throw std::runtime_error("Unknown postprocessor algorithm");
}

boost::shared_ptr<ThresholdImage> GetSegmentationType(std::string name, double pfa, double threshold, double psfWidth) {
	if (name == std::string("glrt"))
		return boost::shared_ptr<ThresholdImage>(new ThresholdImage_GLRT_FFT(pfa, psfWidth));
	
	if (name == std::string("isodata"))
		return boost::shared_ptr<ThresholdImage>(new ThresholdImage_Isodata());
	
	if (name == std::string("triangle"))
		return boost::shared_ptr<ThresholdImage>(new ThresholdImage_Triangle());
	
	if (name == std::string("direct"))
		return boost::shared_ptr<ThresholdImage>(new ThresholdImage_Direct(threshold));
	
	// if we get here then we didn't recognize the preprocessor
	throw std::runtime_error("Unknown segmentation algorithm");
}

boost::shared_ptr<ParticleFinder> GetParticleFinderType(std::string name) {
	if (name == std::string("4way"))
		return boost::shared_ptr<ParticleFinder>(new ParticleFinder_adjacent4());
	
	if (name == std::string("8way"))
		return boost::shared_ptr<ParticleFinder>(new ParticleFinder_adjacent8());
	
	// if we get here then we didn't recognize the particlefinder
	throw std::runtime_error("Unknown particle finding algorithm");
}

boost::shared_ptr<ParticleVerifier> GetParticleVerifierType(std::string name, double psfWidth, double sigma) {
	if (name == std::string("none"))
		return boost::shared_ptr<ParticleVerifier>();
	
	if (name == std::string("removeoverlapping"))
		return boost::shared_ptr<ParticleVerifier>(new ParticleVerifier_RemoveOverlappingParticles(psfWidth));
	
	if (name == std::string("symmetric2dgauss"))
		return boost::shared_ptr<ParticleVerifier>(new ParticleVerifier_SymmetricGaussian(psfWidth, sigma));
	
	if (name == std::string("ellipsoidal2dgauss"))
		return boost::shared_ptr<ParticleVerifier>(new ParticleVerifier_EllipsoidalGaussian(psfWidth, sigma));
	
	// if we get here then we didn't recognize the particleverifier
	throw std::runtime_error("Unknown particle verification algorithm");
}
	

boost::shared_ptr<FitPositions> GetPositionsFitter(std::string name, double psfWidth) {
	double sigma = 1;
	
	if (name == std::string("symmetric2dgauss"))
		return boost::shared_ptr<FitPositions>(new FitPositions_SymmetricGaussian(psfWidth, sigma));
	
	if (name == std::string("symmetric2dgaussfixedwidth"))
		return boost::shared_ptr<FitPositions>(new FitPositions_FixedWidthGaussian(psfWidth, sigma));
	
	if (name == std::string("ellipsoidal2dgauss"))
		return boost::shared_ptr<FitPositions>(new FitPositions_EllipsoidalGaussian(psfWidth, sigma));
	
	if (name == std::string("multiplication"))
		return boost::shared_ptr<FitPositions>(new FitPositionsMultiplication(psfWidth, sigma));
	
	if (name == std::string("centroid"))
		return boost::shared_ptr<FitPositions>(new FitPositionsCentroid(psfWidth));
	
	if (name == std::string("mlewg"))
		return boost::shared_ptr<FitPositions>(new FitPositions_MLEwG(psfWidth));
	
	// if we get here then we didn't recognize the localizing algorithm
	throw std::runtime_error("Unknown localization algorithm");
}

boost::shared_ptr<ImageLoader> GetImageLoader(std::string filePath) {
	// look at the file extension
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

std::string GetOutputPositionsFilePath(std::string dataFilePath) {
	// remove the last 4 characters and replace with "_positions.txt"
	
	std::string outputPath = dataFilePath;
	outputPath.erase(outputPath.length() - 4, 4);
	outputPath += "_positions.txt";
	return outputPath;
}

