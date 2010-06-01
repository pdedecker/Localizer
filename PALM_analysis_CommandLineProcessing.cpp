/*
 *  PALM_analysis_CommandLineProcessing.cpp
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 01/06/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "PALM_analysis_CommandLineProcessing.h"

namespace po = boost::program_options;



// the main function
int main(int argc, char *argv[]) {
	po::options_description desc("Allowed options");
	desc.add_options()
	("help", "produce help message")
	("process", po::value<std::string>(), "Process the CCD files and save them as a converted stack. If this options is present then it takes precedence and no localization will be done! Options are \"subtractaverage\", \"differenceimage\", \"converttophotons\".")
	("averaging", po::value<size_t>()->default_value(0), "subtractaverage: number of frames to average over.")
	("cameramultiplier", po::value<double>(), "converttophotons: camera multiplication factor.")
	("cameraoffset", po::value<double>(), "converttophotons: camera offset.")
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
	
	std::string processorName;
	std::vector<std::string> inputFiles = vm["input-file"].as<std::vector <std::string> >();
	size_t nInputFiles = inputFiles.size();
	
	boost::shared_ptr<CCDImagesProcessor> ccdImagesProcessor;
	
	// we want to process the images, not localize
	// get the processor
	processorName = vm["process"].as<std::string>();
	size_t nFramesAveraging = vm["averaging"].as<size_t>();
	
	double cameraMultiplicationFactor, cameraOffset;
	// check if the required options are available
	if (processorName == std::string("converttophotons")) {
		if (vm.count("cameramultiplier") == 0)
			throw std::runtime_error("A camera multiplier factor needs to be specified when converting to photons");
		else
			cameraMultiplicationFactor = vm["cameramultiplier"].as<double>();
		
		if (vm.count("cameraoffset") == 0)
			throw std::runtime_error("A camera offset needs to be specified when converting to photons");
		else
			cameraOffset = vm["cameraoffset"].as<double>();
	}
	
	// get a progress reporter
	boost::shared_ptr<PALMAnalysisProgressReporter> progressReporter(new PALMAnalysisProgressReporter_stdout);
	
	ccdImagesProcessor = GetCCDImagesProcessor(processorName, progressReporter, nFramesAveraging, cameraMultiplicationFactor, cameraOffset);
	
	boost::shared_ptr<ImageLoader> imageLoader;
	boost::shared_ptr<ImageOutputWriter> imageOutputWriter;
	std::string outputFilePath;
		
	for (size_t i = 0; i < nInputFiles; ++i) {
		try {
			// if there is more than one input file then tell the user which inputfile is being processed
			if (nInputFiles > 1) {
				std::cout << "Now processing " << inputFiles.at(i) << " (" << i + 1 << " of " << nInputFiles << ")\n";
				std::cout.flush();
			}
			
			// get an imageloader and output writer
			imageLoader = GetImageLoader(inputFiles.at(i));
			outputFilePath = GetOutputProcessedImagesFilePath(inputFiles.at(i));
			imageOutputWriter = GetImageOutputWriter(processorName, outputFilePath, COMPRESSION_NONE);
			
			// do the processing
			ccdImagesProcessor->convert_images(imageLoader, imageOutputWriter);
			
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

boost::shared_ptr<CCDImagesProcessor> GetCCDImagesProcessor(std::string name, boost::shared_ptr<PALMAnalysisProgressReporter> progressReporter, size_t nFramesAveraging, double cameraMultiplier, double cameraOffset) {
	if (name == std::string("subtractaverage"))
		return boost::shared_ptr<CCDImagesProcessor> (new CCDImagesProcessorAverageSubtraction(progressReporter, nFramesAveraging));
	
	if (name == std::string("differenceimage"))
		return boost::shared_ptr<CCDImagesProcessor> (new CCDImagesProcessorDifferenceImage(progressReporter));
	
	if (name == std::string("converttophotons"))
		return boost::shared_ptr<CCDImagesProcessor> (new CCDImagesProcessorConvertToPhotons(progressReporter, cameraMultiplier, cameraOffset));
	
	// if we get here then we didn't recognize the processing algorithm
	throw std::runtime_error("Unknown processing algorithm");
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
	
	if (fileExtension == std::string("sif"))
		return boost::shared_ptr<ImageLoader>(new ImageLoaderTIFF(filePath));
	
	// if we get here then we don't recognize the file type
	throw (std::runtime_error("Unknown data file type with extension " + fileExtension));
	
}

boost::shared_ptr<ImageOutputWriter> GetImageOutputWriter(std::string processMethodName, std::string outputFilePath, size_t compression) {
	if (processMethodName == std::string("subtractaverage"))
		return boost::shared_ptr<ImageOutputWriter> (new TIFFImageOutputWriter(outputFilePath, 1, compression, STORAGE_TYPE_FP32));
	
	if (processMethodName == std::string("differenceimage"))
		return boost::shared_ptr<ImageOutputWriter> (new TIFFImageOutputWriter(outputFilePath, 1, compression, STORAGE_TYPE_FP32));
	
	if (processMethodName == std::string("converttophotons"))
		return boost::shared_ptr<ImageOutputWriter> (new TIFFImageOutputWriter(outputFilePath, 1, compression, STORAGE_TYPE_UINT32));
	
	// if we get here then we didn't recognize the processing algorithm
	throw std::runtime_error("Unknown processing algorithm");
}

std::string GetOutputPositionsFilePath(std::string dataFilePath) {
	// remove the last 4 characters and replace with "_positions.txt"
	
	std::string outputPath = dataFilePath;
	outputPath.erase(outputPath.length() - 4, 4);
	outputPath += "_positions.txt";
	return outputPath;
}

std::string GetOutputProcessedImagesFilePath(std::string dataFilePath) {
	// remove the last 4 characters and replace with "_processed.tif"
	
	std::string outputPath = dataFilePath;
	outputPath.erase(outputPath.length() - 4, 4);
	outputPath += "_processed.tif";
	return outputPath;
}