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

#include "PALM_analysis_CommandLineProcessing.h"

namespace po = boost::program_options;



// the main function
int main(int argc, char *argv[]) {
	po::options_description desc("Allowed options");
	desc.add_options()
	("help", "produce help message")
	("process", po::value<std::string>()->default_value("convertfileformat"), "Process the CCD files and save them as a converted stack. Options are \"subtractaverage\", \"differenceimage\", \"converttophotons\", \"convertfileformat\".")
	("outputformat", po::value<std::string>()->default_value("pde"), "Select the output format: \"tiff\", \"compressedtiff\" or \"pde\".")
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
	
	std::string processorName = vm["process"].as<std::string>();
	std::string outputFormat = vm["outputformat"].as<std::string>();
	std::vector<std::string> inputFiles = vm["input-file"].as<std::vector <std::string> >();
	size_t nInputFiles = inputFiles.size();
	int originalStorageFormat;
	
	std::shared_ptr<CCDImagesProcessor> ccdImagesProcessor;
	
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
	std::shared_ptr<ProgressReporter> progressReporter(new ProgressReporter_stdout);
	
	ccdImagesProcessor = GetCCDImagesProcessor(processorName, progressReporter, nFramesAveraging, cameraMultiplicationFactor, cameraOffset);
	
	std::shared_ptr<ImageLoader> imageLoader;
	std::shared_ptr<ImageOutputWriter> imageOutputWriter;
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
			originalStorageFormat = imageLoader->getStorageType();
			outputFilePath = GetOutputProcessedImagesFilePath(inputFiles.at(i), processorName, outputFormat);
			imageOutputWriter = GetImageOutputWriter(processorName, originalStorageFormat, outputFormat, outputFilePath);
			
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

std::shared_ptr<CCDImagesProcessor> GetCCDImagesProcessor(std::string name, std::shared_ptr<ProgressReporter> progressReporter, size_t nFramesAveraging, double cameraMultiplier, double cameraOffset) {
	if (name == std::string("subtractaverage"))
		return std::shared_ptr<CCDImagesProcessor> (new CCDImagesProcessorAverageSubtraction(progressReporter, nFramesAveraging));
	
	if (name == std::string("differenceimage"))
		return std::shared_ptr<CCDImagesProcessor> (new CCDImagesProcessorDifferenceImage(progressReporter));
	
	if (name == std::string("converttophotons"))
		return std::shared_ptr<CCDImagesProcessor> (new CCDImagesProcessorConvertToPhotons(progressReporter, cameraMultiplier, cameraOffset));
	
	if (name == std::string("convertfileformat"))
		return std::shared_ptr<CCDImagesProcessor> (new CCDImagesProcessorConvertToSimpleFileFormat(progressReporter));
	
	// if we get here then we didn't recognize the processing algorithm
	throw std::runtime_error("Unknown processing algorithm");
}

std::shared_ptr<ImageLoader> GetImageLoader(std::string filePath) {
	// look at the file extension
	std::string fileExtension = filePath.substr(filePath.length() - 3, 3);
	
	// convert the extension to lowercase for easy comparison
	std::transform(fileExtension.begin(), fileExtension.end(), fileExtension.begin(), ::tolower);
	
	if (fileExtension == std::string("spe"))
		return std::shared_ptr<ImageLoader>(new ImageLoaderSPE(filePath));
	
	if (fileExtension == std::string("sif"))
		return std::shared_ptr<ImageLoader>(new ImageLoaderAndor(filePath));
	
	if (fileExtension == std::string("his"))
		return std::shared_ptr<ImageLoader>(new ImageLoaderHamamatsu(filePath));
	
	if (fileExtension == std::string("tif"))
		return std::shared_ptr<ImageLoader>(new ImageLoaderTIFF(filePath));
	
	if (fileExtension == std::string("pde"))
		return std::shared_ptr<ImageLoader>(new ImageLoaderPDE(filePath));
	
	// if we get here then we don't recognize the file type
	throw (std::runtime_error("Unknown data file type with extension " + fileExtension));
	
}

std::shared_ptr<ImageOutputWriter> GetImageOutputWriter(std::string processMethodName, int originalStorageFormat, std::string requestedFormat, std::string outputFilePath) {
	size_t compression;
	
	if ((requestedFormat == std::string("tiff")) || (requestedFormat == std::string("compressedtiff"))) {
		if (requestedFormat == std::string("compressedtiff"))
			compression = COMPRESSION_DEFLATE;
		else
			compression = COMPRESSION_NONE;

		if (processMethodName == std::string("subtractaverage"))
			return std::shared_ptr<ImageOutputWriter> (new TIFFImageOutputWriter(outputFilePath, 0, compression, STORAGE_TYPE_FP32));
		
		if (processMethodName == std::string("differenceimage"))
			return std::shared_ptr<ImageOutputWriter> (new TIFFImageOutputWriter(outputFilePath, 0, compression, STORAGE_TYPE_FP32));
		
		if (processMethodName == std::string("converttophotons"))
			return std::shared_ptr<ImageOutputWriter> (new TIFFImageOutputWriter(outputFilePath, 0, compression, STORAGE_TYPE_UINT32));
		
		if (processMethodName == std::string("convertfileformat"))
			return std::shared_ptr<ImageOutputWriter> (new TIFFImageOutputWriter(outputFilePath, 0, compression, originalStorageFormat));
		
		// if we get here then we didn't recognize the processing algorithm
		throw std::runtime_error("Unknown processing algorithm");
	}
	
	if (requestedFormat == std::string("pde")) {
		if (processMethodName == std::string("subtractaverage"))
			return std::shared_ptr<ImageOutputWriter> (new PDEImageOutputWriter(outputFilePath, 0, STORAGE_TYPE_FP32));
		
		if (processMethodName == std::string("differenceimage"))
			return std::shared_ptr<ImageOutputWriter> (new PDEImageOutputWriter(outputFilePath, 0, STORAGE_TYPE_FP32));
		
		if (processMethodName == std::string("converttophotons"))
			return std::shared_ptr<ImageOutputWriter> (new PDEImageOutputWriter(outputFilePath, 0, STORAGE_TYPE_UINT32));
		
		if (processMethodName == std::string("convertfileformat"))
			return std::shared_ptr<ImageOutputWriter> (new PDEImageOutputWriter(outputFilePath, 0, originalStorageFormat));
		
		// if we get here then we didn't recognize the processing algorithm
		throw std::runtime_error("Unknown processing algorithm");
	}
	
	// if we get here then we didn't recognize the output format
	throw std::runtime_error("Unknown output format");
}

std::string GetOutputPositionsFilePath(std::string dataFilePath) {
	// remove the last 4 characters and replace with "_positions.txt"
	
	std::string outputPath = dataFilePath;
	outputPath.erase(outputPath.length() - 4, 4);
	outputPath += "_positions.txt";
	return outputPath;
}

std::string GetOutputProcessedImagesFilePath(std::string dataFilePath, std::string processMethodName, std::string outputFormat) {
	// remove the last 4 characters and replace with "_processed.extension"
	
	std::string outputPath = dataFilePath;
	outputPath.erase(outputPath.length() - 4, 4);
	
	if (processMethodName != std::string("convertfileformat")) {
		outputPath += "_processed.";
	} else {
		outputPath += ".";
	}
	
	if ((outputFormat == std::string("tiff")) || (outputFormat == std::string("compressedtiff")))
		outputPath += "tif";
	else if (outputFormat == std::string("pde"))
		outputPath += "pde";
	
	return outputPath;
}
