/*
 *  PALM_analysis_FileIO.h
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 30/03/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef PALM_ANALYSIS_FILEIO_H
#define PALM_ANALYSIS_FILEIO_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <queue>
#include <list>
#include <string>
#include "PALM_analysis_errors.h"
#include "PALM_analysis_storage.h"
#include "PALM_analysis_defines.h"
#include "PALM_analysis.h"
#include "tiffio.h"
#include "boost/cstdint.hpp"
#include "boost/thread.hpp"
#include <Eigen/Eigen>

#ifdef _WIN32
#include <stdio.h>
#endif

#ifdef WITH_IGOR
#include "XOPStandardHeaders.h"
#include "PALM_analysis_IgorUtilities.h"
#endif

using boost::uint64_t;
using boost::int64_t;
using boost::uint32_t;
using boost::int32_t;
using boost::uint16_t;
using boost::int16_t;
using boost::uint8_t;
using boost::int8_t;

/**
 * Provide an ifstream-like class for Windows that handles large file offset
 */
#ifdef _WIN32
class WindowsFileStream {
public:
    WindowsFileStream() {fileRef = NULL;}
    ~WindowsFileStream();
    
    void open(const char *path_rhs);
    void open(const char *path_rhs, std::ios_base::openmode mode) {open (path_rhs);}    // for compatibility with the standard library
    
    void close();
	
    int fail() {return ferror(fileRef);}
	int good() {return !ferror(fileRef);}
    
    int is_open() {return (fileRef != NULL);}
    
    void get(char & c);
    void read(char *buffer, size_t nBytes);
    void getline(char *buffer, size_t nMax);
    
    uint64_t tellg();
    void seekg(uint64_t pos);
    
    
private:
    FILE *fileRef;
	std::string path;
};
#endif

class ImageLoader {
public:
	ImageLoader();
	virtual ~ImageLoader();
	
	size_t getNImages() const {return nImages;}
	size_t getXSize() const {return xSize;}
	size_t getYSize() const {return ySize;}
	int getStorageType() const {return storage_type;}
	
	/**
	 * readImage explicitly asks for the image at a certain index, but is not reentrant.
	 */
	ImagePtr readImage(const size_t index);	// images are numbered from 0 to N - 1
	
	/*
	 * readNext asks for the next image in the sequence, and is
	 * reentrant, but throws a std::runtime exception if there are no more images
	 * in the sequence. Also, due to its reentrant nature there must be some way
	 * for the caller to know which image was returned. This is returned by reference
	 * in the argument.
	 */
	virtual ImagePtr readNextImage(size_t &index) = 0;
	/*
	 * Get the next image in the file and loop to the begin after the last image.
	 * Not reentrant.
	 */
	ImagePtr readNextImageAndLoop(size_t &index);
	/*
	 * Spool the file so the specified frame will be the one read by readNextImage
	 */
	void spoolTo(size_t index);
	/*
	 * 'Rewind' the file to the beginning so that readNextImage will return the first frame.
	 */
	void rewind() {this->spoolTo(0);}
	
protected:
	virtual void parse_header_information() = 0;
	
	/**
	 * Check if the values parsed from the header (xSize etc)
	 * are reasonable. If not throw an std::runtime error
	 */
	void checkForReasonableValues();
	
	std::string filePath;
#ifdef _WIN32
	WindowsFileStream file;
#else
	std::ifstream file;
#endif
	uint64_t header_length;
	size_t nImages;
	uint64_t xSize;
	uint64_t ySize;
	int storage_type;
	size_t nextImageToRead;
	
	boost::mutex loadImagesMutex;	// a mutex to ensure that we don't try to load two images at once
};

class ImageLoaderSPE : public ImageLoader {
public:
	ImageLoaderSPE(std::string rhs);
	~ImageLoaderSPE();
	
	ImagePtr readNextImage(size_t &index);
	
protected:
	void parse_header_information();
};

class ImageLoaderAndor : public ImageLoader {
public:
	ImageLoaderAndor(std::string rhs);
	~ImageLoaderAndor();
	
	ImagePtr readNextImage(size_t &index);
	
protected:
	void parse_header_information();
};

class ImageLoaderHamamatsu_HeaderStructure {
public:
	uint16 magic;
	uint16 commentLength;
	uint16 xSize;
	uint16 ySize;
	uint16 xBinning;	// uncertain
	uint16 yBinning;	// uncertain
	uint16 storageFormat;
	uint32 nImages;
	uint16 nChannels;
	uint16 channel;		// uncertain
	double timeStamp;
	uint32 marker;
	uint32 misc;		// function unknown
};

class ImageLoaderHamamatsu : public ImageLoader {
public:
	ImageLoaderHamamatsu(std::string rhs);
	~ImageLoaderHamamatsu();
	
	ImagePtr readNextImage(size_t &index);
	
protected:
	void parse_header_information();
};

class ImageLoaderPDE : public ImageLoader {	// loads data from a binary file from a square array consisting of size_ts in row-major order
public:
	ImageLoaderPDE(std::string rhs);
	~ImageLoaderPDE();
	
	ImagePtr readNextImage(size_t &index);
	
protected:
	void parse_header_information();
};

class ImageLoaderTIFF : public ImageLoader {	// loads data from TIFF files using the libtiff library
public:
	ImageLoaderTIFF(std::string rhs);
	~ImageLoaderTIFF();
	
	ImagePtr readNextImage(size_t &index);
	
	/*
	 * The implementation of spoolTo in the base ImageLoader class
	 * is overridden due to the linked list nature of TIFF files
	 */
	void spoolTo(size_t index);
	
protected:
	void parse_header_information();
	
	TIFF* tiff_file;
	std::vector<size_t> directoryIndices;
	
	size_t currentDirectoryIndex;
};

#ifdef WITH_IGOR
class ImageLoaderIgor : public ImageLoader {
public:
	ImageLoaderIgor(std::string waveName);
	~ImageLoaderIgor() {;}
	
	ImagePtr readNextImage(size_t &index);
	
protected:
	void parse_header_information() {;}
	
	waveHndl igor_data_wave;
};
#endif // WITH_IGOR



class ImageOutputWriter {
public:
	ImageOutputWriter();
	virtual ~ImageOutputWriter() {;}
	
	std::string getOutputFilePath() const {return outputFilePath;}
	size_t getNImagesWritten() const {return nImagesWritten;}
	
	virtual void write_image(ImagePtr imageToWrite) = 0;
	
protected:
	
	std::string outputFilePath;
	std::ofstream file;
	
	size_t nImagesWritten;
};


class PDEImageOutputWriter : public ImageOutputWriter {
public:
	PDEImageOutputWriter(const std::string &rhs, int overwrite, uint32_t storageType);
	~PDEImageOutputWriter();
	
	void write_image(ImagePtr imageToWrite);
	
protected:
	void WriteHeader();
	
	uint32_t xSize, ySize;
	uint32_t storageType;
};

struct PDEFormatHeader {
	uint32_t magic;
	uint32_t version;
	uint32_t nImages;
	uint32_t xSize;
	uint32_t ySize;
	uint32_t storageFormat;
};
typedef struct PDEFormatHeader PDEFormatHeader;


class TIFFImageOutputWriter : public ImageOutputWriter {
public:
	TIFFImageOutputWriter(const std::string &rhs, int overwrite, int compression_rhs, int storageType);
	~TIFFImageOutputWriter();
	
	void write_image(ImagePtr imageToWrite);
protected:
	int compression;	// if 1 then don't compress the data, otherwise compress
	int storageType;
	
	TIFF *tiff_file;
};

#ifdef WITH_IGOR
class IgorImageOutputWriter : public ImageOutputWriter {
public:
	// constructor when the wave is specified using a fully qualified path or just a single name
	IgorImageOutputWriter(std::string waveName, size_t nImagesTotal, int overwrite, int storageType);
	// constructor when using the DataFolderAndName type provided by the XOP toolkit
	IgorImageOutputWriter(DataFolderAndName outputDataFolderAndName, size_t nImagesTotal, int overwrite, int storageType);
	
	~IgorImageOutputWriter() {;}
	
	void write_image(ImagePtr new_image);
	
protected:
	int GetIgorStorageType();
	
	size_t nImagesTotal;
	std::string fullPathToWave;
	DataFolderAndName waveDataFolderAndName;
	waveHndl outputWave;
	int overwrite;
	int storageType;
};
#endif // WITH_IGOR

#endif
