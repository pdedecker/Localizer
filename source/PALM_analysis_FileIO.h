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
 with analysis programs such as Igor Pro or Matlab, or to acquire
 data from or control hardware related to an experimental measurement,
 the licensors of this Program grant you additional permission 
 to convey the resulting work.
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
#include "tiffio.h"
#include "boost/cstdint.hpp"
#include "boost/thread.hpp"

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
 Function that writes the contents of a char* buffer containing a single image
 (i.e. as read from a file) to an Image.
 */
template <typename T>
void WriteBufferToImage(char *buffer, ImagePtr imagePtr) {
    T *bufferPtr = reinterpret_cast<T *>(buffer);
    size_t nRows = imagePtr->rows();
    size_t nCols = imagePtr->cols();
    for (size_t j  = 0; j < nCols; j++) {
        for (size_t i = 0; i < nRows; i++) {
            (*imagePtr)(i, j) = static_cast<double>(*bufferPtr);
            ++bufferPtr;
        }
    }
}

/**
 Function that writes the contents of a single image to a char * buffer,
 while converting to the requested number type
 */
template <typename T>
void WriteImageToBuffer(ImagePtr imagePtr, char *buffer) {
	T* bufferPtr = reinterpret_cast<T *>(buffer);
	size_t nRows = imagePtr->rows();
    size_t nCols = imagePtr->cols();
    for (size_t j  = 0; j < nCols; j++) {
        for (size_t i = 0; i < nRows; i++) {
            *bufferPtr = static_cast<T>((*imagePtr)(i, j));
            ++bufferPtr;
        }
    }
}

/**
 Provides a replacement for an fstream class, since the standard fstream classes
 in win32 do not handle file offsets larger than 2 GB.
 The replacement is not complete; in particular the arguments to open, seekg, and seekp
 are ignored. Also, there is no separate concept of a seek and put pointer, instead
 there is only a single one.
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
    
    void write(char *buffer, size_t nBytes);
    
    uint64_t tellg();
    void seekg(uint64_t pos);
    void seekp(uint64_t pos, std::ios_base::seekdir dir);
    
    
private:
    FILE *fileRef;
	std::string path;
};
#endif // _WIN32

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
	 * readNextImage asks for the next image in the sequence, and is
	 * reentrant, but throws a std::runtime exception if there are no more images
	 * in the sequence. Also, due to its reentrant nature there must be some way
	 * for the caller to know which image was returned. This is returned by reference
	 * in the argument. The non-reentrant version below does not require this argument.
	 */
	virtual ImagePtr readNextImage(size_t &index) = 0;
	ImagePtr readNextImage();
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
	uint64_t nextImageToRead;
	
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

class ImageLoaderPDE : public ImageLoader {
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
#ifdef _WIN32
	WindowsFileStream file;
#else
	std::ofstream file;
#endif
	
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
