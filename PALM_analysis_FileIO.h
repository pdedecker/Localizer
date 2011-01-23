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
#include "stdint.h"
#include "boost/thread.hpp"
#include <Eigen/Eigen>

#ifdef _WIN32
#include <stdio.h>
#endif

#ifdef WITH_IGOR
#include "XOPStandardHeaders.h"
#include "PALM_analysis_IgorUtilities.h"
#endif

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
	ImageLoader(const std::string rhs);
	virtual ~ImageLoader();
	
	size_t GetNImages() const {return total_number_of_images;}
	size_t getXSize() const {return x_size;}
	size_t getYSize() const {return y_size;}
	int getStorageType() const {return storage_type;}
	virtual boost::shared_ptr<Eigen::MatrixXd> readImage(const size_t index) = 0;	// images are numbered from 0 to N - 1
	
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
	size_t total_number_of_images;
	uint64_t x_size;
	uint64_t y_size;
	int storage_type;
	
	boost::mutex loadImagesMutex;	// a mutex to ensure that we don't try to load two images at once
};

class ImageLoaderSPE : public ImageLoader {
public:
	ImageLoaderSPE(std::string rhs);
	~ImageLoaderSPE();
	
	boost::shared_ptr<Eigen::MatrixXd> readImage(const size_t index);
	
protected:
	void parse_header_information();
};

class ImageLoaderAndor : public ImageLoader {
public:
	ImageLoaderAndor(std::string rhs);
	~ImageLoaderAndor();
	
	boost::shared_ptr<Eigen::MatrixXd> readImage(const size_t index);
	
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
	
	boost::shared_ptr<Eigen::MatrixXd> readImage(const size_t index);
	
protected:
	void parse_header_information();
};

class ImageLoaderPDE : public ImageLoader {	// loads data from a binary file from a square array consisting of size_ts in row-major order
public:
	ImageLoaderPDE(std::string rhs);
	~ImageLoaderPDE();
	
	boost::shared_ptr<Eigen::MatrixXd> readImage(const size_t index);
	
protected:
	void parse_header_information();
};

class ImageLoaderTIFF : public ImageLoader {	// loads data from TIFF files using the libtiff library
public:
	ImageLoaderTIFF(std::string rhs);
	~ImageLoaderTIFF();
	
	boost::shared_ptr<Eigen::MatrixXd> readImage(const size_t index);
	
protected:
	void parse_header_information();
	
	TIFF* tiff_file;
	std::vector<size_t> directoryIndices;
};

#ifdef WITH_IGOR
class ImageLoaderIgor : public ImageLoader {
public:
	ImageLoaderIgor(std::string waveName);
	~ImageLoaderIgor() {;}
	
	boost::shared_ptr<Eigen::MatrixXd> readImage(const size_t index);
	
protected:
	void parse_header_information() {;}
	
	waveHndl igor_data_wave;
};
#endif // WITH_IGOR



class ImageOutputWriter {
public:
	ImageOutputWriter();
	ImageOutputWriter(const std::string &rhs, int overwrite);
	virtual ~ImageOutputWriter() {;}
	
	std::string get_file_path() const {return file_path;}
	size_t get_n_images_written() const {return n_images_written;}
	
	virtual void write_image(boost::shared_ptr<Eigen::MatrixXd> imageToWrite) = 0;
	
protected:
	
	std::string file_path;
	std::ofstream file;
	
	size_t n_images_written;
	
	std::queue <boost::shared_ptr<Eigen::MatrixXd> > image_buffer;
};


class PDEImageOutputWriter : public ImageOutputWriter {
public:
	PDEImageOutputWriter(const std::string &rhs, int overwrite, uint32_t storageType);
	~PDEImageOutputWriter();
	
	void write_image(boost::shared_ptr<Eigen::MatrixXd> imageToWrite);
	
protected:
	void WriteHeader();
	
	uint32_t x_size, y_size;
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
	
	void write_image(boost::shared_ptr<Eigen::MatrixXd> imageToWrite);
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
	
	void write_image(boost::shared_ptr<Eigen::MatrixXd> new_image);
	
protected:
	size_t nImagesTotal;
	std::string fullPathToWave;
	DataFolderAndName waveDataFolderAndName;
	waveHndl outputWave;
	int overwrite;
	int storageType;
};
#endif // WITH_IGOR

#endif
