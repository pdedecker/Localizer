/*
 *  PALM_analysis_FileIO.h
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 30/03/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef PALM_ANALYSIS_FILEIO
#define PALM_ANALYSIS_FILEIO

#include <fstream>
#include <iostream>
#include <sstream>
#include <queue>
#include <list>
#include "PALM_analysis_errors.h"
#include "PALM_analysis_storage.h"
#include "PALM_analysis_defines.h"
#include "tiffio.h"
#include "stdint.h"
#include "XOPStandardHeaders.h"

using namespace std;

class XOPFileHandler {
public:
	XOPFileHandler() {fileRef = NULL;}
	~XOPFileHandler();
	
	void open(const char *path_rhs);
	void open(const char *path_rhs, std::ios_base::openmode mode) {open (path_rhs);}	// for compatibility with the standard library
	
	void close();
	int fail() {return 0;}	// we will let this class throw exceptions
	// so if a failure occurs then the code will never be able to check for fail() since control will have passed to the error-handling function
	// we keep the function anyway since it matches the standard ifstream classes so the class can be swapped in
	
	int is_open() {return (fileRef != NULL);}
	
	void get(char & c);
	void read(char *buffer, size_t nBytes);
	void getline(char *buffer, size_t nMax);
	
	uint64_t tellg();
	void seekg(uint64_t pos);
	
	
private:
	XOP_FILE_REF fileRef;
	string path;
};

uint16 getUINT16FromCharArray(char *array, size_t offset);
uint32 getUINT32FromCharArray(char *array, size_t offset);

class ImageLoader {
public:
	ImageLoader();
	ImageLoader(const string rhs);
	virtual ~ImageLoader();
	
	size_t get_total_number_of_images() const {return total_number_of_images;}
	size_t getXSize() const {return x_size;}
	size_t getYSize() const {return y_size;}
	boost::shared_ptr<PALMMatrix<double> > get_nth_image(const size_t n);	// images are numbered from 0 to N - 1
	
protected:
	virtual void parse_header_information() = 0;
	virtual void ReadImagesFromDisk(size_t const nStart, size_t const nEnd, vector<boost::shared_ptr<PALMMatrix <double> > > & cache) = 0;
	
	size_t cacheStart;
	size_t cacheEnd;
	
	size_t image_cache_size;
	
	string path;
	// ifstream file;
	XOPFileHandler file;
	uint64_t header_length;
	size_t total_number_of_images;
	uint64_t x_size;
	uint64_t y_size;
	int storage_type;
	vector <boost::shared_ptr<PALMMatrix<double> > > image_cache;	// array of pointer to gsl_matrices containing cached images
};

class ImageLoaderAndor : public ImageLoader {
public:
	ImageLoaderAndor(string rhs, size_t image_cache_size_rhs);
	~ImageLoaderAndor();
	
protected:
	void parse_header_information();
	void ReadImagesFromDisk(size_t const nStart, size_t const nEnd, vector<boost::shared_ptr<PALMMatrix <double> > > & cache);
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
	ImageLoaderHamamatsu(string rhs, size_t image_cache_size);
	~ImageLoaderHamamatsu();
	
protected:
	void parse_header_information();
	void ReadImagesFromDisk(size_t const nStart, size_t const nEnd, vector<boost::shared_ptr<PALMMatrix <double> > > & cache);
};

class ImageLoaderSPE : public ImageLoader {
public:
	ImageLoaderSPE(string rhs, size_t image_cache_size);
	~ImageLoaderSPE();
	
protected:
	void parse_header_information();
	void ReadImagesFromDisk(size_t const nStart, size_t const nEnd, vector<boost::shared_ptr<PALMMatrix <double> > > & cache);
};

class SimpleImageLoader : public ImageLoader {	// loads data from a binary file from a square array consisting of size_ts in row-major order
public:
	SimpleImageLoader(string rhs, size_t image_cache_size);
	~SimpleImageLoader();
	
protected:
	void parse_header_information();
	void ReadImagesFromDisk(size_t const nStart, size_t const nEnd, vector<boost::shared_ptr<PALMMatrix <double> > > & cache);
};

class ImageLoaderTIFF : public ImageLoader {	// loads data from TIFF files using the libtiff library
public:
	ImageLoaderTIFF(string rhs, size_t image_cache_size);
	~ImageLoaderTIFF();
	
protected:
	void parse_header_information();
	void ReadImagesFromDisk(size_t const nStart, size_t const nEnd, vector<boost::shared_ptr<PALMMatrix <double> > > & cache);
	
	TIFF* tiff_file;
	unsigned int bitsPerPixel;
	unsigned int sampleFormat;	// unsigned integer of floating point?
	// 1 for uint, 3 for floating point (same as tiff specification)
	vector<size_t> directoryIndices;	// some images in the TIFF format are not actual data, but low resolution previews of other images and such
	// this means that the image at directory i is not necessarily the i-th experimental frame
	// we correct for this using this array
	
};

class ImageLoaderIgor : public ImageLoader {
public:
	ImageLoaderIgor(string waveName, size_t image_cache_size_rhs);
	~ImageLoaderIgor() {;}
	
protected:
	void parse_header_information() {;}
	void ReadImagesFromDisk(size_t const nStart, size_t const nEnd, vector<boost::shared_ptr<PALMMatrix <double> > > & cache);
	// technically we are not reading from disk, but we keep the name since it corresponds to a virtual method in the base class
	
	waveHndl igor_data_wave;
};



class OutputWriter {
public:
	OutputWriter();
	OutputWriter(const string &rhs, int overwrite);
	virtual ~OutputWriter() {;}
	
	string get_file_path() const {return file_path;}
	size_t get_n_images_written() const {return n_images_written;}
	
	virtual void write_image(boost::shared_ptr<PALMMatrix<double> > new_image) = 0;
	
	virtual int flush_and_close() = 0;
	
protected:
	
	string file_path;
	ofstream file;
	
	size_t n_images_written;
	size_t x_size;
	size_t y_size;
	
	queue <boost::shared_ptr<PALMMatrix<double> > > image_buffer;
};


class SimpleOutputWriter : public OutputWriter {
public:
	SimpleOutputWriter(const string &rhs, int overwrite);
	~SimpleOutputWriter();
	
	void write_image(boost::shared_ptr<PALMMatrix<double> > new_image);
	
	int flush_and_close();
	
protected:
	void flush_cache();
};

class TIFFOutputWriter : public OutputWriter {
public:
	TIFFOutputWriter(const string &rhs, int overwrite);
	~TIFFOutputWriter();
	
	void write_image(boost::shared_ptr<PALMMatrix<double> > new_image);
	
	int flush_and_close();
	
protected:
	void flush_cache();
	int compress;	// if 0 then don't compress the data, otherwise compress
	
	TIFF *tiff_file;
};

class IgorOutputWriter {	// this class will take care of outputting the fitted positions to an Igor wave
public:
	IgorOutputWriter(const string rhs);
	~IgorOutputWriter();
	
	int append_new_positions(boost::shared_ptr<PALMMatrix<double> > positions);
	int write_positions_to_wave();	// to be called when the positions have been uploaded
	// writes the positions to the wave
	// the wave is only created at this point, and WILL overwrite previous waves of the same name
	
private:
	list <boost::shared_ptr<PALMMatrix<double> > > positionsList;
	size_t total_number_of_positions;
	string wave_name;
};


#endif
