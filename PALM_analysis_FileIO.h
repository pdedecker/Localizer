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
	
	// int open_file(const string rhs);
	// int close_current_file();
	
	unsigned long get_total_number_of_images() const {return total_number_of_images;}
	unsigned long get_x_size() const {return x_size;}
	unsigned long get_y_size() const {return y_size;}
	virtual boost::shared_ptr<encap_gsl_matrix> get_nth_image(const unsigned long n) = 0;	// images are numbered from 0 to N - 1
	
protected:
	virtual int parse_header_information() = 0;
	
	unsigned long image_cache_size;
	
	string path;
	// ifstream file;
	XOPFileHandler file;
	uint64_t header_length;
	unsigned long total_number_of_images;
	uint64_t x_size;
	uint64_t y_size;
	int storage_type;
	vector <boost::shared_ptr<encap_gsl_matrix> > image_cache;	// array of pointer to gsl_matrices containing cached images
	vector <unsigned long> images_in_cache;	// keeps track of the indices of the images in the cache
	// if there is no image at a particular location then this is set to -1
};

class ImageLoaderAndor : public ImageLoader {
public:
	ImageLoaderAndor(string rhs);
	ImageLoaderAndor(string rhs, unsigned long image_cache_size_rhs);
	~ImageLoaderAndor();
	
	boost::shared_ptr<encap_gsl_matrix> get_nth_image(const unsigned long n);
	
protected:
	int parse_header_information();
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
	ImageLoaderHamamatsu(string rhs);
	ImageLoaderHamamatsu(string rhs, unsigned long image_cache_size);
	~ImageLoaderHamamatsu();
	
	boost::shared_ptr<encap_gsl_matrix> get_nth_image(const unsigned long n);
	
protected:
	int parse_header_information();
};

class ImageLoaderSPE : public ImageLoader {
public:
	ImageLoaderSPE(string rhs);
	ImageLoaderSPE(string rhs, unsigned long image_cache_size);
	~ImageLoaderSPE();
	
	boost::shared_ptr<encap_gsl_matrix> get_nth_image(const unsigned long n);
	
protected:
	int parse_header_information();
};

class SimpleImageLoader : public ImageLoader {	// loads data from a binary file from a square array consisting of unsigned longs in row-major order
public:
	SimpleImageLoader(string rhs);
	SimpleImageLoader(string rhs, unsigned long image_cache_size);
	~SimpleImageLoader();
	
	boost::shared_ptr<encap_gsl_matrix> get_nth_image(const unsigned long n);
	
protected:
	int parse_header_information();
};

class ImageLoaderTIFF_Igor : public ImageLoader {	// loads data from TIFF files using extensive callbacks to IGOR because the TIFF format is a mess
public:
	ImageLoaderTIFF_Igor(string rhs);
	ImageLoaderTIFF_Igor(string rhs, unsigned long image_cache_size);
	~ImageLoaderTIFF_Igor();
	
	boost::shared_ptr<encap_gsl_matrix> get_nth_image(const unsigned long n);
	
protected:
	int parse_header_information();
	
	string native_Igor_FilePath;	// in Windows format or the native HFS Macintosh format
	waveHndl tiff_cache_wave;
	
};

class ImageLoaderTIFF : public ImageLoader {	// loads data from TIFF files using the libtiff library
public:
	ImageLoaderTIFF(string rhs);
	ImageLoaderTIFF(string rhs, unsigned long image_cache_size);
	~ImageLoaderTIFF();
	
	boost::shared_ptr<encap_gsl_matrix> get_nth_image(const unsigned long n);
	
protected:
	int parse_header_information();
	
	TIFF* tiff_file;
	unsigned int bitsPerPixel;
	unsigned int sampleFormat;	// unsigned integer of floating point?
	// 1 for uint, 3 for floating point (same as tiff specification)
	vector<unsigned long> directoryIndices;	// some images in the TIFF format are not actual data, but low resolution previews of other images and such
	// this means that the image at directory i is not necessarily the i-th experimental frame
	// we correct for this using this array
	
};

class ImageLoaderIgor : public ImageLoader {
public:
	ImageLoaderIgor(string waveName);
	~ImageLoaderIgor() {;}
	
	boost::shared_ptr<encap_gsl_matrix> get_nth_image(const unsigned long n);
	
protected:
	int parse_header_information() {return 0;}
	
	waveHndl igor_data_wave;
};



class OutputWriter {
public:
	OutputWriter();
	OutputWriter(const string &rhs, int overwrite);
	virtual ~OutputWriter() {;}
	
	string get_file_path() const {return file_path;}
	unsigned long get_n_images_written() const {return n_images_written;}
	
	virtual void write_image(boost::shared_ptr<encap_gsl_matrix> new_image) = 0;
	
	virtual int flush_and_close() = 0;
	
protected:
	
	string file_path;
	ofstream file;
	
	unsigned long n_images_written;
	unsigned long x_size;
	unsigned long y_size;
	
	queue <boost::shared_ptr<encap_gsl_matrix> > image_buffer;
};


class SimpleOutputWriter : public OutputWriter {
public:
	SimpleOutputWriter(const string &rhs, int overwrite);
	~SimpleOutputWriter();
	
	void write_image(boost::shared_ptr<encap_gsl_matrix> new_image);
	
	int flush_and_close();
	
protected:
	void flush_cache();
};

class TIFFOutputWriter : public OutputWriter {
public:
	TIFFOutputWriter(const string &rhs, int overwrite);
	~TIFFOutputWriter();
	
	void write_image(boost::shared_ptr<encap_gsl_matrix> new_image);
	
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
	
	int append_new_positions(boost::shared_ptr<encap_gsl_matrix> positions);
	int write_positions_to_wave();	// to be called when the positions have been uploaded
	// writes the positions to the wave
	// the wave is only created at this point, and WILL overwrite previous waves of the same name
	
private:
	list <boost::shared_ptr<encap_gsl_matrix> > positionsList;
	unsigned long total_number_of_positions;
	string wave_name;
};


#endif