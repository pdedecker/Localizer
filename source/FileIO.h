/*
 Copyright 2008-2014 Peter Dedecker.
 
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
#include <map>
#include <utility>
#include <cassert>
#include "Errors.h"
#include "Storage.h"
#include "Defines.h"
#include "tiffio.h"
#include "boost/cstdint.hpp"
#include "boost/thread.hpp"

#ifdef _WIN32
#include <stdio.h>
#endif

#ifdef WITH_IGOR
#include "XOPStandardHeaders.h"
#include "IgorUtilities.h"
#endif

#ifdef WITH_MATLAB
#include "mex.h"
#endif

#ifdef _WIN32
#define OFSTREAM_T WindowsFileStream
#else
#define OFSTREAM_T std::ofstream
#endif

using boost::uint64_t;
using boost::int64_t;
using boost::uint32_t;
using boost::int32_t;
using boost::uint16_t;
using boost::int16_t;
using boost::uint8_t;
using boost::int8_t;

class ImageLoader;
std::shared_ptr<ImageLoader> GetImageLoader(const std::string& data_file_path, int cameraType = -1);

/**
 Returns last modification time of the file pointed to by path, in some non-portable
 but internally consistent int format.
*/
int64_t GetLastModificationTime(const std::string& path);

ImagePtr BufferWithFormatToImage(const std::vector<char>& imageBuffer, int nRows, int nCols, int format, int treatAsRowMajor = 0);
void ImageToBufferWithFormat(ImagePtr image, int format, std::vector<char>& imageBuffer, int treatAsRowMajor = 0);
size_t NBytesInImage(int nRows, int nCols, int format);

/**
 Function that writes the contents of a std::vector<char> containing a single image
 (i.e. as read from a file) to an Image.
 */
template <typename T>
void CopyBufferToImage(const std::vector<char>& buffer, ImagePtr imagePtr, int treatAsRowMajor = 0) {
    size_t nRows = imagePtr->rows();
    size_t nCols = imagePtr->cols();
    size_t nBytesRequired = nRows * nCols * sizeof(T);
    assert(nBytesRequired == buffer.size());
    const T* bufferPtr = reinterpret_cast<const T*>(buffer.data());
	
	if (treatAsRowMajor == 0) {
		for (size_t j  = 0; j < nCols; j++) {
			for (size_t i = 0; i < nRows; i++) {
				(*imagePtr)(i, j) = static_cast<double>(*bufferPtr);
				++bufferPtr;
			}
		}
	} else {
		for (size_t i = 0; i < nRows; i++) {
			for (size_t j  = 0; j < nCols; j++) {
				(*imagePtr)(i, j) = static_cast<double>(*bufferPtr);
				++bufferPtr;
			}
		}
	}
}

/**
 Function that writes the contents of a single image to a char * buffer,
 while converting to the requested number type. Buffer will be resized to
 fit the data if needed.
 */
template <typename T>
void CopyImageToBuffer(ImagePtr imagePtr, std::vector<char>& buffer, int treatAsRowMajor = 0) {
	size_t nRows = imagePtr->rows();
    size_t nCols = imagePtr->cols();
    size_t nBytesRequired = nRows * nCols * sizeof(T);
    if (buffer.size() != nBytesRequired)
        buffer.resize(nBytesRequired);
    T* bufferPtr = reinterpret_cast<T*>(buffer.data());
	
	if (treatAsRowMajor == 0) {
		for (size_t j  = 0; j < nCols; j++) {
			for (size_t i = 0; i < nRows; i++) {
				*bufferPtr = static_cast<T>((*imagePtr)(i, j));
				++bufferPtr;
			}
		}
	} else {
		for (size_t i = 0; i < nRows; i++) {
			for (size_t j  = 0; j < nCols; j++) {
				*bufferPtr = static_cast<T>((*imagePtr)(i, j));
				++bufferPtr;
			}
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
    uint64_t tellp();
    void seekg(uint64_t pos, std::ios_base::seekdir dir = std::ios_base::beg);
    void seekp(uint64_t pos, std::ios_base::seekdir dir = std::ios_base::beg);
    
    
private:
    FILE *fileRef;
	std::string path;
};
#endif // _WIN32

template <typename T> void WriteBinaryValue(OFSTREAM_T& file, T value) {
    char* valPtr = reinterpret_cast<char*>(&value);
    file.write(valPtr, sizeof(T));
}

class ImageLoader {
public:
	ImageLoader();
	virtual ~ImageLoader();
	
	virtual int getNImages() const {return nImages;}
	virtual int getXSize() const {return xSize;}
	virtual int getYSize() const {return ySize;}
	virtual int getStorageType() const {return storage_type;}
	virtual int getFileType() const = 0;
	
	/**
	 * readImage explicitly asks for the image at a certain index, but is not reentrant.
	 */
	ImagePtr readImage(const int index);	// images are numbered from 0 to N - 1
	
	/*
	 * readNextImage asks for the next image in the sequence, and is
	 * reentrant, but throws a std::runtime exception if there are no more images
	 * in the sequence. Also, due to its reentrant nature there must be some way
	 * for the caller to know which image was returned. This is returned by reference
	 * in the argument. The non-reentrant version below does not require this argument.
	 */
	virtual ImagePtr readNextImage(int &indexOfImageThatWasRead) = 0;
	ImagePtr readNextImage();
	/*
	 * Get the next image in the file and loop to the begin after the last image.
	 * Not reentrant.
	 */
	ImagePtr readNextImageAndLoop(int &indexOfImageThatWasRead);
	/*
	 * Spool the file so the specified frame will be the one read by readNextImage
	 */
	virtual void spoolTo(int index);
	/*
	 * 'Rewind' the file to the beginning so that readNextImage will return the first frame.
	 */
	void rewind() {this->spoolTo(0);}
	
protected:
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
	int nImages;
	int xSize;
	int ySize;
	int storage_type;
	int nextImageToRead;
	
	boost::mutex loadImagesMutex;	// a mutex to ensure that we don't try to load two images at once
};

/**
 * Wraps around an existing ImageLoader, allowing it to appear as an ImageLoader with a lower number of images
 * or a cropped ROI.
 */
class ImageLoaderWrapper : public ImageLoader {
public:
    ImageLoaderWrapper(std::shared_ptr<ImageLoader> baseImageLoader);
    ~ImageLoaderWrapper() {;}
    
    virtual int getNImages() const;
    virtual int getXSize() const;
    virtual int getYSize() const;
    virtual int getStorageType() const;
	virtual int getFileType() const;
    
    virtual ImagePtr readNextImage(int &indexOfImageThatWasRead);
    virtual void spoolTo(int index);
    
    void setImageRange(int firstImageToInclude, int lastImageToInclude);
    void setROI(int minX, int maxX, int minY, int maxY);
    
private:
    std::shared_ptr<ImageLoader> _baseImageLoader;
    int _firstImageToInclude;
    int _lastImageToInclude;
    int _minX, _maxX;
    int _minY, _maxY;
    bool _haveCustomROI;
};

class ImageLoaderSPE : public ImageLoader {
public:
	ImageLoaderSPE(std::string rhs);
	~ImageLoaderSPE();
	
	ImagePtr readNextImage(int &indexOfImageThatWasRead);
	
	int getFileType() const {return CAMERA_TYPE_WINSPEC;}
	
private:
	void parse_header_information();
};

class ImageLoaderAndor : public ImageLoader {
public:
	ImageLoaderAndor(std::string rhs);
	~ImageLoaderAndor();
	
	ImagePtr readNextImage(int &indexOfImageThatWasRead);
	
	int getFileType() const {return CAMERA_TYPE_ANDOR;}
	
private:
	void parse_header_information();
};

class ImageLoaderHamamatsu : public ImageLoader {
public:
	class ImageOffsets {
	public:
		std::vector<uint64_t> offsets;
		int64_t modificationTime;
		uint64_t xSize;
		uint64_t ySize;
	};
	
	ImageLoaderHamamatsu(std::string rhs);
	~ImageLoaderHamamatsu();
	
	ImagePtr readNextImage(int &indexOfImageThatWasRead);
	
	int getFileType() const {return CAMERA_TYPE_HAMAMATSU;}
	
private:
	void parse_header_information();
	
	static std::map<std::string, ImageLoaderHamamatsu::ImageOffsets> _offsetsMap;
	std::vector<uint64_t> _offsets;
};

class ImageLoaderPDE : public ImageLoader {
public:
	ImageLoaderPDE(std::string rhs);
	~ImageLoaderPDE();
	
	ImagePtr readNextImage(int &indexOfImageThatWasRead);
	
	int getFileType() const {return CAMERA_TYPE_PDE;}
	
private:
	void parse_header_information();
};

class ImageLoaderTIFF : public ImageLoader {	// loads data from TIFF files using the libtiff library
public:
    class ImageOffsets {
    public:
		std::vector<uint64_t> offsets;
		int64_t modificationTime;
    };
    
	ImageLoaderTIFF(std::string rhs);
	~ImageLoaderTIFF();
	
	ImagePtr readNextImage(int &indexOfImageThatWasRead);
	
	int getFileType() const {return CAMERA_TYPE_TIFF;}
	
private:
	void parse_header_information();
    void _extractSampleFormat();
	
	TIFF* tiff_file;
    static std::map<std::string, ImageLoaderTIFF::ImageOffsets> _offsetsMap;
	std::vector<uint64_t> _directoryOffsets;
};

class ImageLoaderMultiFileTIFF : public ImageLoader {	// loads data from TIFF files using the libtiff library
public:
	ImageLoaderMultiFileTIFF(std::string rhs);
	~ImageLoaderMultiFileTIFF() {;}
	
	ImagePtr readNextImage(int &indexOfImageThatWasRead);
	
	int getFileType() const {return CAMERA_TYPE_MULTIFILE_TIFF;}
	
private:
	std::string getFilePathForImageAtIndex(int index);
	bool imageFileAtIndexExists(int index);
	std::pair<int, int> findFirstAndLastValidImageIndices(int knownValidImageIndex);
	
	std::string baseFilePath;
	std::string extension;	// includes the '.'
	int nDigitsInNumber;
	int firstImageIndex;
};

#ifdef WITH_IGOR
class ImageLoaderIgor : public ImageLoader {
public:
	ImageLoaderIgor(std::string waveName);
	~ImageLoaderIgor() {;}
	
	ImagePtr readNextImage(int &indexOfImageThatWasRead);
	
	int getFileType() const {return CAMERA_TYPE_IGOR_WAVE;}
	
private:
    waveHndl igor_data_wave;
};
#endif // WITH_IGOR

#ifdef WITH_MATLAB
class ImageLoaderMatlab : public ImageLoader {
public:
	ImageLoaderMatlab(mxArray* matlabArray);
	~ImageLoaderMatlab() {;}
	
	ImagePtr readNextImage(int &indexOfImageThatWasRead);
	
	int getFileType() const {return CAMERA_TYPE_MATLAB_MATRIX;}
	
private:
    mxArray* _matlabArray;
};
#endif // WITH_MATLAB

class ImageOutputWriter {
public:
	ImageOutputWriter();
	virtual ~ImageOutputWriter() {;}
	
	std::string getOutputFilePath() const {return outputFilePath;}
	size_t getNImagesWritten() const {return nImagesWritten;}
	
	virtual void write_image(ImagePtr imageToWrite) = 0;
	
protected:
	
	std::string outputFilePath;
	OFSTREAM_T file;
	
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

void TIFFSampleFormatAndBitsPerSampleForFormat(const int dataFormat, int& sampleFormat, int& bitsPerSample);

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

class MultiFileTIFFImageOutputWriter : public ImageOutputWriter {
public:
	MultiFileTIFFImageOutputWriter(const std::string &baseOutputFilePath_rhs, int overwrite_rhs, int compression_rhs, int storageType_rhs);
	// baseOutputFilePath must be the full path to the output base name
	// so if we want files such as /folder/base0000.tif then baseOutputFilePath is /folder/base
	
	void write_image(ImagePtr imageToWrite);
protected:
	std::string baseOutputFilePath;
	int overwrite;
	int compression;	// if 1 then don't compress the data, otherwise compress
	int storageType;
};

class LocalizerTIFFImageOutputWriter : public ImageOutputWriter {
public:
    class TIFFIFDOnDisk {
    public:
        uint64_t ifdOffset;
        int nRows;
        int nCols;
        uint64_t dataOffset;
        uint64_t nextIFDFieldOffset;
    };
    
    LocalizerTIFFImageOutputWriter(const std::string &rhs, int overwrite, int compression, int storageType);
	~LocalizerTIFFImageOutputWriter();
	
	void write_image(ImagePtr imageToWrite);
    
private:
    std::pair<LocalizerTIFFImageOutputWriter::TIFFIFDOnDisk, std::vector<char> > _constructIFD(ImagePtr image, uint64_t ifdWillBeAtThisOffset, bool isBigTiff, bool reuseExistingData = false, LocalizerTIFFImageOutputWriter::TIFFIFDOnDisk existingIFDOnDisk = LocalizerTIFFImageOutputWriter::TIFFIFDOnDisk()) const;
    void _writeTag(char*& bufferPtr, int tagID, uint64_t count, uint64_t value, bool isBigTiff) const;
    void _convertToBigTiff();
    template <typename T> void _storeInBuffer(char*& bufferPtr, T value) const {
        void *valuePtr = reinterpret_cast<void*>(&value);
        memcpy(bufferPtr, valuePtr, sizeof(T));
        bufferPtr += sizeof(T);
    }
    void _writeTiffHeader();
    void _writeBigTiffHeader();
    void _touchupOffsets();
    
    std::vector<TIFFIFDOnDisk> _writtenIFDs;
    int _storageType;
    bool _isBigTiff;
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
    waveHndl getWave() const {return outputWave;}
	
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
