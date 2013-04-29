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

#include "PALM_analysis_MatrixRecycler.h"
#include "PALM_analysis_FileIO.h"

#include <sys/types.h>
#include <sys/stat.h>

int64_t GetLastModificationTime(const std::string& path) {
	int64_t modTime = 1;
	
#ifdef _WIN32
	struct _stat buffer;
	int err = _stat(path.c_str(), &buffer);
	if ((err != 0) && (errno != 0)) {
		// for whatever reason _stat() will return -1 but not set errno,
		// which leads me to believe that this -1 value does not indicate an error.
		// the microsoft docs say that _stat() should return non-zero for error, and should
		// set errno as well. So I'm seeing undocumented behavior.
		char errorMessage[100];
		sprintf(errorMessage, "Return value %d from _stat(), errno is %d", err, errno);
		if (errno != 0) {
			throw std::runtime_error(errorMessage);
		}
	}
	
	modTime = static_cast<int64_t>(buffer.st_mtime);
	
#else	// not _WIN32
	struct stat buffer;
	int err = stat(path.c_str(), &buffer);
	if (err != 0)
		throw std::runtime_error("non-zero return from _stat()");
	
	modTime = static_cast<int64_t>(buffer.st_mtime);
	
#endif	// _WIN32
	
	return modTime;
}

#ifdef _WIN32
WindowsFileStream::~WindowsFileStream() {
    if (this->fileRef != NULL) {
        fclose(this->fileRef);
        this->fileRef = NULL;
    }
}

void WindowsFileStream::open(const char* path_rhs) {
    int err;
    
	if (this->fileRef != NULL) {
		fclose(this->fileRef);
		this->fileRef = NULL;
	}
    err = fopen_s(&fileRef, path_rhs, "rb");
    if (err != 0) {
		std::string error;
		std::stringstream ss;
        ss << "Error " << err << " returned using open() on the image file at \"" << path_rhs << "\"";
		error = ss.str();
        throw ERROR_READING_FILE_DATA(error);
    }
    
    path = path_rhs;
}

void WindowsFileStream::close() {
    int err;
    if (this->fileRef != NULL) {
        err = fclose(this->fileRef);
        if (err != 0) {
			std::string error;
			std::stringstream ss;
			ss << "Error " << err << " returned using close() on the image file at \"" << path << "\"";
            error = ss.str();
            throw ERROR_READING_FILE_DATA(error);
        }
        this->fileRef = NULL;
    }
}

void WindowsFileStream::get(char& c) {
	assert (this->fileRef != NULL);
	if (this->fileRef == NULL) {
		throw ERROR_READING_FILE_DATA(std::string("\"get\" was called on a NULL FILE*"));
	}
    c = fgetc(this->fileRef);
}

void WindowsFileStream::read(char *buffer, size_t nBytes) {
	assert (this->fileRef != NULL);
    int itemsRead;
    itemsRead = fread(buffer, nBytes, 1, this->fileRef);
    if (itemsRead != 1) {	// signals that not everything was read
		std::string error;
		error = "Could not read a full buffer from the file at \"";
		error += path;
		error += "\"";
        throw ERROR_READING_FILE_DATA(error);
    }
}

void WindowsFileStream::getline(char *buffer, size_t nMax) {
	assert (this->fileRef != NULL);
    fgets(buffer, nMax, this->fileRef);
    if (ferror(fileRef) != 0) {
		std::string error;
		error = "Error returned using getline() on the image file at \"";
		error += path;
		error += "\"";
        throw ERROR_READING_FILE_DATA(error);
    }
}

void WindowsFileStream::write(char *buffer, size_t nBytes) {
    assert (this->fileRef != NULL);
    
    int err = fwrite(buffer, nBytes, 1, this->fileRef);
    if (ferror(fileRef) != 0) {
		std::string error;
		error = "Error returned using write() on the image file at \"";
		error += path;
		error += "\"";
        throw ERROR_READING_FILE_DATA(error);
    }
}

uint64_t WindowsFileStream::tellg() {
	assert (this->fileRef != NULL);
    uint64_t pos;
    
    pos = _ftelli64(this->fileRef);
    return pos;
}

void WindowsFileStream::seekg(uint64_t pos) {
	assert (this->fileRef != NULL);
    int err;
    
    err = _fseeki64(this->fileRef, pos, SEEK_SET);
    if (err != 0) {
		std::string error;
       	error = "Error returned using seekg() on the image file at \"";
		error += path;
		error += "\"";
        throw ERROR_READING_FILE_DATA(error);
    }
}

void WindowsFileStream::seekp(uint64_t pos, std::ios_base::seekdir dir) {
    assert (this->fileRef != NULL);
    int err;
    
    err = _fseeki64(this->fileRef, pos, SEEK_SET);
    if (err != 0) {
		std::string error;
       	error = "Error returned using seekg() on the image file at \"";
		error += path;
		error += "\"";
        throw ERROR_READING_FILE_DATA(error);
    }
}
#endif // _WIN32


ImageLoader::ImageLoader() {
	this->nImages = 0;
	this->xSize = 0;
	this->ySize = 0;
	this->header_length = 0;
	this->nextImageToRead = 0;
}

ImageLoader::~ImageLoader() {
	if (file.is_open() == 1)
		file.close();
}

void ImageLoader::checkForReasonableValues() {
	if ((this->xSize > kMaxImageDimension) || (this->ySize > kMaxImageDimension)) {
		throw std::runtime_error("the reported frame size is unreasonably large");
	}
	
	if (this->nImages > kMaxNFrames) {
		throw std::runtime_error("the reported number of frames is unreasonably large");
	}
}

ImagePtr ImageLoader::readImage(size_t index) {
	this->spoolTo(index);
	return this->readNextImage(index);
}

ImagePtr ImageLoader::readNextImage() {
	size_t dummy;
	return this->readNextImage(dummy);
}

ImagePtr ImageLoader::readNextImageAndLoop(size_t &index) {
	this->spoolTo(index % this->nImages);
	return this->readNextImage(index);
}

void ImageLoader::spoolTo(size_t index) {
	if (index >= nImages)
		throw IMAGE_INDEX_BEYOND_N_IMAGES(std::string("Requested more images than there are in the file"));
	
	// don't do any work if it's not required
	if (this->nextImageToRead != index)
		this->nextImageToRead = index;
}

ImageLoaderSPE::ImageLoaderSPE(std::string rhs) {
	this->filePath = rhs;
	
	header_length = 4100;
	
	file.open(this->filePath.c_str(), std::ios::binary | std::ios::in);
	if (file.fail() == 1) {
		std::string error ("Unable to open the file at ");
		error += this->filePath;
		throw CANNOT_OPEN_FILE(error);
	}
	
	parse_header_information();
}

ImageLoaderSPE::~ImageLoaderSPE() {
	if (file.is_open() == 1)
		file.close();
}

void ImageLoaderSPE::parse_header_information() {
	// warning: only safe as long as both the writing
	// and reading systems are little-endian
	this->xSize = 0;
	this->ySize = 0;
	this->nImages = 0;
	this->storage_type = 0;
	
	file.seekg(42);
	file.read((char *)&(this->xSize), 2);
	
	file.seekg(656);
	file.read((char *)&(this->ySize), 2);
	
	file.seekg(1446);
	file.read((char *)&(this->nImages), 4);
	
	file.seekg(108);
	file.read((char *)&(this->storage_type), 2);
	
	switch (storage_type) {
		case 0:
			storage_type = STORAGE_TYPE_FP32;
			break;
		case 1:
			storage_type = STORAGE_TYPE_UINT32;
			break;
		case 2:
			storage_type = STORAGE_TYPE_INT16;
			break;
		case 3:
			storage_type = STORAGE_TYPE_UINT16;
			break;
		default:
			std::string error("Unable to determine the storage type used in ");
			error += this->filePath;
			throw CANNOT_DETERMINE_SPE_STORAGE_TYPE(error);
			break;
	}
	
	// was there an error sometime during this procedure that would have caused the reading to fail?
	if (file.fail() != 0) {
		std::string error;
		error = "Error parsing the header information in \"";
		error += this->filePath;
		error += "\" assuming the SPE format";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	file.seekg(0);
	this->checkForReasonableValues();
}

ImagePtr ImageLoaderSPE::readNextImage(size_t &indexOfImageThatWasRead) {
	uint64_t offset;
	
	uint64_t imageSize, pixelSize;
	
	// determine how big we have to make the single image buffer and the offset
	switch(storage_type) {
		case STORAGE_TYPE_FP32:	// 4 byte float
		case STORAGE_TYPE_UINT32:	// 4-byte long
			pixelSize = 4;
			break;
		case STORAGE_TYPE_INT16:	// 2 byte signed short
		case STORAGE_TYPE_UINT16:	// 2 byte unsigned short
			pixelSize = 2;
			break;
		default:
			std::string error("Unable to determine the storage type used in ");
			error += this->filePath;
			throw CANNOT_DETERMINE_SPE_STORAGE_TYPE(error);
			break;
	}
	
	imageSize = pixelSize * xSize * ySize;
	
	boost::scoped_array<char> single_image_buffer(new char[imageSize]);
	ImagePtr image(GetRecycledMatrix((int)xSize, (int)ySize), FreeRecycledMatrix);
	
	{
		boost::lock_guard<boost::mutex> locker(loadImagesMutex);
		if (this->nextImageToRead >= this->nImages)
			throw std::runtime_error("requested more images than there are in the file");
		
		offset = this->header_length + this->nextImageToRead * imageSize;
		
		file.seekg(offset);
		file.read((char *)single_image_buffer.get(), imageSize);
		if (file.fail() != 0) {
			std::string error;
			error = "Error trying to read image data from \"";
			error += this->filePath;
			error += "\" assuming the SPE format";
			throw ERROR_READING_FILE_DATA(error);
		}
		indexOfImageThatWasRead = this->nextImageToRead;
		this->nextImageToRead += 1;
	}
	
	switch(storage_type) {
		case STORAGE_TYPE_FP32:	// 4-byte float
            CopyBufferToImage<float>(single_image_buffer.get(), image);
			break;
            
		case STORAGE_TYPE_UINT32:	// 4-byte long
			CopyBufferToImage<uint32_t>(single_image_buffer.get(), image);
			break;
            
		case STORAGE_TYPE_INT16:	// 2-byte signed short
			CopyBufferToImage<int16_t>(single_image_buffer.get(), image);
			break;
			
		case STORAGE_TYPE_UINT16: // 2-byte unsigned short
			CopyBufferToImage<uint16_t>(single_image_buffer.get(), image);
			break;
            
		default:
			std::string error("Unable to determine the storage type used in ");
			error += this->filePath;
			throw CANNOT_DETERMINE_SPE_STORAGE_TYPE(error);
			break;
	}
	
	return image;
}

ImageLoaderAndor::ImageLoaderAndor(std::string rhs) {
	this->filePath = rhs;
	
	file.open(this->filePath.c_str(), std::ios::binary | std::ios::in);
	if (file.fail() == 1) {
		std::string error ("Error opening the file at ");
		error += this->filePath;
		throw CANNOT_OPEN_FILE(error);
	}
	
	parse_header_information();
}

ImageLoaderAndor::~ImageLoaderAndor() {
	if (file.is_open() == 1) {
		file.close();
	}
}

void ImageLoaderAndor::parse_header_information() {
	size_t temp;
	size_t xBinning, yBinning;
	size_t frameXSize, frameYSize;
	size_t xStart, yStart, xEnd, yEnd;
	int result;
	boost::scoped_array<char> headerBuffer(new char[60001]);
	boost::scoped_array<char> singleLineBuffer(new char[4096]);
	std::string singleLine;
	
	storage_type = STORAGE_TYPE_FP32;
	
	this->file.read(headerBuffer.get(), 60001);
	if (this->file.good() != 1)
		throw std::runtime_error(std::string("Error encountered assuming the Andor format on the file at ") + this->filePath);
	
	headerBuffer[60000] = '\0';
	// look for and replace any intermediate nul characters
	for (int i = 0; i < 60000; ++i) {	// make sure that there are no intermediate NULL characters
		if (headerBuffer[i] == '\0') {
			headerBuffer[i] = '1';
		}
	}
	
	std::stringstream ss(headerBuffer.get(), std::ios::in);
	ss.getline(singleLineBuffer.get(), 4096);
	singleLine = singleLineBuffer.get();
	
	if (singleLine.find("Andor Technology Multi-Channel File") == std::string::npos) {
		throw std::runtime_error(std::string("the file at ") + this->filePath + "does not appear to be an Andor data file");
	}
	
	for (size_t i = 0;; ++i) {
		ss.getline(singleLineBuffer.get(), 4096);
		singleLine = singleLineBuffer.get();
		if ((ss.eof() == 1) || (ss.fail() == 1))
			throw std::runtime_error(std::string("premature end-of-file encountered assuming the Andor format on the file at ") + this->filePath);
		if (singleLine.find("Pixel number65") != std::string::npos) {
			// this is the first line containing info required to read the data
			break;
		}
	}
	
	// as usual Visual Studio is failing for whatever reason on these conversions
	// (no, a %Iu format specifier doesn't work either)
	// so we'll have to handle them the ugly way. Thanks again, Microsoft!
#ifdef WIN32
	int wTemp, wFrameYSize, wFrameXSize, WNImages;
	result = sscanf(singleLine.c_str(), "Pixel number%d 1 %d %d 1 %d", &wTemp, &wFrameYSize, &wFrameXSize, &WNImages);
	frameYSize = wFrameYSize; frameXSize = wFrameXSize; this->nImages = WNImages;
#else
	result = sscanf(singleLine.c_str(), "Pixel number%zu 1 %zu %zu 1 %zu", &temp, &frameYSize, &frameXSize, &this->nImages);
#endif
	if (result != 4)
		throw std::runtime_error(std::string("an error occured parsing the file assuming the Andor format"));
	
	ss.getline(singleLineBuffer.get(), 4096);
	singleLine = singleLineBuffer.get();
	
#ifdef WIN32
	int wXStart, wYStart, wYEnd, wXEnd, wXBinning, wYBinning;
	result = sscanf(singleLine.c_str(), "%d %d %d %d %d %d %d", &wTemp, &wXStart, &wYEnd, &wXEnd, &wYStart, &wXBinning, &wYBinning);
	xStart = wXStart; yEnd = wYEnd; xEnd = wXEnd; yStart = wYStart; xBinning = wXBinning; yBinning = wYBinning;
#else
	result = sscanf(singleLine.c_str(), "%zu %zu %zu %zu %zu %zu %zu", &temp, &xStart, &yEnd, &xEnd, &yStart, &xBinning, &yBinning);
#endif
	if (result != 7)
		throw std::runtime_error(std::string("an error occured parsing the file assuming the Andor format"));
	
	this->xSize = (xEnd - xStart + 1) / xBinning;
	this->ySize = (yEnd - yStart + 1) / yBinning;	// integer division
	
	// now there are some lines that may contain timestamps. There are as many lines as there are images
	for (size_t i = 0;  i < nImages; i++) {
		ss.getline(singleLineBuffer.get(), 4096);
	}
	
	this->header_length = ss.tellg();
	
	// the line after this may be a single line containing just a zero without any spaces
	// in that case add it to the header offset
	ss.getline(singleLineBuffer.get(), 4096);
	singleLine = singleLineBuffer.get();
	if (singleLine == "0") {
		this->header_length = ss.tellg();
	}
	
	// did some error happen while reading the file?
	if (ss.bad() == 1) {
		throw ERROR_READING_FILE_DATA(std::string("Error parsing the header information in \"") + this->filePath + "\" assuming the Andor format");
	}
	
	this->checkForReasonableValues();
}

ImagePtr ImageLoaderAndor::readNextImage(size_t &indexOfImageThatWasRead) {
	uint64_t offset;
	int nBytesInImage = xSize * ySize * sizeof(float);
	
	boost::scoped_array<char> single_image_buffer(new char[nBytesInImage]);
	ImagePtr image (GetRecycledMatrix((int)xSize, (int)ySize), FreeRecycledMatrix);
	
	{
		boost::lock_guard<boost::mutex> locker(loadImagesMutex);
		if (this->nextImageToRead >= nImages)
			throw IMAGE_INDEX_BEYOND_N_IMAGES(std::string("Requested more images than there are in the file"));
		
		offset = this->header_length + this->nextImageToRead * nBytesInImage;
		
		file.seekg(offset);
		file.read((char *)single_image_buffer.get(), nBytesInImage);
		if (file.fail() != 0) {
			std::string error;
			error = "Error trying to read image data from \"";
			error += this->filePath;
			error += "\" assuming the Andor format";
			throw ERROR_READING_FILE_DATA(error);
		}
		indexOfImageThatWasRead = this->nextImageToRead;
		this->nextImageToRead += 1;
	}
	
	CopyBufferToImage<float>(single_image_buffer.get(), image);
	
	return image;
}

std::map<std::string, ImageLoaderHamamatsu::ImageOffsets> ImageLoaderHamamatsu::_offsetsMap;

ImageLoaderHamamatsu::ImageLoaderHamamatsu(std::string rhs) {
	this->filePath = rhs;
	
	file.open(this->filePath.c_str(), std::ios::binary | std::ios::in);
	if (file.fail() == 1) {
		std::string error ("Unable to open the file at ");
		error += this->filePath;
		throw CANNOT_OPEN_FILE(error);
	}
	
	parse_header_information();
}

ImageLoaderHamamatsu::~ImageLoaderHamamatsu() {
	if (file.is_open() == 1)
		file.close();
}


void ImageLoaderHamamatsu::parse_header_information() {
	
	// do we already have information on the offsets of the images for this file?
	if (_offsetsMap.count(this->filePath) != 0) {
		// we have information, now check if the timestamp is still okay
		if (GetLastModificationTime(this->filePath) != _offsetsMap[this->filePath].modificationTime) {
			// looks like the file was modified
			// delete the offset information, it will be recreated below
			_offsetsMap.erase(this->filePath);
		}
	}
	
	// now create the offsets map if needed (if it did not exist or was deleted above)
	if (_offsetsMap.count(this->filePath) == 0) {
		ImageLoaderHamamatsu::ImageOffsets imageOffsets;
		ImageLoaderHamamatsu_HeaderStructure header;
		
		file.seekg(0);
		
		int nImagesRead = 0;
		int nImagesInFile = -1;
		int thisImageXSize, thisImageYSize;
		uint64_t nextHeaderOffset = 0;
		do {
			file.seekg(nextHeaderOffset);
			
			file.read((char *)&(header.magic), sizeof(header.magic));
			file.read((char *)&(header.commentLength), sizeof(header.commentLength));
			file.read((char *)&(header.xSize), sizeof(header.xSize));
			file.read((char *)&(header.ySize), sizeof(header.ySize));
			file.read((char *)&(header.xBinning), sizeof(header.xBinning));
			file.read((char *)&(header.yBinning), sizeof(header.yBinning));
			file.read((char *)&(header.storageFormat), sizeof(header.storageFormat));
			file.read((char *)&(header.nImages), sizeof(header.nImages));
			file.read((char *)&(header.nChannels), sizeof(header.nChannels));
			file.read((char *)&(header.channel), sizeof(header.channel));
			file.read((char *)&(header.timeStamp), sizeof(header.timeStamp));
			file.read((char *)&(header.marker), sizeof(header.marker));
			file.read((char *)&(header.misc), sizeof(header.misc));
			
			if (nImagesInFile == -1)
				nImagesInFile = header.nImages;
			thisImageXSize = header.xSize;
			thisImageYSize = header.ySize;
			
			if (header.storageFormat != 2) {	// not UINT16
				std::string error;
				error = "The file at \"";
				error += this->filePath;
				error += "\" specifies that it doesn't use UINT16 for storage. This usually means that the manufacturer's software corrupted the file.";
				throw ERROR_READING_FILE_DATA(error);
			}
			
			// was there an error reading the file?
			if (file.fail() != 0) {
				std::string error;
				error = "Error parsing the header information in \"";
				error += this->filePath;
				error += "\" assuming the Hamamatsu format";
				throw ERROR_READING_FILE_DATA(error);
			}
			
			imageOffsets.offsets.push_back(nextHeaderOffset + header.commentLength + 64);
			nextHeaderOffset += header.commentLength + 64 + thisImageXSize * thisImageYSize * sizeof(uint16_t);
			
			nImagesRead += 1;
			
			#ifdef WITH_IGOR
			if ((nImagesRead % 20) == 0) {
				int abort = SpinProcess();
				if (abort)
					throw USER_ABORTED("user abort");
			}
			#endif
		} while (nImagesRead < nImagesInFile);
		
		// store the obtained offsets and information
		imageOffsets.xSize = thisImageXSize;
		imageOffsets.ySize = thisImageYSize;
		imageOffsets.modificationTime = GetLastModificationTime(this->filePath);
		
		_offsetsMap[this->filePath] = imageOffsets;
		
		file.seekg(0);
	}
	
	this->nImages = _offsetsMap[this->filePath].offsets.size();
	this->xSize = _offsetsMap[this->filePath].xSize;
	this->ySize = _offsetsMap[this->filePath].ySize;
	this->_offsets = _offsetsMap[this->filePath].offsets;
	this->storage_type = STORAGE_TYPE_UINT16;
	
	this->checkForReasonableValues();
}

ImagePtr ImageLoaderHamamatsu::readNextImage(size_t &indexOfImageThatWasRead) {
	uint64_t offset;
	uint64_t imageSize = xSize * ySize * 2; // assume a 16-bit format
	
	boost::scoped_array<char> single_image_buffer(new char[imageSize]);
	ImagePtr image (GetRecycledMatrix((int)xSize, (int)ySize), FreeRecycledMatrix);
	
	{
		boost::lock_guard<boost::mutex> locker(loadImagesMutex);
		if (this->nextImageToRead >= nImages)
			throw IMAGE_INDEX_BEYOND_N_IMAGES(std::string("Requested more images than there are in the file"));
		
		offset = _offsets.at(nextImageToRead);
		
		file.seekg(offset);
		file.read(single_image_buffer.get(), imageSize);
		if (file.fail() != 0) {
			std::string error;
			error = "Error reading image data from \"";
			error += this->filePath;
			error += "\" assuming the Hamamatsu format";
			throw ERROR_READING_FILE_DATA(error);
		}
		indexOfImageThatWasRead = this->nextImageToRead;
		this->nextImageToRead += 1;
	}
	
	CopyBufferToImage<uint16_t>(single_image_buffer.get(), image);
	
	return image;
}

ImageLoaderPDE::ImageLoaderPDE(std::string rhs) {
	this->filePath = rhs;
	
	file.open(this->filePath.c_str(), std::ios::binary | std::ios::in);
	if (file.fail() == 1) {
		std::string error ("Unable to open the file at ");
		error += this->filePath;
		throw CANNOT_OPEN_FILE(error);
	}
	
	parse_header_information();
}

ImageLoaderPDE::~ImageLoaderPDE() {
	if (file.is_open() == 1) {
		file.close();
	}
}

void ImageLoaderPDE::parse_header_information() {
	
	PDEFormatHeader header;
	
	file.seekg(0);
	file.read((char *)&(header.magic), sizeof(header.magic));
	file.read((char *)&(header.version), sizeof(header.version));
	file.read((char *)&(header.nImages), sizeof(header.nImages));
	file.read((char *)&(header.xSize), sizeof(header.xSize));
	file.read((char *)&(header.ySize), sizeof(header.ySize));
	file.read((char *)&(header.storageFormat), sizeof(header.storageFormat));
	
	if (file.fail() != 0) {
		std::string error;
		error = "Error parsing the header information in \"";
		error += this->filePath;
		error += "\" assuming the simple image format";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	if (header.magic != 27)
		throw std::runtime_error("The data is either not in the PDE format or of a different endianness");
	if (header.version != 1)
		throw std::runtime_error("Unsupported version of the PDE format.");
	
	this->nImages = header.nImages;
	this->xSize = header.xSize;
	this->ySize = header.ySize;
	this->storage_type = header.storageFormat;
	this->header_length = 24;
	
	this->checkForReasonableValues();
}

ImagePtr ImageLoaderPDE::readNextImage(size_t &indexOfImageThatWasRead) {
	size_t nPixels = this->xSize * this->ySize;
	uint64_t offset, imageSize;
	
	ImagePtr image(GetRecycledMatrix((int)this->xSize, (int)this->ySize), FreeRecycledMatrix);
	
	switch (this->storage_type) {
		case STORAGE_TYPE_INT8:
		case STORAGE_TYPE_UINT8:
			imageSize = nPixels;
			break;
		case STORAGE_TYPE_INT16:
		case STORAGE_TYPE_UINT16:
			imageSize = nPixels * 2;
			break;
		case STORAGE_TYPE_INT32:
		case STORAGE_TYPE_UINT32:
		case STORAGE_TYPE_FP32:
			imageSize = nPixels * 4;
			break;
		case STORAGE_TYPE_FP64:
			imageSize = nPixels * 8;
			break;
		default:
			throw std::runtime_error("The data file does not appear to contain a recognized storage type");
			break;
	}
	
	boost::scoped_array<char> buffer(new char[imageSize]);
	
	{
		boost::lock_guard<boost::mutex> locker(loadImagesMutex);
		if (this->nextImageToRead >= this->nImages)
			throw IMAGE_INDEX_BEYOND_N_IMAGES(std::string("Requested more images than there are in the file"));
		
		offset = this->header_length + this->nextImageToRead * imageSize;
		
		file.seekg(offset);
		file.read(buffer.get(), imageSize);
		if (file.fail() != 0) {
			std::string error;
			error = "Error reading image data from \"";
			error += this->filePath;
			error += "\" assuming the simple image format";
			throw ERROR_READING_FILE_DATA(error);
		}
		indexOfImageThatWasRead = this->nextImageToRead;
		this->nextImageToRead += 1;
	}
	
	
	switch (this->storage_type) {
        case STORAGE_TYPE_INT8:
            CopyBufferToImage<int8_t>(buffer.get(), image, 1);
            break;
            
        case STORAGE_TYPE_UINT8:
            CopyBufferToImage<uint8_t>(buffer.get(), image, 1);
            break;
            
        case STORAGE_TYPE_INT16:
            CopyBufferToImage<int16_t>(buffer.get(), image, 1);
            break;
            
        case STORAGE_TYPE_UINT16:
            CopyBufferToImage<uint16_t>(buffer.get(), image, 1);
            break;
            
		case STORAGE_TYPE_UINT32:
            CopyBufferToImage<uint32_t>(buffer.get(), image, 1);
            break;
            
		case STORAGE_TYPE_FP32:
            CopyBufferToImage<float>(buffer.get(), image, 1);
            break;
            
		case STORAGE_TYPE_FP64:
            CopyBufferToImage<double>(buffer.get(), image, 1);
            break;
            
		default:
			throw std::runtime_error("The data file does not appear to contain a recognized storage type");
			break;
	}
	
	return image;
}

ImageLoaderTIFF::ImageLoaderTIFF(std::string rhs) {
	TIFFSetWarningHandler(NULL);
	this->filePath = rhs;
	
	this->currentDirectoryIndex = 0;
	
	tiff_file = NULL;
	tiff_file = TIFFOpen(this->filePath.c_str(), "rm");
	if (tiff_file == NULL) {
		std::string error ("Unable to open the file at ");
		error += this->filePath;
		throw CANNOT_OPEN_FILE(error);
	}
	
	parse_header_information();
}


ImageLoaderTIFF::~ImageLoaderTIFF() {
	
	if (tiff_file != NULL) {
		TIFFClose(tiff_file);
	}
	
}


void ImageLoaderTIFF::parse_header_information() {
	int result;
	uint16_t result_uint16;
	uint32_t result_uint32;
	int isInt;
	size_t bitsPerPixel;
	
	// how many images are there in the file?
	// because there is a possibility of a file with different
	// subfile types (e.g. Zeiss LSM files), we need to go over all of them
	// and look at which subfiles are the ones we're interested in
	result = TIFFSetDirectory(tiff_file, 0);
	if (result != 1)
		throw std::runtime_error("Error reading from the file");
	
	for (size_t index = 0; ; ++index) {
		result = TIFFGetField(tiff_file, TIFFTAG_SUBFILETYPE, &result_uint32);
		
		if (((result_uint32 == FILETYPE_REDUCEDIMAGE) || (result_uint32 == FILETYPE_MASK)) && (result == 1)) {
			// this is one of the subtypes that we don't support
			// do nothing with it
		} else {
			// if we're here then the image is appropriate, store its index
			directoryIndices.push_back(index);
		}
		
		if (TIFFReadDirectory(tiff_file) != 1) {
			// there are no more directories in the file,
			break;
		}
	}
	
	
	this->nImages = directoryIndices.size();
	
	// obtain some properties from the first image in the file
	// we assume that these properties are the for all other valid subfiles
	result = TIFFSetDirectory(tiff_file, directoryIndices[0]);
	if (result != 1)
		throw std::runtime_error("Error reading from the file");
	
	// is the image in grayscale format?
	result = TIFFGetField(tiff_file, TIFFTAG_PHOTOMETRIC, &result_uint16);
	if (result != 1) {
		std::string error;
		error = "The image at\"";
		error += this->filePath;
		error += "\" is not a grayscale image";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	if ((result_uint16 != 0) && (result_uint16 != 1)) {	// not a grayscale image
		std::string error;
		error = "The image at\"";
		error += this->filePath;
		error += "\" is not a grayscale image";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	// is it a binary image?
	result = TIFFGetField(tiff_file, TIFFTAG_BITSPERSAMPLE, &result_uint16);
	if (result != 1) {
		std::string error;
		error = "The image at\"";
		error += this->filePath;
		error += "\" is not a grayscale image";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	if (result_uint16 < 4) {	// 4 is the minimum number of bits allowed for grayscale images in the tiff specification, so this is a bilevel image
		std::string error;
		error = "The image at\"";
		error += this->filePath;
		error += "\" is not a grayscale image";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	bitsPerPixel = (unsigned int)result_uint16;
	
	// is the data in unsigned integer or floating point format?
	result = TIFFGetField(tiff_file, TIFFTAG_SAMPLEFORMAT, &result_uint16);
	if (result != 1) {	// if the field does not exist then we assume that it is integer format
		result_uint16 = 1;
	}
	
	switch (result_uint16) {
		case 1:
			isInt = 1;
			break;
		case 3:
			isInt = 0;
			break;
		default:
			std::string error;
			error = "The SampleFormat of the image at\"";
			error += this->filePath;
			error += "\" is unknown";
			throw ERROR_READING_FILE_DATA(error);
			break;
	}
	
	if (isInt == 1) {
		switch (bitsPerPixel) {
			case 4:
				storage_type = STORAGE_TYPE_UINT4;
				break;
			case 8:
				storage_type = STORAGE_TYPE_UINT8;
				break;
			case 16:
				storage_type = STORAGE_TYPE_UINT16;
				break;
			case 32:
				storage_type = STORAGE_TYPE_UINT32;
				break;
			default:
				std::string error;
				error = "The SampleFormat of the image at\"";
				error += this->filePath;
				error += "\" is unknown";
				throw ERROR_READING_FILE_DATA(error);
				break;
		}
	} else {	// the image is not in integer format but is a floating point
		switch (bitsPerPixel) {
			case 32:
				storage_type = STORAGE_TYPE_FP32;
				break;
			case 64:
				storage_type = STORAGE_TYPE_FP64;
				break;
			default:
				std::string error;
				error = "The SampleFormat of the image at\"";
				error += this->filePath;
				error += "\" is unknown";
				throw ERROR_READING_FILE_DATA(error);
				break;
		}
	}
	
	// what is the x size?
	result = TIFFGetField(tiff_file, TIFFTAG_IMAGEWIDTH, &result_uint32);
	if (result != 1) {
		std::string error;
		error = "The image at\"";
		error += this->filePath;
		error += "\" does not specify a width";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	xSize = (size_t)(result_uint32);
	
	// what is the y size?
	result = TIFFGetField(tiff_file, TIFFTAG_IMAGELENGTH, &result_uint32);
	if (result != 1) {
		std::string error;
		error = "The image at\"";
		error += this->filePath;
		error += "\" does not specify a height";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	ySize = (size_t)(result_uint32);
	
	result = TIFFSetDirectory(tiff_file, directoryIndices[0]);
	if (result != 1) {
		std::string error;
		error = "Unable to set the directory for the image at\"";
		error += this->filePath;
		error += "\"";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	this->checkForReasonableValues();
}

ImagePtr ImageLoaderTIFF::readNextImage(size_t &indexOfImageThatWasRead) {
	int result;
	
	boost::lock_guard<boost::mutex> lock(loadImagesMutex);
	
	std::shared_ptr<void> single_scanline_buffer (_TIFFmalloc(TIFFScanlineSize(tiff_file)), _TIFFfree);
	if (single_scanline_buffer.get() == NULL) {
		throw std::bad_alloc();
	}
	
	if (this->nextImageToRead >= this->nImages)
		throw IMAGE_INDEX_BEYOND_N_IMAGES(std::string("Requested more images than there are in the file"));
	
	// advance the active TIFF directory to that needed to read the next image
	while (this->currentDirectoryIndex != directoryIndices.at(this->nextImageToRead)) {
		result = TIFFReadDirectory(this->tiff_file);
		this->currentDirectoryIndex += 1;
		if (result != 1) {
			std::string error;
			error = "Unable to set the directory for the image at\"";
			error += this->filePath;
			error += "\"";
			throw ERROR_READING_FILE_DATA(error);
		}
	}
	
	ImagePtr image(GetRecycledMatrix((int)xSize, (int)ySize), FreeRecycledMatrix);
	
	for (size_t j = 0; j < ySize; ++j) {
		result = TIFFReadScanline(tiff_file, single_scanline_buffer.get(), j, 0);	// sample is ignored
		if (result != 1) {
			std::string error;
			error = "Unable to read a scanline from the image at\"";
			error += this->filePath;
			error += "\"";
			throw ERROR_READING_FILE_DATA(error);
		}
		
		switch (storage_type) {	// handle the different possibilities (floating, integer) and variable sizes
			case STORAGE_TYPE_UINT4:
			{
				char *charPtr = (char *)single_scanline_buffer.get();
				for (size_t k = 0; k < xSize; k+=2) {
					(*image)(k, j) = (double)(0x0000000F & (*charPtr));
					(*image)(k + 1, j) = (double)((0x000000F0 & (*charPtr)) / 16);			  
					++charPtr;
				}
				break;
			}
			case STORAGE_TYPE_UINT8:
			{
				char *charPtr = (char *)single_scanline_buffer.get();
				for (size_t k = 0; k < xSize; ++k) {
					(*image)(k, j) = (double)(*charPtr);
					++charPtr;
				}
				break;
			}
			case STORAGE_TYPE_UINT16:
			{
				uint16_t *uint16tPtr = (uint16_t *)single_scanline_buffer.get();
				for (size_t k = 0; k < xSize; ++k) {
					(*image)(k, j) = (double)(*uint16tPtr);
					++uint16tPtr;
				}
				break;
			}
			case STORAGE_TYPE_UINT32:
			{
				uint32_t *uint32tPtr = (uint32_t *)single_scanline_buffer.get();
				for (size_t k = 0; k < xSize; ++k) {
					(*image)(k, j) = (double)(*uint32tPtr);
					++uint32tPtr;
				}
				break;
			}
			case STORAGE_TYPE_FP32:
			{
				float *floatPtr = (float *)single_scanline_buffer.get();
				for (size_t k = 0; k < xSize; ++k) {
					(*image)(k, j) = (double)(*floatPtr);
					++floatPtr;
				}
				break;
			}
			case STORAGE_TYPE_FP64:
			{
				double *doublePtr = (double *)single_scanline_buffer.get();
				for (size_t k = 0; k < xSize; ++k) {
					(*image)(k, j) = (double)(*doublePtr);
					++doublePtr;
				}
				break;
			}
			default:
				std::string error;
				error = "Invalid floating point data size for the image at\"";
				error += this->filePath;
				error += "\"";
				throw ERROR_READING_FILE_DATA(error);
				break;
		}
	}
	
	indexOfImageThatWasRead = this->nextImageToRead;
	this->nextImageToRead += 1;
	
	return image;
}

void ImageLoaderTIFF::spoolTo(size_t index) {
	int result;
	
	// don't do any work unless it is required
	if (this->nextImageToRead != index) {
		result = TIFFSetDirectory(this->tiff_file, this->directoryIndices.at(index));
		if (result != 1) {
			std::string error;
			error = "Unable to set the directory for the image at\"";
			error += this->filePath;
			error += "\"";
			throw ERROR_READING_FILE_DATA(error);
		}
		
		this->nextImageToRead = index;
		this->currentDirectoryIndex = this->directoryIndices.at(index);
	}
}

ImageLoaderMultiFileTIFF::ImageLoaderMultiFileTIFF(std::string filePath) {
	// assume that filePath consists of something such as the following
	// /path/to/folder/XXXX123.YYY where X is character (but the last one is not a number)
	// and Y is the extension. "/path/to/folder/XXXX" is the base file path, YYY is the
	// extension, and the number of digits in between is assumed to be constant for every
	// file in the sequence
	
	size_t stringLength = filePath.length();
	size_t extensionStartIndex = filePath.rfind('.');
	size_t lengthOfExtension;
	if (extensionStartIndex == std::string::npos) {
		// no extension. Maybe this is OK but it probably is an error.
		// pretend that everything is fine.
		extension.clear();
	} else {
		extension = filePath.substr(extensionStartIndex);
	}
	lengthOfExtension = extension.length();
	
	// now determine how many digits there are that indicate the position of each file in the
	// frame sequence. Assume that all digits that precede the extension are part of this index.
	int nIndexDigits = 0;
	for (int i = extensionStartIndex - 1; i >= 0; --i) {
		if (isdigit(filePath[i])) {
			nIndexDigits += 1;
		} else {
			break;
		}
	}
	if (nIndexDigits == 0) {
		std::string error("Unable to determine the numbering format for the multi-tiff file at \"");
		error += filePath;
		error += "\"";
		throw ERROR_READING_FILE_DATA(error);
	}
	if (nIndexDigits > 9) {	// arbitrary limit
		std::string error("The numeric index of the images at \"");
		error += filePath;
		error += "\" appears to contain too many digits to be realistic";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	nDigitsInNumber = nIndexDigits;
	
	// everything coming before the index digits is the base file path
	this->baseFilePath = filePath.substr(0, stringLength - lengthOfExtension - nDigitsInNumber);
	std::string indexOfThisImageStr = filePath.substr(stringLength - lengthOfExtension - nDigitsInNumber, nDigitsInNumber);
	
	// try to load the provided image to get the metadata
	ImageLoaderTIFF imageLoaderTIFF(filePath);
	this->xSize = imageLoaderTIFF.getXSize();
	this->ySize = imageLoaderTIFF.getYSize();
	this->storage_type = imageLoaderTIFF.getStorageType();
	
	// now figure out how many images there are in the dataset
	// get the index of the current image
	int indexOfThisImage;
	int nItemsMatched = sscanf(indexOfThisImageStr.c_str(), "%d", &indexOfThisImage);
	if (nItemsMatched != 1) {
		std::string error("Error parsing the numeric index for the image at \"");
		error += filePath;
		error += "\"";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	std::pair<int, int> firstAndLastIndices = findFirstAndLastValidImageIndices(indexOfThisImage);
	this->firstImageIndex = firstAndLastIndices.first;
	this->nImages = firstAndLastIndices.second - firstAndLastIndices.first + 1;
	
	this->checkForReasonableValues();
}

ImagePtr ImageLoaderMultiFileTIFF::readNextImage(size_t &indexOfImageThatWasRead) {
	boost::lock_guard<boost::mutex> lock(loadImagesMutex);
	
	if (this->nextImageToRead >= this->nImages)
		throw IMAGE_INDEX_BEYOND_N_IMAGES(std::string("Requested more images than there are in the file"));
	
	std::string filePath = getFilePathForImageAtIndex(this->nextImageToRead);
	ImageLoaderTIFF imageLoaderTIFF(filePath);
	size_t dummy;
	ImagePtr image = imageLoaderTIFF.readNextImage(dummy);
	
	indexOfImageThatWasRead = this->nextImageToRead;
	this->nextImageToRead += 1;
	
	return image;
}

std::string ImageLoaderMultiFileTIFF::getFilePathForImageAtIndex(int index) {
	int imageIndexOnDisk = this->firstImageIndex + index;
	
	boost::scoped_array<char> imageIndexStr(new char[this->nDigitsInNumber + 1]);
	char formatString[10];
	sprintf(formatString, "%%0%dd", nDigitsInNumber);
	sprintf(imageIndexStr.get(), formatString, imageIndexOnDisk);
	
	std::string filePath = this->baseFilePath;
	filePath += imageIndexStr.get();
	filePath += this->extension;
	return filePath;
}

bool ImageLoaderMultiFileTIFF::imageFileAtIndexExists(int index) {
	std::string thisFilePath = getFilePathForImageAtIndex(index);
	
	bool exists = true;
	try {
		ImageLoaderTIFF imageLoaderTIFF(thisFilePath);
	}
	catch (...) {
		exists = false;
	}
	
	return exists;
}

std::pair<int, int> ImageLoaderMultiFileTIFF::findFirstAndLastValidImageIndices(int knownValidImageIndex) {
	std::pair<int, int> firstAndLastIndices;
	
	assert(knownValidImageIndex >= 0);
	
	int trialIndex;
	int lower, upper, mid;
	
	// find the first image index that is valid
	if (imageFileAtIndexExists(0)) {
		// check if the first image is simply at index zero
		firstAndLastIndices.first = 0;
	} else {
		// find an invalid index less than knownValidImageIndex
		trialIndex = knownValidImageIndex;
		for (int delta = -1; ; delta *= 2) {
			trialIndex = trialIndex + delta;
			if (trialIndex < 0) {
				trialIndex = -1;	// the image at -1 certainly doesn't exist
				break;
			}
			if (!imageFileAtIndexExists(trialIndex))
				break;
		}
		
		// now we have a known bad and a known good index. Find the first index that is valid.
		lower = trialIndex;
		upper = knownValidImageIndex;
		for (;;) {
			if (upper - lower == 1) {
				firstAndLastIndices.first = upper;
				break;
			}
			
			mid = (lower + upper) / 2;
			if (imageFileAtIndexExists(mid)) {
				upper = mid;
			} else {
				lower = mid;
			}
		}
	}
	
	// now find the last index that is valid. First find an index that is invalid.
	trialIndex = knownValidImageIndex + 1;
	for (int delta = 1; ; delta *= 2) {
		trialIndex = trialIndex + delta;
		if (!imageFileAtIndexExists(trialIndex))
			break;
	}
	
	lower = knownValidImageIndex; upper = trialIndex;
	for (;;) {
		if (upper - lower == 1) {
			firstAndLastIndices.second = lower;
			break;
		}
		
		mid = (lower + upper) / 2;
		if (imageFileAtIndexExists(mid)) {
			lower = mid;
		} else {
			upper = mid;
		}
	}
	
	return firstAndLastIndices;
}

#ifdef WITH_IGOR
ImageLoaderIgor::ImageLoaderIgor(std::string waveName) {
	int waveType;
	
	// try to get images from an Igor wave
	
	this->igor_data_wave = FetchWaveUsingFullPath(waveName);
	
	size_t DimensionSizes[MAX_DIMENSIONS + 1];
	int numDimensions;
	int result;
	
	result = MDGetWaveDimensions(igor_data_wave, &numDimensions, (CountInt *)DimensionSizes);
	if (result != 0) {
		throw result;
	}
	if ((numDimensions != 2) && (numDimensions != 3)) {
		throw INCOMPATIBLE_DIMENSIONING;
	}
	
	xSize = (size_t)DimensionSizes[0];
	ySize = (size_t)DimensionSizes[1];
	nImages = (size_t)DimensionSizes[2];
	
	// special case: if the wave contains only a single image
    // then it is usually two-dimensional, that is, DimensionSizes[2] == 0
	// in that case nImages is still one
	if (DimensionSizes[2] == 0) {
		nImages = 1;
	}
	
	waveType = WaveType(igor_data_wave);
    // do not handle complex waves
    if (waveType & NT_CMPLX)
        throw NO_COMPLEX_WAVE;
    
	switch (waveType) {
		case NT_I8:
			storage_type = STORAGE_TYPE_INT8;
			break;
		case NT_I16:
			storage_type = STORAGE_TYPE_INT16;
			break;
		case NT_I32:
			storage_type = STORAGE_TYPE_INT32;
			break;
		case (NT_I8 | NT_UNSIGNED):
			storage_type = STORAGE_TYPE_UINT8;
			break;
		case (NT_I16 | NT_UNSIGNED):
			storage_type = STORAGE_TYPE_UINT16;
			break;
		case (NT_I32 | NT_UNSIGNED):
			storage_type = STORAGE_TYPE_UINT32;
			break;
		case NT_FP32:
			storage_type = STORAGE_TYPE_FP32;
			break;
		case NT_FP64:
			storage_type = STORAGE_TYPE_FP64;
			break;
		default:
			storage_type = STORAGE_TYPE_FP64;
	}
	
	this->checkForReasonableValues();
}

ImagePtr ImageLoaderIgor::readNextImage(size_t &index) {
	int err;
	int waveType = WaveType(this->igor_data_wave);
	size_t waveDataOffset;
	char* startOfWaveData;
	size_t nPixels = this->xSize * this->ySize;
	
	// allocate a new image
	ImagePtr image (GetRecycledMatrix((int)xSize, (int)ySize), FreeRecycledMatrix);
	
	// get a pointer to the data in the wave
	err = MDAccessNumericWaveData(this->igor_data_wave, kMDWaveAccessMode0, (BCInt*)&waveDataOffset);
	if (err != 0) {
		throw err;
	}
	startOfWaveData = ((char*)(*this->igor_data_wave) + waveDataOffset);
	
	{
		boost::lock_guard<boost::mutex> lock(this->loadImagesMutex);
		if (this->nextImageToRead >= nImages)
			throw IMAGE_INDEX_BEYOND_N_IMAGES(std::string("Requested more images than there are in the file"));
		
		index = this->nextImageToRead;
		this->nextImageToRead += 1;
	}
	
	// get a pointer to the image
	double* imagePtr = image->data();
	
	switch (waveType) {
		case NT_FP32:
		{
			float* floatPtr = (float*)startOfWaveData;
			floatPtr += nPixels * index;
			for (size_t i = 0; i < nPixels; ++i) {
				imagePtr[i] = floatPtr[i];
			}
			break;
		}
		case NT_FP64:
		{
			double* doublePtr = (double*)startOfWaveData;
			doublePtr += nPixels * index;
			memcpy(imagePtr, doublePtr, nPixels * sizeof(double));
			break;
		}
		case NT_I8:
		{
			int8_t* int8Ptr = (int8_t*)startOfWaveData;
			int8Ptr += nPixels * index;
			for (size_t i = 0; i < nPixels; ++i) {
				imagePtr[i] = int8Ptr[i];
			}
			break;
		}
		case NT_I16:
		{
			int16_t* int16Ptr = (int16_t*)startOfWaveData;
			int16Ptr += nPixels * index;
			for (size_t i = 0; i < nPixels; ++i) {
				imagePtr[i] = int16Ptr[i];
			}
			break;
		}
		case NT_I32:
		{
			int32_t* int32Ptr = (int32_t*)startOfWaveData;
			int32Ptr += nPixels * index;
			for (size_t i = 0; i < nPixels; ++i) {
				imagePtr[i] = int32Ptr[i];
			}
			break;
		}
		case NT_I8 | NT_UNSIGNED:
		{
			uint8_t* uint8Ptr = (uint8_t*)startOfWaveData;
			uint8Ptr += nPixels * index;
			for (size_t i = 0; i < nPixels; ++i) {
				imagePtr[i] = uint8Ptr[i];
			}
			break;
		}
		case NT_I16 | NT_UNSIGNED:
		{
			uint16_t* uint16Ptr = (uint16_t*)startOfWaveData;
			uint16Ptr += nPixels * index;
			for (size_t i = 0; i < nPixels; ++i) {
				imagePtr[i] = uint16Ptr[i];
			}
			break;
		}
		case NT_I32 | NT_UNSIGNED:
		{
			uint32_t* uint32Ptr = (uint32_t*)startOfWaveData;
			uint32Ptr += nPixels * index;
			for (size_t i = 0; i < nPixels; ++i) {
				imagePtr[i] = uint32Ptr[i];
			}
			break;
		}
		default:
			throw std::runtime_error("Unknown or unsupported wavetype");
			break;
	}
	
	return image;
}
#endif // WITH_IGOR

#ifdef WITH_MATLAB
ImageLoaderMatlab::ImageLoaderMatlab(mxArray* matlabArray) {
	_matlabArray = matlabArray;
	
	mwSize nDims = mxGetNumberOfDimensions(matlabArray);
	if ((nDims < 2) || (nDims > 3))
		throw std::runtime_error("Expected 2D or 3D matrix");
	
	const mwSize* dimensionSizes = mxGetDimensions(matlabArray);
	xSize = static_cast<size_t>(dimensionSizes[0]);
	ySize = static_cast<size_t>(dimensionSizes[1]);
	if (nDims == 3)
		nImages = static_cast<size_t>(dimensionSizes[2]);
	else
		nImages = 1;
	
	this->storage_type = static_cast<int>(mxGetClassID(matlabArray));
	this->checkForReasonableValues();
}

ImagePtr ImageLoaderMatlab::readNextImage(size_t &index) {
	size_t nPixels = this->xSize * this->ySize;
	
	{
		boost::lock_guard<boost::mutex> lock(this->loadImagesMutex);
		if (this->nextImageToRead >= nImages)
			throw IMAGE_INDEX_BEYOND_N_IMAGES(std::string("Requested more images than there are in the file"));
		
		index = this->nextImageToRead;
		this->nextImageToRead += 1;
	}
	
	// allocate a new image
	ImagePtr image (GetRecycledMatrix((int)xSize, (int)ySize), FreeRecycledMatrix);
	
	char* dataPtr = reinterpret_cast<char*>(mxGetData(_matlabArray));
	
	switch (mxGetClassID(_matlabArray)) {
		case mxINT8_CLASS:
		{
			size_t offset = index * nPixels * sizeof(int8_t);
			CopyBufferToImage<int8_t>(dataPtr + offset, image);
			break;
		}
		case mxUINT8_CLASS:
		{
			size_t offset = index * nPixels * sizeof(uint8_t);
			CopyBufferToImage<uint8_t>(dataPtr + offset, image);
			break;
		}
		case mxINT16_CLASS:
		{
			size_t offset = index * nPixels * sizeof(int16_t);
			CopyBufferToImage<int16_t>(dataPtr + offset, image);
			break;
		}
		case mxUINT16_CLASS:
		{
			size_t offset = index * nPixels * sizeof(uint16_t);
			CopyBufferToImage<uint16_t>(dataPtr + offset, image);
			break;
		}
		case mxINT32_CLASS:
		{
			size_t offset = index * nPixels * sizeof(int32_t);
			CopyBufferToImage<int32_t>(dataPtr + offset, image);
			break;
		}
		case mxUINT32_CLASS:
		{
			size_t offset = index * nPixels * sizeof(uint32_t);
			CopyBufferToImage<uint32_t>(dataPtr + offset, image);
			break;
		}
		case mxINT64_CLASS:
		{
			size_t offset = index * nPixels * sizeof(int64_t);
			CopyBufferToImage<int64_t>(dataPtr + offset, image);
			break;
		}
		case mxUINT64_CLASS:
		{
			size_t offset = index * nPixels * sizeof(uint64_t);
			CopyBufferToImage<uint64_t>(dataPtr + offset, image);
			break;
		}
		case mxSINGLE_CLASS:
		{
			size_t offset = index * nPixels * sizeof(float);
			CopyBufferToImage<float>(dataPtr + offset, image);
			break;
		}
		case mxDOUBLE_CLASS:
		{
			size_t offset = index * nPixels * sizeof(double);
			CopyBufferToImage<double>(dataPtr + offset, image);
			break;
		}
		default:
			throw std::runtime_error("Unknown or unsupported matrix class ID");
			break;
	}
	
	return image;
}
#endif // WITH_MATLAB

ImageOutputWriter::ImageOutputWriter() {
	outputFilePath.assign("");
	nImagesWritten = 0;
}

PDEImageOutputWriter::PDEImageOutputWriter(const std::string &rhs,int overwrite, uint32_t storageType_rhs) {
	// if overwrite is non-zero then we overwrite any file that exists at the output path
	// if it is set to zero then we throw an error and abort instead of overwriting
	outputFilePath = rhs;
	
	if (overwrite == 0) {
		std::ifstream input_test;
		input_test.open(outputFilePath.c_str(), std::ios::in | std::ios::binary);
		if (input_test.good() == 1) {
			std::string error("The output file at ");
			error += this->outputFilePath;
			error += " already exists.";
			throw OUTPUT_FILE_ALREADY_EXISTS(error);	// escape without overwriting
		}
	}
	
	file.open(outputFilePath.c_str(), std::ios::binary | std::ios::out | std::ios::trunc);	// DANGER: OVERWRITING THE FILE
	if (file.fail() != 0) {
		std::string error("Cannot create an output file at ");
		error += this->outputFilePath;
		throw CANNOT_OPEN_OUTPUT_FILE(error);
	}
	
	this->xSize = 0;
	this->ySize = 0;
	this->nImagesWritten = 0;
	this->storageType = storageType_rhs;
	
	// write the header out in advance
	// before closing the file, after all images have been written,
	// we will come back to this and overwrite it with the correct
	// values
	this->WriteHeader();
}


PDEImageOutputWriter::~PDEImageOutputWriter() {
	if (file.is_open() != 0) {
		WriteHeader();
		file.close();
	}
}

void PDEImageOutputWriter::WriteHeader() {
	PDEFormatHeader header;
	
	header.magic = 27;
	header.version = 1;
	header.nImages = this->nImagesWritten;
	header.xSize = this->xSize;
	header.ySize = this->ySize;
	header.storageFormat = this->storageType;
	
	if (file.is_open()) {
		file.seekp(0, std::ios_base::beg);
		file.write((char *)&(header.magic), sizeof(header.magic));
		file.write((char *)&(header.version), sizeof(header.version));
		file.write((char *)&(header.nImages), sizeof(header.nImages));
		file.write((char *)&(header.xSize), sizeof(header.xSize));
		file.write((char *)&(header.ySize), sizeof(header.ySize));
		file.write((char *)&(header.storageFormat), sizeof(header.storageFormat));
		if (file.fail())
			throw std::runtime_error("Error trying to write to the PDE file");
	}
}

void PDEImageOutputWriter::write_image(ImagePtr imageToWrite) {
	
	// determine the size of the frames
	size_t currentXSize = imageToWrite->rows();
	size_t currentYSize = imageToWrite->cols();
	size_t nPixels = currentXSize * currentYSize;
	
	if (this->nImagesWritten == 0) {
		this->xSize = currentXSize;
		this->ySize = currentYSize;
	} else {
		if ((currentXSize != this->xSize) || (currentYSize != this->ySize)) {
			throw std::runtime_error("Tried to write an image with different dimensions to an SPE file");
		}
	}
	
	int bytesPerPixel;
	switch (this->storageType) {
		case STORAGE_TYPE_INT8:
		case STORAGE_TYPE_UINT8:
			bytesPerPixel = 1;
			break;
		case STORAGE_TYPE_INT16:
		case STORAGE_TYPE_UINT16:
			bytesPerPixel = 2;
			break;
		case STORAGE_TYPE_INT32:
		case STORAGE_TYPE_UINT32:
		case STORAGE_TYPE_FP32:
			bytesPerPixel = 4;
			break;
		case STORAGE_TYPE_FP64:
			bytesPerPixel = 8;
			break;
		default:
			throw std::runtime_error("Unsupport file type requested in the PDE output format");
	}
	
	int nBytesToWrite = nPixels * bytesPerPixel;
	boost::scoped_array<char> buffer(new char[nBytesToWrite]);
	
	switch (this->storageType) {
		case STORAGE_TYPE_INT8:
			CopyImageToBuffer<int8_t>(imageToWrite, buffer.get(), 1);
			break;
		case STORAGE_TYPE_UINT8:
			CopyImageToBuffer<uint8_t>(imageToWrite, buffer.get(), 1);
			break;
		case STORAGE_TYPE_INT16:
			CopyImageToBuffer<int16_t>(imageToWrite, buffer.get(), 1);
			break;
		case STORAGE_TYPE_UINT16:
			CopyImageToBuffer<uint16_t>(imageToWrite, buffer.get(), 1);
			break;
		case STORAGE_TYPE_INT32:
			CopyImageToBuffer<int32_t>(imageToWrite, buffer.get(), 1);
			break;
		case STORAGE_TYPE_UINT32:
			CopyImageToBuffer<uint32_t>(imageToWrite, buffer.get(), 1);
			break;
		case STORAGE_TYPE_FP32:
			CopyImageToBuffer<float>(imageToWrite, buffer.get(), 1);
			break;
		case STORAGE_TYPE_FP64:
			CopyImageToBuffer<double>(imageToWrite, buffer.get(), 1);
			break;
	}
	
	this->file.write(buffer.get(), nBytesToWrite);
	++this->nImagesWritten;
}

TIFFImageOutputWriter::TIFFImageOutputWriter(const std::string &rhs,int overwrite, int compression_rhs, int storageType_rhs) {
	// if overwrite is non-zero then we overwrite any file that exists at the output path
	// if it is set to zero then we throw an error and abort instead of overwriting
	
	TIFFSetWarningHandler(NULL);
	this->outputFilePath = rhs;
	this->compression = compression_rhs;
	this->storageType = storageType_rhs;
	
	if (overwrite == 0) {
		std::ifstream input_test;
		input_test.open(outputFilePath.c_str(), std::ios::in | std::ios::binary);
		if (input_test.good() == 1) {
			std::string error("The output file at ");
			error += this->outputFilePath;
			error += " already exists.";
			throw OUTPUT_FILE_ALREADY_EXISTS(error);	// escape without overwriting
		}
		input_test.close();
	}
	
	tiff_file = TIFFOpen(outputFilePath.c_str(), "w");
	if (tiff_file == NULL) {
		std::string error("Cannot create an output file at ");
		error += this->outputFilePath;
		throw CANNOT_OPEN_OUTPUT_FILE(error);
	}
	
	nImagesWritten = 0;
}


TIFFImageOutputWriter::~TIFFImageOutputWriter() {
	if (tiff_file != NULL) {
		TIFFClose(tiff_file);
		tiff_file = NULL;
	}
}



void TIFFImageOutputWriter::write_image(ImagePtr imageToWrite) {
	
	size_t xSize = imageToWrite->rows();
	size_t ySize = imageToWrite->cols();
	
	int result, sampleFormat, bitsPerSample;
	uint16_t current_uint16;
	uint32_t current_uint32;
	
	// determine the output storage type
	switch (this->storageType) {
		case STORAGE_TYPE_INT4:
		case STORAGE_TYPE_UINT4:
		case STORAGE_TYPE_INT8:
			bitsPerSample = 8;
			sampleFormat = SAMPLEFORMAT_INT;
			break;
		case STORAGE_TYPE_UINT8:
			bitsPerSample = 8;
			sampleFormat = SAMPLEFORMAT_UINT;
			break;
		case STORAGE_TYPE_INT16:
			bitsPerSample = 16;
			sampleFormat = SAMPLEFORMAT_INT;
			break;
		case STORAGE_TYPE_UINT16:
			bitsPerSample = 16;
			sampleFormat = SAMPLEFORMAT_UINT;
			break;
		case STORAGE_TYPE_INT32:
			bitsPerSample = 32;
			sampleFormat = SAMPLEFORMAT_INT;
			break;
		case STORAGE_TYPE_UINT32:
			bitsPerSample = 32;
			sampleFormat = SAMPLEFORMAT_UINT;
			break;
		case STORAGE_TYPE_INT64:
			bitsPerSample = 64;
			sampleFormat = SAMPLEFORMAT_INT;
			break;
		case STORAGE_TYPE_UINT64:
			bitsPerSample = 64;
			sampleFormat = SAMPLEFORMAT_UINT;
			break;
		case STORAGE_TYPE_FP32:
			bitsPerSample = 32;
			sampleFormat = SAMPLEFORMAT_IEEEFP;
			break;
		case STORAGE_TYPE_FP64:
			bitsPerSample = 64;
			sampleFormat = SAMPLEFORMAT_IEEEFP;
			break;
		default:
			throw std::runtime_error("Unknown storage type requested for TIFF output");
			break;
	}
	
	// make a scoped_array that will act as a single scanline buffer
	// make it a buffer of chars equal to the total number of bytes required
	boost::scoped_array<char> scanLine(new char[xSize * bitsPerSample / 8]);
	
	// make sure that all the image tags have the correct values
	result = TIFFSetField(tiff_file, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
	if (result != 1) {
		std::string error;
		error = "Unable to set the photometric type for the image at\"";
		error += outputFilePath;
		error += "\"";
		throw ERROR_WRITING_FILE_DATA(error);
	}
	
	current_uint32 = xSize;
	result = TIFFSetField(tiff_file, TIFFTAG_IMAGEWIDTH, current_uint32);
	if (result != 1) {
		std::string error;
		error = "Unable to set the image width for the image at\"";
		error += outputFilePath;
		error += "\"";
		throw ERROR_WRITING_FILE_DATA(error);
	}
	
	current_uint32 = ySize;
	result = TIFFSetField(tiff_file, TIFFTAG_IMAGELENGTH, current_uint32);
	if (result != 1) {
		std::string error;
		error = "Unable to set the image height for the image at\"";
		error += outputFilePath;
		error += "\"";
		throw ERROR_WRITING_FILE_DATA(error);
	}
	
	result = TIFFSetField(tiff_file, TIFFTAG_SAMPLEFORMAT, sampleFormat);
	if (result != 1) {
		std::string error;
		error = "Unable to set the SampleFormat for the image at\"";
		error += outputFilePath;
		error += "\"";
		throw ERROR_WRITING_FILE_DATA(error);
	}
	
	result = TIFFSetField(tiff_file, TIFFTAG_BITSPERSAMPLE, bitsPerSample);
	if (result != 1) {
		std::string error;
		error = "Unable to set the BitsPerSample for the image at\"";
		error += outputFilePath;
		error += "\"";
		throw ERROR_WRITING_FILE_DATA(error);
	}
	
	result = TIFFSetField(tiff_file, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
	if (result != 1) {
		std::string error;
		error = "Unable to set the SubFileType for the image at\"";
		error += outputFilePath;
		error += "\"";
		throw ERROR_WRITING_FILE_DATA(error);
	}
	
	current_uint16 = (uint16_t)nImagesWritten;
	result = TIFFSetField(tiff_file, TIFFTAG_PAGENUMBER, current_uint16);
	if (result != 1) {
		std::string error;
		error = "Unable to set the PageNumber for the image at\"";
		error += outputFilePath;
		error += "\"";
		throw ERROR_WRITING_FILE_DATA(error);
	}
	
	result = TIFFSetField(tiff_file, TIFFTAG_COMPRESSION, compression);
	if (result != 1) {
		std::string error;
		error = "Unable to set the compression method for the image at\"";
		error += outputFilePath;
		error += "\"";
		throw ERROR_WRITING_FILE_DATA(error);
	}
	
	size_t offset = 0;
	for (size_t j = 0; j < ySize; j++) {
		offset = 0;
		
		switch (this->storageType) {
			case STORAGE_TYPE_INT4:
			case STORAGE_TYPE_UINT4:
			case STORAGE_TYPE_INT8:
			{
				int8_t* charPtr = (int8_t *)scanLine.get();
				for (size_t i = 0; i < xSize; i++) {
					charPtr[i] = (int8_t)(*imageToWrite)(i, j);
				}
				break;
			}
			case STORAGE_TYPE_UINT8:
			{
				uint8_t* uCharPtr = (uint8_t *)scanLine.get();
				for (size_t i = 0; i < xSize; i++) {
					uCharPtr[i] = (uint8_t)(*imageToWrite)(i, j);
				}
				break;
			}
			case STORAGE_TYPE_INT16:
			{
				int16_t* int16Ptr = (int16_t *)scanLine.get();
				for (size_t i = 0; i < xSize; i++) {
					int16Ptr[i] = (int16_t)(*imageToWrite)(i, j);
				}
				break;
			}
			case STORAGE_TYPE_UINT16:
			{
				uint16_t* uInt16Ptr = (uint16_t *)scanLine.get();
				for (size_t i = 0; i < xSize; i++) {
					uInt16Ptr[i] = (uint16_t)(*imageToWrite)(i, j);
				}
				break;
			}
			case STORAGE_TYPE_INT32:
			{
				int32_t* int32Ptr = (int32_t *)scanLine.get();
				for (size_t i = 0; i < xSize; i++) {
					int32Ptr[i] = (int32_t)(*imageToWrite)(i, j);
				}
				break;
			}
			case STORAGE_TYPE_UINT32:
			{
				uint32_t* uInt32Ptr = (uint32_t *)scanLine.get();
				for (size_t i = 0; i < xSize; i++) {
					uInt32Ptr[i] = (uint32_t)(*imageToWrite)(i, j);
				}
				break;
			};
			case STORAGE_TYPE_INT64:
			{
				int64_t* int64Ptr = (int64_t *)scanLine.get();
				for (size_t i = 0; i < xSize; i++) {
					int64Ptr[i] = (int64_t)(*imageToWrite)(i, j);
				}
				break;
			}
			case STORAGE_TYPE_UINT64:
			{
				uint64_t* uInt64Ptr = (uint64_t *)scanLine.get();
				for (size_t i = 0; i < xSize; i++) {
					uInt64Ptr[i] = (uint64_t)(*imageToWrite)(i, j);
				}
				break;
			}
			case STORAGE_TYPE_FP32:
			{
				float* floatPtr = (float *)scanLine.get();
				for (size_t i = 0; i < xSize; i++) {
					floatPtr[i] = (float)(*imageToWrite)(i, j);
				}
				break;
			}
			case STORAGE_TYPE_FP64:
			{
				double* doublePtr = (double *)scanLine.get();
				for (size_t i = 0; i < xSize; i++) {
					doublePtr[i] = (double)(*imageToWrite)(i, j);
				}
				break;
			}
		}
		
		result = TIFFWriteScanline(tiff_file, (char *)scanLine.get(), j);
		if (result != 1) {
			std::string error;
			error = "There was an error writing a scanline for the image at\"";
			error += outputFilePath;
			error += "\"";
			throw ERROR_WRITING_FILE_DATA(error);
		}
		
	}
	
	result = TIFFWriteDirectory(tiff_file);
	if (result != 1) {
		std::string error;
		error = "Unable to write a directory for the image at\"";
		error += outputFilePath;
		error += "\"";
		throw ERROR_WRITING_FILE_DATA(error);
	}
	
	++this->nImagesWritten;
}

MultiFileTIFFImageOutputWriter::MultiFileTIFFImageOutputWriter(const std::string &baseOutputFilePath_rhs, int overwrite_rhs, int compression_rhs, int storageType_rhs) {
	baseOutputFilePath = baseOutputFilePath_rhs;
	overwrite = overwrite_rhs;
	compression = compression_rhs;
	storageType = storageType_rhs;
}

void MultiFileTIFFImageOutputWriter::write_image(ImagePtr imageToWrite) {
	char imageIndexStr[50];
	sprintf(imageIndexStr, "%06d", static_cast<int>(this->nImagesWritten));
	
	std::string thisFileName = this->baseOutputFilePath + std::string(imageIndexStr) + std::string(".tif");
	
	TIFFImageOutputWriter singleFileOutputWriter(thisFileName, this->overwrite, this->compression, this->storageType);
	singleFileOutputWriter.write_image(imageToWrite);
	
	++this->nImagesWritten;
}

#ifdef WITH_IGOR
IgorImageOutputWriter::IgorImageOutputWriter(std::string waveName_rhs, size_t nImages_rhs, int overwrite_rhs, int storageType_rhs) {
	
	this->overwrite = overwrite_rhs;
	this->outputWave = NULL;
	this->fullPathToWave = waveName_rhs;
	this->nImagesTotal = nImages_rhs;
	this->nImagesWritten = 0;
    
    // not all storage types are supported in Igor
    // silently replace these with supported types
    switch (storageType_rhs) {
        case STORAGE_TYPE_INT4:
        case STORAGE_TYPE_UINT4:
            this->storageType = STORAGE_TYPE_INT8;
            break;
        case STORAGE_TYPE_INT64:
            this->storageType = STORAGE_TYPE_INT32;
            break;
        case STORAGE_TYPE_UINT64:
            this->storageType = STORAGE_TYPE_UINT32;
            break;
        default:
            this->storageType = storageType_rhs;
    }
	
	// check if the output wave already exists
	if (this->overwrite != 1) {
		waveHndl wave = FetchWaveUsingFullPath(waveName_rhs);
		if (wave != NULL) {
			std::string error ("the output wave ");
			error += waveName_rhs;
			error += " already exists (and the /O flag was not specified)";
		}
	}
}

IgorImageOutputWriter::IgorImageOutputWriter(DataFolderAndName outputDataFolderAndName_rhs, size_t nImages_rhs, int overwrite_rhs, int storageType_rhs) {
	
	this->overwrite = overwrite_rhs;
	this->outputWave = NULL;
	this->fullPathToWave.assign("");
	this->nImagesTotal = nImages_rhs;
	this->nImagesWritten = 0;
	this->storageType = storageType_rhs;
	this->waveDataFolderAndName = outputDataFolderAndName_rhs;
	
	// check if the output wave already exists
	if (this->overwrite != 1) {
		waveHndl wave = FetchWaveFromDataFolder(this->waveDataFolderAndName.dfH, this->waveDataFolderAndName.name);
		if (wave != NULL) {
			throw std::runtime_error("the requested output wave already exists (and the /O flag was not specified)");
		}
	}
}

void IgorImageOutputWriter::write_image(ImagePtr imageToWrite) {
	
	int result;
	size_t xSize = imageToWrite->rows();
	size_t ySize = imageToWrite->cols();
	size_t nPixels = xSize * ySize;
	int storage = this->GetIgorStorageType();
	
	if (this->outputWave == NULL) {
		// the outputwave has not been created yet, do it now
		CountInt dimensionSizes[MAX_DIMENSIONS + 1];
		
		dimensionSizes[0] = xSize;
		dimensionSizes[1] = ySize;
		dimensionSizes[2] = this->nImagesTotal;
		dimensionSizes[3] = 0;
		
		// the way to make the wave depends on whether this object was constructed with a full path
		// or with a DataFolderAndName argument
		if (this->fullPathToWave.length() != 0) {
			this->outputWave = MakeWaveUsingFullPath(this->fullPathToWave, dimensionSizes, storage, this->overwrite);
		} else {
			result = MDMakeWave(&(this->outputWave), this->waveDataFolderAndName.name, this->waveDataFolderAndName.dfH, dimensionSizes, storage, this->overwrite);
			if (result != 0)
				throw result;
		}
	}
	
	// check that we are not trying to write too many images
	// which would otherwise trigger an out-of-bounds memory
	// access
    if (this->nImagesWritten >= this->nImagesTotal)
		throw std::runtime_error("Writing too many images to the IgorImageOutputWriter");
	
	// the strategy for writing the data depends on the storage type
	size_t waveDataOffset;
	char *waveDataPtr;
	
	result = MDAccessNumericWaveData(this->outputWave, kMDWaveAccessMode0, (BCInt*)&waveDataOffset);
	if (result != 0)
		throw result;
	
	waveDataPtr = (char *)((char*)(*this->outputWave) + waveDataOffset);
	double* imagePtr = imageToWrite->data();
	
	switch (storage) {
		case NT_I8:
		{
			int8_t* int8Ptr = (int8_t*)waveDataPtr;
			int8Ptr += nPixels * this->nImagesWritten;
			for (size_t i = 0; i < nPixels; ++i) {
				int8Ptr[i] = imagePtr[i];
			}
			break;
		}
		case NT_I16:
		{
			int16_t* int16Ptr = (int16_t*)waveDataPtr;
			int16Ptr += nPixels * this->nImagesWritten;
			for (size_t i = 0; i < nPixels; ++i) {
				int16Ptr[i] = imagePtr[i];
			}
			break;
		}
		case NT_I32:
		{
			int32_t* int32Ptr = (int32_t*)waveDataPtr;
			int32Ptr += nPixels * this->nImagesWritten;
			for (size_t i = 0; i < nPixels; ++i) {
				int32Ptr[i] = imagePtr[i];
			}
			break;
		}
		case NT_I8 | NT_UNSIGNED:
		{
			uint8_t* uint8Ptr = (uint8_t*)waveDataPtr;
			uint8Ptr += nPixels * this->nImagesWritten;
			for (size_t i = 0; i < nPixels; ++i) {
				uint8Ptr[i] = imagePtr[i];
			}
			break;
		}
		case NT_I16 | NT_UNSIGNED:
		{
			uint16_t* uint16Ptr = (uint16_t*)waveDataPtr;
			uint16Ptr += nPixels * this->nImagesWritten;
			for (size_t i = 0; i < nPixels; ++i) {
				uint16Ptr[i] = imagePtr[i];
			}
			break;
		}
		case NT_I32 | NT_UNSIGNED:
		{
			uint32_t* uint32Ptr = (uint32_t*)waveDataPtr;
			uint32Ptr += nPixels * this->nImagesWritten;
			for (size_t i = 0; i < nPixels; ++i) {
				uint32Ptr[i] = imagePtr[i];
			}
			break;
		}
		case NT_FP32:
		{
			float* floatPtr = (float*)waveDataPtr;
			floatPtr += nPixels * this->nImagesWritten;
			for (size_t i = 0; i < nPixels; ++i) {
				floatPtr[i] = imagePtr[i];
			}
			break;
		}
		case NT_FP64:
		{
			double* doublePtr = (double*)waveDataPtr;
			doublePtr += nPixels * this->nImagesWritten;
			memcpy(doublePtr, imagePtr, nPixels * sizeof(double));
			break;
		}
	}
	
	++nImagesWritten;
}

int IgorImageOutputWriter::GetIgorStorageType() {
	int storage;
	
	switch (this->storageType) {
		case STORAGE_TYPE_INT4:
		case STORAGE_TYPE_UINT4:
		case STORAGE_TYPE_INT8:
			storage = NT_I8;
			break;
		case STORAGE_TYPE_UINT8:
			storage = NT_I8 | NT_UNSIGNED;
			break;
		case STORAGE_TYPE_INT16:
			storage = NT_I16;
			break;
		case STORAGE_TYPE_UINT16:
			storage = NT_I16 | NT_UNSIGNED;
			break;
		case STORAGE_TYPE_INT32:
			storage = NT_I32;
			break;
		case STORAGE_TYPE_UINT32:
			storage = NT_I32 | NT_UNSIGNED;
			break;
        case STORAGE_TYPE_INT64:
            storage = NT_I32;
			break;
		case STORAGE_TYPE_UINT64:
			storage = NT_I32 | NT_UNSIGNED;
			break;
		case STORAGE_TYPE_FP32:
			storage = NT_FP32;
			break;
		case STORAGE_TYPE_FP64:
			storage = NT_FP64;
			break;
	}
	
	return storage;
}

#endif // WITH_IGOR
