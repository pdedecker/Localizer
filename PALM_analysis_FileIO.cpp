/*
 *  PALM_analysis_FileIO.cpp
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 30/03/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "PALM_analysis_FileIO.h"

uint16 getUINT16FromCharArray(char *array, size_t offset) {
	char byteReader1, byteReader2;
	uint16 result;
	
	byteReader1 = array[offset];
	byteReader2 = array[offset + 1];
	
	result = 0xFF & byteReader2;
	result *= 256;
	result = result | (0x000000FF & byteReader1);
	
	return result;
}

uint32 getUINT32FromCharArray(char *array, size_t offset) {
	char byteReader1, byteReader2, byteReader3, byteReader4;
	uint32 result;
	
	byteReader1 = array[offset];
	byteReader2 = array[offset + 1];
	byteReader3 = array[offset + 2];
	byteReader4 = array[offset + 3];
	
	result = 0xFF & byteReader4;
	result *= 256;
	result = result | (0x000000FF & byteReader3);
	result *= 256;
	result = result | (0x000000FF & byteReader2);
	result *= 256;
	result = result | (0x000000FF & byteReader1);
	
	return result;
}

#ifdef WIN32
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
    int err;
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
    int err;
    fgets(buffer, nMax, this->fileRef);
    if (ferror(fileRef) != 0) {
		std::string error;
		error = "Error returned using getline() on the image file at \"";
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
#endif


ImageLoader::ImageLoader() {
	total_number_of_images = 0;
	x_size = 0;
	y_size = 0;
	header_length = 0;
}

ImageLoader::~ImageLoader() {
	if (file.is_open() == 1)
		file.close();
}

boost::shared_ptr<ublas::matrix<double> > ImageLoader::get_nth_image(const size_t n) {
	std::vector<boost::shared_ptr<ublas::matrix <double> > > images;
	if (n >= total_number_of_images)
		throw IMAGE_INDEX_BEYOND_N_IMAGES(std::string("Requested more images than there are in the file"));
	
	size_t firstImageToLoad = n;
	size_t lastImageToLoad = n;
	
	images = ReadImagesFromDisk(firstImageToLoad, lastImageToLoad);
	
	return images.at(0);
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
	long current_bytes = 0;	// assume that this is a four-byte variable
	char header_buffer[1500];
	char byte_reader1, byte_reader2, byte_reader3, byte_reader4;	// we want to only read 1 byte at a time
	
	// read the entire header into a buffer
	file.read(header_buffer, 1500);
	
	byte_reader1 = header_buffer[42];
	byte_reader2 = header_buffer[43];
	current_bytes = 0x000000FF & byte_reader2;
	current_bytes *= 256;
	current_bytes = current_bytes | (0x000000FF & byte_reader1);
	
	x_size = current_bytes;
	current_bytes = 0;
	
	byte_reader1 = header_buffer[656];
	byte_reader2 = header_buffer[657];
	current_bytes = 0x000000FF & byte_reader2;
	current_bytes *= 256;
	current_bytes = current_bytes | (0x000000FF & byte_reader1);
	
	y_size = current_bytes;
	current_bytes = 0;
	
	byte_reader1 = header_buffer[1446];
	byte_reader2 = header_buffer[1447];
	byte_reader3 = header_buffer[1448];
	byte_reader4 = header_buffer[1449];
	current_bytes = 0x000000FF  & byte_reader4;
	current_bytes *= 256;
	current_bytes = current_bytes | (0x000000FF & byte_reader3);
	current_bytes *= 256;
	current_bytes = current_bytes | (0x000000FF & byte_reader2);
	current_bytes *= 256;
	current_bytes = current_bytes | (0x000000FF & byte_reader1);
	
	total_number_of_images = current_bytes;
	current_bytes = 0;
	
	byte_reader1 = header_buffer[108];
	byte_reader2 = header_buffer[109];
	current_bytes = 0x000000FF & byte_reader2;
	current_bytes *= 256;
	current_bytes = current_bytes | (0x000000FF & byte_reader1);
	
	storage_type = (int)current_bytes;
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
	current_bytes = 0;
	
	// was there an error sometime during this procedure that would have caused the reading to fail?
	if (file.fail() != 0) {
		std::string error;
		error = "Error parsing the header information in \"";
		error += this->filePath;
		error += "\" assuming the SPE format";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	file.seekg(0);
}

std::vector<boost::shared_ptr<ublas::matrix <double> > > ImageLoaderSPE::ReadImagesFromDisk(size_t const nStart, size_t const nEnd) {	
	uint64_t offset;
	long current_long = 0;
	float current_float = 0;
	short current_short = 0;
	unsigned short current_unsigned_short = 0;
	std::string error;
	boost::shared_ptr<ublas::matrix<double> > image;
	std::vector<boost::shared_ptr<ublas::matrix <double> > > requestedImages;
	
	uint64_t n_bytes_in_single_image;
	uint64_t cache_offset;
	
	loadImagesMutex.lock();
	
	// determine how big we have to make the single image buffer
	switch(storage_type) {
		case STORAGE_TYPE_FP32:	// 4 byte float
			n_bytes_in_single_image = x_size * y_size * 4;
			break;
		case STORAGE_TYPE_UINT32:	// 4-byte long
			n_bytes_in_single_image = x_size * y_size * 4;
			break;
		case STORAGE_TYPE_INT16:	// 2 byte signed short
			n_bytes_in_single_image = x_size * y_size * 2;
			break;
		case STORAGE_TYPE_UINT16:	// 2 byte unsigned short
			n_bytes_in_single_image = x_size * y_size * 2;
			break;
		default:
			loadImagesMutex.unlock();
			std::string error("Unable to determine the storage type used in ");
			error += this->filePath;
			throw CANNOT_DETERMINE_SPE_STORAGE_TYPE(error);
			break;
	}
	
	boost::scoped_array<float> single_image_buffer_float(new float[x_size * y_size]);
	boost::scoped_array<char> single_image_buffer(new char[n_bytes_in_single_image]);
	
	for (uint64_t i = nStart; i <= nEnd; i++) {
		
		image = boost::shared_ptr<ublas::matrix<double> >(new ublas::matrix<double>(x_size, y_size));
		
		switch(storage_type) {
			case STORAGE_TYPE_FP32:	// 4 byte float
				offset = header_length + i * (x_size) * (y_size) * 4;
				break;
			case STORAGE_TYPE_UINT32:	// 4-byte long
				offset = header_length + i * (x_size) * (y_size) * 4;
				break;
			case STORAGE_TYPE_INT16:	// 2 byte signed short
				offset = header_length + i * (x_size) * (y_size) * 2;
				break;
			case STORAGE_TYPE_UINT16:	// 2 byte unsigned short
				offset = header_length + i * (x_size) * (y_size) * 2;
				break;
			default:
				loadImagesMutex.unlock();
				std::string error("Unable to determine the storage type used in ");
				error += this->filePath;
				throw CANNOT_DETERMINE_SPE_STORAGE_TYPE(error);
				break;
		}
		
		
		file.seekg(offset);
		cache_offset = 0;
		
		file.read((char *)single_image_buffer.get(), n_bytes_in_single_image);
		if (file.fail() != 0) {
			std::string error;
			error = "Error trying to read image data from \"";
			error += this->filePath;
			error += "\" assuming the SPE format";
			loadImagesMutex.unlock();
			throw ERROR_READING_FILE_DATA(error);
		}
		
		switch(storage_type) {
			case STORAGE_TYPE_FP32:	// 4-byte float
				// this is currently only safe on little-endian systems!
				for (size_t j  = 0; j < y_size; j++) {
					for (size_t i = 0; i < x_size; i++) {
						current_float = single_image_buffer_float[cache_offset];
						
						(*image)(i, j) = (double)current_float;
						
						cache_offset++;
					}
				}
				break;
				
			case STORAGE_TYPE_UINT32:	// 4-byte long
				for (size_t j  = 0; j < y_size; j++) {
					for (size_t i = 0; i < x_size; i++) {
						current_long = 0x000000FF & single_image_buffer[cache_offset + 3];	// little endian
						current_long *= 256;
						current_long = current_long | (0x000000FF & single_image_buffer[cache_offset + 2]);
						current_long *= 256;
						current_long = current_long | (0x000000FF & single_image_buffer[cache_offset + 1]);
						current_long *= 256;
						current_long = current_long | (0x000000FF & single_image_buffer[cache_offset + 0]);
						
						(*image)(i, j) = (double)current_long;
						
						cache_offset += 4;
					}
				}
				break;
				
			case STORAGE_TYPE_INT16:	// 2-byte signed short
				for (size_t j  = 0; j < y_size; j++) {
					for (size_t i = 0; i < x_size; i++) {
						current_short = 0x000000FF & single_image_buffer[cache_offset + 1];	// little endian
						current_short *= 256;
						current_short = current_short | (0x000000FF & single_image_buffer[cache_offset + 0]);
						
						(*image)(i, j) = (double)current_short;
						
						cache_offset += 2;
					}
				}
				break;
				
			case STORAGE_TYPE_UINT16: // 2-byte unsigned short
				for (size_t j  = 0; j < y_size; j++) {
					for (size_t i = 0; i < x_size; i++) {
						current_unsigned_short = 0x000000FF & single_image_buffer[cache_offset + 1];	// little endian
						current_unsigned_short *= 256;
						current_unsigned_short = current_unsigned_short | (0x000000FF & single_image_buffer[cache_offset + 0]);
						
						(*image)(i, j) = (double)current_unsigned_short;
						
						cache_offset += 2;
						
						//	current_unsigned_short = 0;
					}
				}
				break;
				
			default:
				loadImagesMutex.unlock();
				std::string error("Unable to determine the storage type used in ");
				error += this->filePath;
				throw CANNOT_DETERMINE_SPE_STORAGE_TYPE(error);
				break;
				
		}
		
		requestedImages.push_back(image);
	}
	
	file.seekg(0);
	
	loadImagesMutex.unlock();
	
	return requestedImages;
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
	
	result = sscanf(singleLine.c_str(), "Pixel number%zu 1 %zu %zu 1 %zu", &temp, &frameYSize, &frameXSize, &this->total_number_of_images);
	if (result != 4)
		throw std::runtime_error(std::string("an error occured parsing the file assuming the Andor format"));
	
	ss.getline(singleLineBuffer.get(), 4096);
	singleLine = singleLineBuffer.get();
	
	
	result = sscanf(singleLine.c_str(), "%zu %zu %zu %zu %zu %zu %zu", &temp, &xStart, &yEnd, &xEnd, &yStart, &xBinning, &yBinning);
	if (result != 7)
		throw std::runtime_error(std::string("an error occured parsing the file assuming the Andor format"));
	
	this->x_size = (xEnd - xStart + 1) / xBinning;
	this->y_size = (yEnd - yStart + 1) / yBinning;	// integer division
	
	// now there are some lines that may contain timestamps. There are as many lines as there are images
	for (size_t i = 0;  i < total_number_of_images; i++) {
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
}

std::vector<boost::shared_ptr<ublas::matrix <double> > > ImageLoaderAndor::ReadImagesFromDisk(size_t const nStart, size_t const nEnd) {	
	uint64_t offset;	// off_t is the size of the file pointer used by the OS
	float current_float = 0;
	
	boost::shared_ptr<ublas::matrix<double> > image;
	std::vector<boost::shared_ptr<ublas::matrix <double> > > requestedImages;
	
	loadImagesMutex.lock();
	
	boost::scoped_array<float> single_image_buffer(new float[x_size * y_size]);
	uint64_t cache_offset;
	
	for (uint64_t i = nStart; i <= nEnd; i++) {
		image = boost::shared_ptr<ublas::matrix<double> >(new ublas::matrix<double>(x_size, y_size));
		
		offset = header_length + i * (x_size) * (y_size) * sizeof(float);
		file.seekg(offset);
		
		cache_offset = 0;
		
		file.read((char *)single_image_buffer.get(), (x_size * y_size * sizeof(float)));
		if (file.fail() != 0) {
			std::string error;
			error = "Error trying to read image data from \"";
			error += this->filePath;
			error += "\" assuming the Andor format";
			loadImagesMutex.unlock();
			throw ERROR_READING_FILE_DATA(error);
		}
		
		
		// this is currently only safe on little-endian systems!
		for (size_t j  = 0; j < y_size; j++) {
			for (size_t i = 0; i < x_size; i++) {
				current_float = single_image_buffer[cache_offset];
				(*image)(i, j) = (double)current_float;
				cache_offset++;
			}
		}
		
		requestedImages.push_back(image);
	}
	
	loadImagesMutex.unlock();
	
	return requestedImages;
}

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
	char headerBuffer[64];	// too large, but we'd better be safe
	ImageLoaderHamamatsu_HeaderStructure header;
	
	storage_type = STORAGE_TYPE_UINT16;
	
	// first we load a set of data
	file.seekg(0);
	file.read(headerBuffer, 63);
	headerBuffer[63] = '\0';
	
	
	// from the char array build up a complete header structure
	header.magic = getUINT16FromCharArray(headerBuffer, 0);
	header.commentLength = getUINT16FromCharArray(headerBuffer, 2);
	header.xSize = getUINT16FromCharArray(headerBuffer, 4);
	header.ySize = getUINT16FromCharArray(headerBuffer, 6);
	header.xBinning = getUINT16FromCharArray(headerBuffer, 8);
	header.yBinning = getUINT16FromCharArray(headerBuffer, 10);
	header.storageFormat = getUINT16FromCharArray(headerBuffer, 12);
	header.nImages = getUINT32FromCharArray(headerBuffer, 14);
	header.nChannels = getUINT16FromCharArray(headerBuffer, 18);
	header.channel = getUINT16FromCharArray(headerBuffer, 20);
	// skip the timestamp for now
	header.marker = getUINT32FromCharArray(headerBuffer, 28);
	header.misc = getUINT32FromCharArray(headerBuffer, 32);
	
	if (header.storageFormat != 2) {	// not UINT16
		std::string error;
		error = "The file at \"";
		error += this->filePath;
		error += "\" specifies that it doesn't use UINT16 for storage. This usually means that the manufacturer's software corrupted the file.";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	header_length = header.commentLength + 64;
	x_size = header.xSize;
	y_size = header.ySize;
	total_number_of_images = header.nImages;
	
	// was there an error reading the file?
	if (file.fail() != 0) {
		std::string error;
		error = "Error parsing the header information in \"";
		error += this->filePath;
		error += "\" assuming the Hamamatsu format";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	file.seekg(0);
}


std::vector<boost::shared_ptr<ublas::matrix <double> > > ImageLoaderHamamatsu::ReadImagesFromDisk(size_t const nStart, size_t const nEnd) {	
	uint64_t offset;	// off_t is the size of the file pointer used by the OS
	
	boost::shared_ptr<ublas::matrix<double> > image;
	std::vector<boost::shared_ptr<ublas::matrix <double> > > requestedImages;
	
	loadImagesMutex.lock();
	
	unsigned int current_uint;
	size_t n_bytes_per_image = x_size * y_size * 2;
	boost::scoped_array<char> single_image_buffer(new char[n_bytes_per_image]);
	uint64_t cache_offset;
	
	for (uint64_t i = nStart; i <= nEnd; i++) {
		image = boost::shared_ptr<ublas::matrix<double> >(new ublas::matrix<double>(x_size, y_size));
		
		offset = (i + 1) * header_length + i * (x_size) * (y_size) * 2;	// assume a 16-bit format
		file.seekg(offset);
		
		file.read(single_image_buffer.get(), n_bytes_per_image);
		if (file.fail() != 0) {
			std::string error;
			error = "Error reading image data from \"";
			error += this->filePath;
			error += "\" assuming the Hamamatsu format";
			loadImagesMutex.unlock();
			throw ERROR_READING_FILE_DATA(error);
		}
		
		cache_offset = 0;
		
		for (size_t j  = 0; j < y_size; j++) {
			for (size_t i = 0; i < x_size; i++) {
				current_uint = 0x000000FF & single_image_buffer[cache_offset + 1];	// little endian
				current_uint *= 256;
				current_uint = current_uint | (0x000000FF & single_image_buffer[cache_offset]);
				
				(*image)(i, j) = (double)current_uint;
				
				cache_offset += 2;	// TWO bytes per value (UINT16)
			}
		}
		
		
		requestedImages.push_back(image);
	}
	
	file.seekg(0);
	
	loadImagesMutex.unlock();
	
	return requestedImages;
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
	file.seekg(0);
	
	PDEFormatHeader header;
	
	file.read((char *)&header, sizeof(header));
	
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
		throw std::runtime_error("Unsupported version of the PDE format. Perhaps you need to upgrade?");
	
	this->total_number_of_images = header.nImages;
	this->x_size = header.xSize;
	this->y_size = header.ySize;
	this->storage_type = header.storageFormat;
	this->header_length = sizeof(header);
}

std::vector<boost::shared_ptr<ublas::matrix <double> > > ImageLoaderPDE::ReadImagesFromDisk(size_t const nStart, size_t const nEnd) {
	uint64_t offset;
	boost::shared_ptr<ublas::matrix<double> > image;
	std::vector<boost::shared_ptr<ublas::matrix <double> > > requestedImages;
	size_t n_pixels = this->x_size * this->y_size;
	
	loadImagesMutex.lock();
	
	for (uint64_t i = nStart; i <= nEnd; i++) {
		image = boost::shared_ptr<ublas::matrix<double> >(new ublas::matrix<double>(this->x_size, this->y_size));
		
		switch (this->storage_type) {
			case STORAGE_TYPE_UINT16:
			{
				offset = header_length + i * n_pixels * sizeof(uint16_t);
				boost::scoped_array<uint16_t> buffer(new uint16_t[n_pixels]);
				file.seekg(offset);
				file.read((char *)buffer.get(), n_pixels * sizeof(uint16_t));
				offset = 0;
				for (ublas::matrix<double>::array_type::iterator it = image->data().begin(); it != image->data().end(); ++it) {
					*it = buffer[offset];
					++offset;
				}
				break;
			}
			case STORAGE_TYPE_UINT32:
			{
				offset = header_length + i * n_pixels * sizeof(uint32_t);
				boost::scoped_array<uint32_t> buffer(new uint32_t[n_pixels]);
				file.seekg(offset);
				file.read((char *)buffer.get(), n_pixels * sizeof(uint32_t));
				offset = 0;
				for (ublas::matrix<double>::array_type::iterator it = image->data().begin(); it != image->data().end(); ++it) {
					*it = buffer[offset];
					++offset;
				}
				break;
			}
			case STORAGE_TYPE_FP32:
			{
				offset = header_length + i * n_pixels * sizeof(float);
				boost::scoped_array<float> buffer(new float[n_pixels]);
				file.seekg(offset);
				file.read((char *)buffer.get(), n_pixels * sizeof(float));
				offset = 0;
				for (ublas::matrix<double>::array_type::iterator it = image->data().begin(); it != image->data().end(); ++it) {
					*it = buffer[offset];
					++offset;
				}
				break;
			}
			case STORAGE_TYPE_FP64:
			{
				offset = header_length + i * n_pixels * sizeof(double);
				boost::scoped_array<double> buffer(new double[n_pixels]);
				file.seekg(offset);
				file.read((char *)buffer.get(), n_pixels * sizeof(double));
				offset = 0;
				for (ublas::matrix<double>::array_type::iterator it = image->data().begin(); it != image->data().end(); ++it) {
					*it = buffer[offset];
					++offset;
				}
				break;
			}
			default:
				throw std::runtime_error("The data file does appear to contain a recognized storage type");
				break;
		}
		
		if (file.fail() != 0) {
			std::string error;
			error = "Error reading image data from \"";
			error += this->filePath;
			error += "\" assuming the simple image format";
			loadImagesMutex.unlock();
			throw ERROR_READING_FILE_DATA(error);
		}
		
		requestedImages.push_back(image);
	}
	
	file.seekg(0);
	
	loadImagesMutex.unlock();
	
	return requestedImages;
}

ImageLoaderTIFF::ImageLoaderTIFF(std::string rhs) {
	this->filePath = rhs;
	
	tiff_file = NULL;
	
	tiff_file = TIFFOpen(this->filePath.c_str(), "r");
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
	size_t index = 0;
	
	
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
	
	x_size = (size_t)(result_uint32);
	
	// what is the y size?
	result = TIFFGetField(tiff_file, TIFFTAG_IMAGELENGTH, &result_uint32);
	if (result != 1) {
		std::string error;
		error = "The image at\"";
		error += this->filePath;
		error += "\" does not specify a height";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	y_size = (size_t)(result_uint32);
	
	
	// how many images are there in the file?
	// the first image needs to be treated as a special case
	result = TIFFGetField(tiff_file, TIFFTAG_SUBFILETYPE, &result_uint32);
	if (((result_uint32 == FILETYPE_REDUCEDIMAGE) || (result_uint32 == FILETYPE_MASK)) && (result == 1)) {
		++index;
	} else {
		directoryIndices.push_back(index);
		++index;
	}
	
	while (TIFFReadDirectory(tiff_file) != 0) {
		result = TIFFGetField(tiff_file, TIFFTAG_SUBFILETYPE, &result_uint32);
		if (((result_uint32 == FILETYPE_REDUCEDIMAGE) || (result_uint32 == FILETYPE_MASK)) && (result == 1)) {
			++index;
			continue;
		}
		
		directoryIndices.push_back(index);
		
		++index;
	}
	
	
	total_number_of_images = directoryIndices.size();
	
	result = TIFFSetDirectory(tiff_file, 0);
	if (result != 1) {
		std::string error;
		error = "Unable to set the directory to '0' for the image at\"";
		error += this->filePath;
		error += "\"";
		throw ERROR_READING_FILE_DATA(error);
	}
}

std::vector<boost::shared_ptr<ublas::matrix <double> > > ImageLoaderTIFF::ReadImagesFromDisk(size_t const nStart, size_t const nEnd) {
	char *single_scanline_buffer;
	char *scanline_pointer;
	uint16_t current_uint16;
	uint32_t current_uint32;
	float current_float;
	double current_double;
	uint16_t *uint16Ptr;
	uint32_t *uint32Ptr;
	float *floatPtr;
	double *doublePtr;
	boost::shared_ptr<ublas::matrix<double> > image;
	std::vector<boost::shared_ptr<ublas::matrix <double> > > requestedImages;
	int result;
	
	loadImagesMutex.lock();
	
	single_scanline_buffer = (char *)_TIFFmalloc(TIFFScanlineSize(tiff_file));
	if (single_scanline_buffer == NULL) {
		throw std::bad_alloc();
	}
	
	for (size_t i = nStart; i <= nEnd; i++) {
		try {
			image = boost::shared_ptr<ublas::matrix<double> >(new ublas::matrix<double>(x_size, y_size));
		}
		catch (std::bad_alloc) {
			_TIFFfree(single_scanline_buffer);
			loadImagesMutex.unlock();
			throw std::bad_alloc();
		}
		
		result = TIFFSetDirectory(tiff_file, directoryIndices.at(i));
		if (result != 1) {
			_TIFFfree(single_scanline_buffer);
			std::string error;
			error = "Unable to set the directory to '0' for the image at\"";
			error += this->filePath;
			error += "\"";
			loadImagesMutex.unlock();
			throw ERROR_READING_FILE_DATA(error);
		}
		
		for (size_t j = 0; j < y_size; ++j) {
			result = TIFFReadScanline(tiff_file, single_scanline_buffer, j, 0);	// sample is ignored
			if (result != 1) {
				_TIFFfree(single_scanline_buffer);
				std::string error;
				error = "Unable to read a scanline from the image at\"";
				error += this->filePath;
				error += "\"";
				loadImagesMutex.unlock();
				throw ERROR_READING_FILE_DATA(error);
			}
			
			switch (storage_type) {	// handle the different possibilities (floating, integer) and variable sizes
				case STORAGE_TYPE_UINT4:
					scanline_pointer = single_scanline_buffer;
					for (size_t k = 0; k < x_size; ++k) {
						if ((k % 2) == 0) {	// this is an even pixel, we use only the first 4 bits
							current_uint16 = 0x0000000F & (*scanline_pointer);
							(*image)(k, j) = (double)current_uint16;
						} else {	// this is an odd pixel, use the last 4 bits and increment the scanline_pointer
							current_uint16 = 0x000000F0 & (*scanline_pointer);
							(*image)(k, j) = (double)current_uint16;
							scanline_pointer += 1;
						}
					}
					break;
					
				case STORAGE_TYPE_UINT8:
					scanline_pointer = single_scanline_buffer;
					for (size_t k = 0; k < x_size; ++k) {
						current_uint16 = (uint16_t)(*scanline_pointer);
						(*image)(k, j) = (double)current_uint16;
						scanline_pointer += 1;
					}
					break;
					
				case STORAGE_TYPE_UINT16:
					uint16Ptr = (uint16_t*)single_scanline_buffer;
					for (size_t k = 0; k < x_size; ++k) {
						current_uint16 = (*uint16Ptr);
						(*image)(k, j) = (double)current_uint16;
						uint16Ptr += 1;
					}
					break;
					
				case STORAGE_TYPE_UINT32:
					uint32Ptr = (uint32_t*)single_scanline_buffer;
					for (size_t k = 0; k < x_size; ++k) {
						current_uint32 = (*uint32Ptr);
						(*image)(k, j) = (double)current_uint32;
						uint32Ptr += 1;
					}
					break;
					
				case STORAGE_TYPE_FP32:
					floatPtr = (float *)single_scanline_buffer;
					for (size_t k = 0; k < x_size; ++k) {
						current_float = *floatPtr;
						(*image)(k, j) = (double)current_float;
						floatPtr += 1;
					}
					break;
					
				case STORAGE_TYPE_FP64:
					doublePtr = (double *)single_scanline_buffer;
					for (size_t k = 0; k < x_size; ++k) {
						current_double = *doublePtr;
						(*image)(k, j) = current_double;
						doublePtr += 1;
					}
					break;
					
				default:
					_TIFFfree(single_scanline_buffer);
					std::string error;
					error = "Invalid floating point data size for the image at\"";
					error += this->filePath;
					error += "\"";
					loadImagesMutex.unlock();
					throw ERROR_READING_FILE_DATA(error);
					break;
			}
		}
		
		requestedImages.push_back(image);
	}
	
	_TIFFfree(single_scanline_buffer);
	result = TIFFSetDirectory(tiff_file, 0);
	if (result != 1) {
		std::string error;
		error = "Invalid to set the directory to '0' for the image at\"";
		error += this->filePath;
		error += "\"";
		loadImagesMutex.unlock();
		throw ERROR_READING_FILE_DATA(error);
	}
	
	loadImagesMutex.unlock();
	
	return requestedImages;
}

#ifdef WITH_IGOR
ImageLoaderIgor::ImageLoaderIgor(std::string waveName) {
	int waveType;
	
	// try to get images from an Igor wave
	
	this->igor_data_wave = FetchWaveUsingFullPath(waveName);
	
	long DimensionSizes[MAX_DIMENSIONS + 1];
	long numDimensions;
	int result;
	
	result = MDGetWaveDimensions(igor_data_wave, &numDimensions, DimensionSizes);
	if (result != 0) {
		throw result;
	}
	if ((numDimensions != 2) && (numDimensions != 3)) {
		throw INCOMPATIBLE_DIMENSIONING;
	}
	
	x_size = (size_t)DimensionSizes[0];
	y_size = (size_t)DimensionSizes[1];
	total_number_of_images = (size_t)DimensionSizes[2];
	
	// special case: if the wave contains only a single image then it is usually two-dimensional, that is, DimensionSizes[2] == 0
	// in that case total_number_of_images is still one
	if (DimensionSizes[2] == 0) {
		total_number_of_images = 1;
	}
	
	waveType = WaveType(igor_data_wave);
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
}

std::vector<boost::shared_ptr<ublas::matrix <double> > > ImageLoaderIgor::ReadImagesFromDisk(size_t const nStart, size_t const nEnd) {
	double value[2];
	long indices[3];
	int result;
	
	boost::shared_ptr<ublas::matrix<double> > image;
	std::vector<boost::shared_ptr<ublas::matrix <double> > > requestedImages;
	
	// no mutex locking is required since these calls are all threadsafe
	
	for (size_t n = nStart; n <= nEnd; ++n) {
		indices[2] = n;
		image = boost::shared_ptr<ublas::matrix<double> > (new ublas::matrix<double>(x_size, y_size));
		
		for (size_t i = 0; i < x_size; i++) {
			for (size_t j  = 0; j < y_size; j++) {
				indices[0] = (long)i;
				indices[1] = (long)j;
				
				result = MDGetNumericWavePointValue(igor_data_wave, indices, value);
				if (result != 0) {
					throw result;
				}
				
				(*image)(i, j) = value[0];
			}
		}
		
		requestedImages.push_back(image);
		
	}
	
	return requestedImages;
}
#endif // WITH_IGOR



ImageOutputWriter::ImageOutputWriter() {
	file_path.assign("");
	n_images_written = 0;
}

ImageOutputWriter::ImageOutputWriter(const std::string &rhs, int overwrite) {
	file_path.assign("");
	n_images_written = 0;
}



PDEImageOutputWriter::PDEImageOutputWriter(const std::string &rhs,int overwrite, uint32_t storageType_rhs) {
	// if overwrite is non-zero then we overwrite any file that exists at the output path
	// if it is set to zero then we throw an error and abort instead of overwriting
	file_path = rhs;
	int header_length = 4 * sizeof(uint32_t);
	
	if (overwrite == 0) {
		std::ifstream input_test;
		input_test.open(file_path.c_str(), std::ios::in | std::ios::binary);
		input_test.close();
		if (input_test.fail() == 0) {
			std::string error("The output file at ");
			error += this->file_path;
			throw OUTPUT_FILE_ALREADY_EXISTS(error);	// escape without overwriting
		} else {
			file.open(file_path.c_str(), std::ios::binary | std::ios::out);
		}
	}
	
	if (overwrite != 0) {
		file.open(file_path.c_str(), std::ios::binary | std::ios::out | std::ios::trunc);	// DANGER: OVERWRITING THE FILE
	}
	
	if (file.fail() != 0) {
		std::string error("Cannot create an output file at ");
		error += this->file_path;
		throw CANNOT_OPEN_OUTPUT_FILE(error);
	}
	
	PDEFormatHeader header;
	header.magic = 0;
	header.version = 0;
	header.nImages = 0;
	header.xSize = 0;
	header.ySize = 0;
	header.storageFormat = 0;
	
	this->x_size = 0;
	this->y_size = 0;
	this->n_images_written = 0;
	this->storageType = storageType_rhs;
	
	// the first 3 * 4 bytes of the file should be written in advance, we will fill them in later
	// by convention they are x_size, y_size, and n_images
	file.write((char *)&header, sizeof(header));
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
	header.nImages = this->n_images_written;
	header.xSize = this->x_size;
	header.ySize = this->y_size;
	header.storageFormat = this->storageType;
	
	if (file.is_open()) {
		file.seekp(0, std::ios_base::beg);
		file.write((char *)&header, sizeof(header));
		if (file.fail())
			throw std::runtime_error("Error trying to write the header");
	}
}

void PDEImageOutputWriter::write_image(boost::shared_ptr<ublas::matrix<double> > imageToWrite) {
	
	// determine the size of the frames
	size_t x_size = imageToWrite->size1();
	size_t y_size = imageToWrite->size2();
	size_t n_pixels = x_size * y_size;
	
	size_t offset = 0;
	
	if (this->n_images_written == 0) {
		this->x_size = x_size;
		this->y_size = y_size;
	} else {
		assert((x_size == this->x_size) && (y_size == this->y_size));
	}
	
	switch (this->storageType) {
		case STORAGE_TYPE_UINT16:
		{
			boost::scoped_array<uint16_t> buffer(new uint16_t[n_pixels]);
			for (ublas::matrix<double>::array_type::const_iterator it = imageToWrite->data().begin(); it != imageToWrite->data().end(); ++it) {
				buffer[offset] = (uint16_t)(*it);
				++offset;
			}
			this->file.write((char *)buffer.get(), n_pixels * sizeof(uint16_t));
			break;
		}
		case STORAGE_TYPE_UINT32:
		{
			boost::scoped_array<uint32_t> buffer(new uint32_t[n_pixels]);
			for (ublas::matrix<double>::array_type::const_iterator it = imageToWrite->data().begin(); it != imageToWrite->data().end(); ++it) {
				buffer[offset] = (uint32_t)(*it);
				++offset;
			}
			this->file.write((char *)buffer.get(), n_pixels * sizeof(uint32_t));
			break;
		}
		case STORAGE_TYPE_FP32:
		{
			boost::scoped_array<float> buffer(new float[n_pixels]);
			for (ublas::matrix<double>::array_type::const_iterator it = imageToWrite->data().begin(); it != imageToWrite->data().end(); ++it) {
				buffer[offset] = (float)(*it);
				++offset;
			}
			this->file.write((char *)buffer.get(), n_pixels * sizeof(float));
			break;
		}
		case STORAGE_TYPE_FP64:
		{
			boost::scoped_array<double> buffer(new double[n_pixels]);
			for (ublas::matrix<double>::array_type::const_iterator it = imageToWrite->data().begin(); it != imageToWrite->data().end(); ++it) {
				buffer[offset] = (double)(*it);
				++offset;
			}
			this->file.write((char *)buffer.get(), n_pixels * sizeof(double));
			break;
		}
		default:
			throw std::runtime_error("Unsupport file type requested in the simple output format");
	}
	++this->n_images_written;
}

TIFFImageOutputWriter::TIFFImageOutputWriter(const std::string &rhs,int overwrite, int compression_rhs, int storageType_rhs) {
	// if overwrite is non-zero then we overwrite any file that exists at the output path
	// if it is set to zero then we throw an error and abort instead of overwriting
	
	this->file_path = rhs;
	this->compression = compression_rhs;
	this->storageType = storageType_rhs;
	
	if (overwrite == 0) {
		std::ifstream input_test;
		input_test.open(file_path.c_str(), std::ios::in | std::ios::binary);
		input_test.close();
		if (input_test.fail() == 0) {
			std::string error("The output file at ");
			error += this->file_path;
			error += " already exists";
			throw OUTPUT_FILE_ALREADY_EXISTS(error);	// escape without overwriting
		}
	}
	
	tiff_file = TIFFOpen(file_path.c_str(), "w");
	if (tiff_file == NULL) {
		std::string error("Cannot create an output file at ");
		error += this->file_path;
		throw CANNOT_OPEN_OUTPUT_FILE(error);
	}
	
	n_images_written = 0;
}


TIFFImageOutputWriter::~TIFFImageOutputWriter() {
	if (tiff_file != NULL) {
		TIFFClose(tiff_file);
		tiff_file = NULL;
	}
}



void TIFFImageOutputWriter::write_image(boost::shared_ptr<ublas::matrix<double> > imageToWrite) {
	
	size_t x_size = imageToWrite->size1();
	size_t y_size = imageToWrite->size2();
	
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
	boost::scoped_array<char> scanLine(new char[x_size * bitsPerSample]);
	
	// make sure that all the image tags have the correct values
	result = TIFFSetField(tiff_file, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
	if (result != 1) {
		std::string error;
		error = "Unable to set the photometric type for the image at\"";
		error += file_path;
		error += "\"";
		throw ERROR_WRITING_FILE_DATA(error);
	}
	
	current_uint32 = x_size;
	result = TIFFSetField(tiff_file, TIFFTAG_IMAGEWIDTH, current_uint32);
	if (result != 1) {
		std::string error;
		error = "Unable to set the image width for the image at\"";
		error += file_path;
		error += "\"";
		throw ERROR_WRITING_FILE_DATA(error);
	}
	
	current_uint32 = y_size;
	result = TIFFSetField(tiff_file, TIFFTAG_IMAGELENGTH, current_uint32);
	if (result != 1) {
		std::string error;
		error = "Unable to set the image height for the image at\"";
		error += file_path;
		error += "\"";
		throw ERROR_WRITING_FILE_DATA(error);
	}
	
	result = TIFFSetField(tiff_file, TIFFTAG_SAMPLEFORMAT, sampleFormat);
	if (result != 1) {
		std::string error;
		error = "Unable to set the SampleFormat for the image at\"";
		error += file_path;
		error += "\"";
		throw ERROR_WRITING_FILE_DATA(error);
	}
	
	result = TIFFSetField(tiff_file, TIFFTAG_BITSPERSAMPLE, bitsPerSample);
	if (result != 1) {
		std::string error;
		error = "Unable to set the BitsPerSample for the image at\"";
		error += file_path;
		error += "\"";
		throw ERROR_WRITING_FILE_DATA(error);
	}
	
	result = TIFFSetField(tiff_file, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
	if (result != 1) {
		std::string error;
		error = "Unable to set the SubFileType for the image at\"";
		error += file_path;
		error += "\"";
		throw ERROR_WRITING_FILE_DATA(error);
	}
	
	current_uint16 = (uint16_t)n_images_written;
	result = TIFFSetField(tiff_file, TIFFTAG_PAGENUMBER, current_uint16);
	if (result != 1) {
		std::string error;
		error = "Unable to set the PageNumber for the image at\"";
		error += file_path;
		error += "\"";
		throw ERROR_WRITING_FILE_DATA(error);
	}
	
	result = TIFFSetField(tiff_file, TIFFTAG_COMPRESSION, compression);
	if (result != 1) {
		std::string error;
		error = "Unable to set the compression method for the image at\"";
		error += file_path;
		error += "\"";
		throw ERROR_WRITING_FILE_DATA(error);
	}
	
	size_t offset = 0;
	for (size_t j = 0; j < y_size; j++) {
		offset = 0;
		
		switch (this->storageType) {
			case STORAGE_TYPE_INT4:
			case STORAGE_TYPE_UINT4:
			case STORAGE_TYPE_INT8:
				for (size_t i = 0; i < x_size; i++) {
					int8_t* buffer = (int8_t *)scanLine.get();
					buffer[offset] = (int8_t)(*imageToWrite)(i, j);
					offset++;
				}
				break;
			case STORAGE_TYPE_UINT8:
				for (size_t i = 0; i < x_size; i++) {
					uint8_t *buffer = (uint8_t *)scanLine.get();
					buffer[offset] = (uint8_t)(*imageToWrite)(i, j);
					offset++;
				}
				break;
			case STORAGE_TYPE_INT16:
				for (size_t i = 0; i < x_size; i++) {
					int16_t *buffer = (int16_t *)scanLine.get();
					buffer[offset] = (int16_t)(*imageToWrite)(i, j);
					offset++;
				}
				break;
			case STORAGE_TYPE_UINT16:
				for (size_t i = 0; i < x_size; i++) {
					uint16_t *buffer = (uint16_t *)scanLine.get();
					buffer[offset] = (uint16_t)(*imageToWrite)(i, j);
					offset++;
				}
				break;
			case STORAGE_TYPE_INT32:
				for (size_t i = 0; i < x_size; i++) {
					int32_t *buffer = (int32_t *)scanLine.get();
					buffer[offset] = (int32_t)(*imageToWrite)(i, j);
					offset++;
				}
				break;
			case STORAGE_TYPE_UINT32:
				for (size_t i = 0; i < x_size; i++) {
					uint32_t *buffer = (uint32_t *)scanLine.get();
					buffer[offset] = (uint32_t)(*imageToWrite)(i, j);
					offset++;
				}
				break;
			case STORAGE_TYPE_INT64:
				for (size_t i = 0; i < x_size; i++) {
					int64_t *buffer = (int64_t *)scanLine.get();
					buffer[offset] = (int64_t)(*imageToWrite)(i, j);
					offset++;
				}
				break;
			case STORAGE_TYPE_UINT64:
				for (size_t i = 0; i < x_size; i++) {
					uint64_t *buffer = (uint64_t *)scanLine.get();
					buffer[offset] = (uint64_t)(*imageToWrite)(i, j);
					offset++;
				}
				break;
			case STORAGE_TYPE_FP32:
				for (size_t i = 0; i < x_size; i++) {
					float *buffer = (float *)scanLine.get();
					buffer[offset] = (float)(*imageToWrite)(i, j);
					offset++;
				}
				break;
			case STORAGE_TYPE_FP64:
				for (size_t i = 0; i < x_size; i++) {
					double *buffer = (double *)scanLine.get();
					buffer[offset] = (double)(*imageToWrite)(i, j);
					offset++;
				}
				break;
		}
		
		result = TIFFWriteScanline(tiff_file, (char *)scanLine.get(), j);
		if (result != 1) {
			std::string error;
			error = "There was an error writing a scanline for the image at\"";
			error += file_path;
			error += "\"";
			throw ERROR_WRITING_FILE_DATA(error);
		}
		
	}
	
	result = TIFFWriteDirectory(tiff_file);
	if (result != 1) {
		std::string error;
		error = "Unable to write a directory for the image at\"";
		error += file_path;
		error += "\"";
		throw ERROR_WRITING_FILE_DATA(error);
	}
	
	++this->n_images_written;
}

#ifdef WITH_IGOR
IgorImageOutputWriter::IgorImageOutputWriter(std::string waveName_rhs, size_t nImages_rhs, int overwrite_rhs, int storageType_rhs) {
	this->outputWave = NULL;
	this->waveName = waveName_rhs;
	this->nImagesTotal = nImages_rhs;
	this->n_images_written = 0;
	this->overwrite = overwrite_rhs;
	this->storageType = storageType_rhs;
}

void IgorImageOutputWriter::write_image(boost::shared_ptr<ublas::matrix<double> > imageToWrite) {
	long indices[MAX_DIMENSIONS + 1];
	int result;
	double value[2];
	
	size_t x_size = imageToWrite->size1();
	size_t y_size = imageToWrite->size2();
	
	if (this->outputWave == NULL) {
		// the outputwave has not been created yet, do it now
		long dimensionSizes[MAX_DIMENSIONS + 1];
		int storage;
		
		dimensionSizes[0] = x_size;
		dimensionSizes[1] = y_size;
		dimensionSizes[2] = this->nImagesTotal;
		dimensionSizes[3] = 0;
		
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
				storage = NT_I32;	// todo: not yet supported in Igor
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
			default:
				throw std::runtime_error("Unsupported output format in IgorImageOutputWriter");
		}
		
		this->outputWave = MakeWaveUsingFullPath(this->waveName, dimensionSizes, storage, this->overwrite);
	}
	
	indices[2] = n_images_written;
	
	for (size_t i = 0; i < x_size; i++) {
		for (size_t j  = 0; j < y_size; j++) {
			indices[0] = i;
			indices[1] = j;
			
			value[0] = (*imageToWrite)(i,j);
			result = MDSetNumericWavePointValue(outputWave, indices, value);
			if (result != 0) {
				throw result;
			}
		}
	}
	
	++n_images_written;
}
#endif // WITH_IGOR
