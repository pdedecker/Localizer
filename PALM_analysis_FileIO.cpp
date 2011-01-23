/*
 *  PALM_analysis_FileIO.cpp
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 30/03/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "PALM_analysis_FileIO.h"

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

void ImageLoader::checkForReasonableValues() {
	if ((this->x_size > kMaxImageDimension) || (this->y_size > kMaxImageDimension)) {
		throw std::runtime_error("the reported frame size is unreasonably large");
	}
	
	if (this->total_number_of_images > kMaxNFrames) {
		throw std::runtime_error("the reported number of frames is unreasonably large");
	}
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
	this->checkForReasonableValues();
}

boost::shared_ptr<Eigen::MatrixXd> ImageLoaderSPE::readImage(const size_t index) {
	if (index >= total_number_of_images)
		throw IMAGE_INDEX_BEYOND_N_IMAGES(std::string("Requested more images than there are in the file"));
	
	uint64_t offset;
	uint32_t *currentUint32t = 0;
	float *currentFloat = 0;
	int16_t *currentInt16t = 0;
	uint16_t *currentUint16t = 0;
	boost::shared_ptr<Eigen::MatrixXd> image;
	
	uint64_t n_bytes_in_single_image;
	
	// determine how big we have to make the single image buffer and the offset
	switch(storage_type) {
		case STORAGE_TYPE_FP32:	// 4 byte float
		case STORAGE_TYPE_UINT32:	// 4-byte long
			n_bytes_in_single_image = x_size * y_size * 4;
			offset = header_length + index * (x_size) * (y_size) * 4;
			break;
		case STORAGE_TYPE_INT16:	// 2 byte signed short
		case STORAGE_TYPE_UINT16:	// 2 byte unsigned short
			n_bytes_in_single_image = x_size * y_size * 2;
			offset = header_length + index * (x_size) * (y_size) * 2;
			break;
		default:
			std::string error("Unable to determine the storage type used in ");
			error += this->filePath;
			throw CANNOT_DETERMINE_SPE_STORAGE_TYPE(error);
			break;
	}
	
	boost::scoped_array<char> single_image_buffer(new char[n_bytes_in_single_image]);
	image = boost::shared_ptr<Eigen::MatrixXd>(new Eigen::MatrixXd((int)x_size, (int)y_size));
	
	{
		boost::lock_guard<boost::mutex> locker(loadImagesMutex);
		file.seekg(offset);
		file.read((char *)single_image_buffer.get(), n_bytes_in_single_image);
		if (file.fail() != 0) {
			std::string error;
			error = "Error trying to read image data from \"";
			error += this->filePath;
			error += "\" assuming the SPE format";
			throw ERROR_READING_FILE_DATA(error);
		}
	}
	
	switch(storage_type) {
		case STORAGE_TYPE_FP32:	// 4-byte float
			currentFloat = (float *)single_image_buffer.get();
			for (size_t j  = 0; j < y_size; j++) {
				for (size_t i = 0; i < x_size; i++) {
					(*image)(i, j) = *currentFloat;
					++currentFloat;
				}
			}
			break;
			
		case STORAGE_TYPE_UINT32:	// 4-byte long
			currentUint32t = (uint32_t *)single_image_buffer.get();
			for (size_t j  = 0; j < y_size; j++) {
				for (size_t i = 0; i < x_size; i++) {
					(*image)(i, j) = *currentUint32t;
					++currentUint32t;
				}
			}
			break;
			
		case STORAGE_TYPE_INT16:	// 2-byte signed short
			currentInt16t = (int16_t *) single_image_buffer.get();
			for (size_t j  = 0; j < y_size; j++) {
				for (size_t i = 0; i < x_size; i++) {
					(*image)(i, j) = *currentInt16t;
					++currentInt16t;				
				}
			}
			break;
			
		case STORAGE_TYPE_UINT16: // 2-byte unsigned short
			currentUint16t = (uint16_t *)single_image_buffer.get();
			for (size_t j  = 0; j < y_size; j++) {
				for (size_t i = 0; i < x_size; i++) {
					(*image)(i, j) = *currentUint16t;
					++currentUint16t;
				}
			}
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
	frameYSize = wFrameYSize; frameXSize = wFrameXSize; this->total_number_of_images = WNImages;
#else
	result = sscanf(singleLine.c_str(), "Pixel number%zu 1 %zu %zu 1 %zu", &temp, &frameYSize, &frameXSize, &this->total_number_of_images);
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
	
	this->checkForReasonableValues();
}

boost::shared_ptr<Eigen::MatrixXd> ImageLoaderAndor::readImage(const size_t index) {
	if (index >= total_number_of_images)
		throw IMAGE_INDEX_BEYOND_N_IMAGES(std::string("Requested more images than there are in the file"));
	
	uint64_t offset;	// off_t is the size of the file pointer used by the OS
	float current_float = 0;
	uint64_t cache_offset = 0;
	
	boost::scoped_array<float> single_image_buffer(new float[x_size * y_size]);
	boost::shared_ptr<Eigen::MatrixXd> image (new Eigen::MatrixXd((int)x_size, (int)y_size));
	
	offset = header_length + index * (x_size) * (y_size) * sizeof(float);
	
	{
		boost::lock_guard<boost::mutex> locker(loadImagesMutex);
		file.seekg(offset);
		file.read((char *)single_image_buffer.get(), (x_size * y_size * sizeof(float)));
		if (file.fail() != 0) {
			std::string error;
			error = "Error trying to read image data from \"";
			error += this->filePath;
			error += "\" assuming the Andor format";
			throw ERROR_READING_FILE_DATA(error);
		}
	}
	
	
	// this is currently only safe on little-endian systems!
	for (size_t j  = 0; j < y_size; j++) {
		for (size_t i = 0; i < x_size; i++) {
			current_float = single_image_buffer[cache_offset];
			(*image)(i, j) = (double)current_float;
			cache_offset++;
		}
	}
	
	return image;
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
	
	ImageLoaderHamamatsu_HeaderStructure header;
	
	this->storage_type = STORAGE_TYPE_UINT16;
	
	file.seekg(0);
	
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
	
	if (header.storageFormat != 2) {	// not UINT16
		std::string error;
		error = "The file at \"";
		error += this->filePath;
		error += "\" specifies that it doesn't use UINT16 for storage. This usually means that the manufacturer's software corrupted the file.";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	this->header_length = header.commentLength + 64;
	this->x_size = header.xSize;
	this->y_size = header.ySize;
	this->total_number_of_images = header.nImages;
	
	// was there an error reading the file?
	if (file.fail() != 0) {
		std::string error;
		error = "Error parsing the header information in \"";
		error += this->filePath;
		error += "\" assuming the Hamamatsu format";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	file.seekg(0);
	this->checkForReasonableValues();
}


boost::shared_ptr<Eigen::MatrixXd> ImageLoaderHamamatsu::readImage(const size_t index) {
	if (index >= total_number_of_images)
		throw IMAGE_INDEX_BEYOND_N_IMAGES(std::string("Requested more images than there are in the file"));
	
	uint64_t offset;	// off_t is the size of the file pointer used by the OS
	size_t n_bytes_per_image = x_size * y_size * 2;
	offset = (index + 1) * header_length + index * (x_size) * (y_size) * 2;	// assume a 16-bit format
	
	boost::scoped_array<char> single_image_buffer(new char[n_bytes_per_image]);
	boost::shared_ptr<Eigen::MatrixXd> image (new Eigen::MatrixXd((int)x_size, (int)y_size));
	
	{
		boost::lock_guard<boost::mutex> locker(loadImagesMutex);
		file.seekg(offset);
		file.read(single_image_buffer.get(), n_bytes_per_image);
		if (file.fail() != 0) {
			std::string error;
			error = "Error reading image data from \"";
			error += this->filePath;
			error += "\" assuming the Hamamatsu format";
			throw ERROR_READING_FILE_DATA(error);
		}
	}
	
	uint16_t *uint16tPtr = (uint16_t *)single_image_buffer.get();
	
	for (size_t j  = 0; j < y_size; j++) {
		for (size_t i = 0; i < x_size; i++) {
			(*image)(i, j) = (double)(*uint16tPtr);
			++uint16tPtr;
		}
	}
	
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
	
	this->total_number_of_images = header.nImages;
	this->x_size = header.xSize;
	this->y_size = header.ySize;
	this->storage_type = header.storageFormat;
	this->header_length = sizeof(header);
	
	this->checkForReasonableValues();
}

boost::shared_ptr<Eigen::MatrixXd> ImageLoaderPDE::readImage(const size_t index) {
	if (index >= total_number_of_images)
		throw IMAGE_INDEX_BEYOND_N_IMAGES(std::string("Requested more images than there are in the file"));
	
	boost::shared_ptr<Eigen::MatrixXd> image;
	size_t n_pixels = this->x_size * this->y_size;
	size_t offset, imageSize;
	
	image = boost::shared_ptr<Eigen::MatrixXd>(new Eigen::MatrixXd((int)this->x_size, (int)this->y_size));
	
	switch (this->storage_type) {
		case STORAGE_TYPE_UINT16:
			imageSize = n_pixels * 2;
			offset = header_length + index * n_pixels * sizeof(uint16_t);
			break;
		case STORAGE_TYPE_UINT32:
		case STORAGE_TYPE_FP32:
			imageSize = n_pixels * 4;
			offset = header_length + index * n_pixels * sizeof(uint32_t);
			break;
		case STORAGE_TYPE_FP64:
			imageSize = n_pixels * 8;
			offset = header_length + index * n_pixels * sizeof(double);
			break;
		default:
			throw std::runtime_error("The data file does not appear to contain a recognized storage type");
			break;
	}
	
	boost::scoped_array<char> buffer(new char[imageSize]);
	
	{
		boost::lock_guard<boost::mutex> locker(loadImagesMutex);
		file.seekg(offset);
		file.read(buffer.get(), imageSize);
		if (file.fail() != 0) {
			std::string error;
			error = "Error reading image data from \"";
			error += this->filePath;
			error += "\" assuming the simple image format";
			throw ERROR_READING_FILE_DATA(error);
		}
	}
	
	
	switch (this->storage_type) {
		case STORAGE_TYPE_UINT16:
		{
			uint16_t *currentUint16t = (uint16_t *)buffer.get();
			for (size_t i = 0; i < this->x_size; ++i) {
				for (size_t j = 0; j < this->y_size; ++j) {
					(*image)(i, j) = *currentUint16t;
					++currentUint16t;
				}
			}
			break;
		}
		case STORAGE_TYPE_UINT32:
		{
			uint32_t *currentUint32t = (uint32_t *)buffer.get();
			for (size_t i = 0; i < this->x_size; ++i) {
				for (size_t j = 0; j < this->y_size; ++j) {
					(*image)(i, j) = *currentUint32t;
					++currentUint32t;
				}
			}
			break;
		}
		case STORAGE_TYPE_FP32:
		{
			float *currentFloat = (float *)buffer.get();
			for (size_t i = 0; i < this->x_size; ++i) {
				for (size_t j = 0; j < this->y_size; ++j) {
					(*image)(i, j) = *currentFloat;
					++currentFloat;
				}
			}
			break;
		}
		case STORAGE_TYPE_FP64:
		{
			double *currentDouble = (double *)buffer.get();
			for (size_t i = 0; i < this->x_size; ++i) {
				for (size_t j = 0; j < this->y_size; ++j) {
					(*image)(i, j) = *currentDouble;
					++currentDouble;
				}
			}
			break;
		}
		default:
			throw std::runtime_error("The data file does appear to contain a recognized storage type");
			break;
	}
	
	return image;
}

ImageLoaderTIFF::ImageLoaderTIFF(std::string rhs) {
	TIFFSetWarningHandler(NULL);
	this->filePath = rhs;
	
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
	
	
	this->total_number_of_images = directoryIndices.size();
	
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

boost::shared_ptr<Eigen::MatrixXd> ImageLoaderTIFF::readImage(const size_t index) {
	if (index >= total_number_of_images)
		throw IMAGE_INDEX_BEYOND_N_IMAGES(std::string("Requested more images than there are in the file"));
	
	boost::shared_ptr<Eigen::MatrixXd> image;
	int result;
	
	boost::lock_guard<boost::mutex> lock(loadImagesMutex);
	
	boost::shared_ptr<void> single_scanline_buffer (_TIFFmalloc(TIFFScanlineSize(tiff_file)), _TIFFfree);
	if (single_scanline_buffer.get() == NULL) {
		throw std::bad_alloc();
	}
	
	image = boost::shared_ptr<Eigen::MatrixXd>(new Eigen::MatrixXd((int)x_size, (int)y_size));
	
	result = TIFFSetDirectory(tiff_file, directoryIndices[index]);
	if (result != 1) {
		std::string error;
		error = "Unable to set the directory for the image at\"";
		error += this->filePath;
		error += "\"";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	for (size_t j = 0; j < y_size; ++j) {
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
				for (size_t k = 0; k < x_size; k+=2) {
					(*image)(k, j) = (double)(0x0000000F & (*charPtr));
					(*image)(k + 1, j) = (double)((0x000000F0 & (*charPtr)) / 16);			  
					++charPtr;
				}
				break;
			}
			case STORAGE_TYPE_UINT8:
			{
				char *charPtr = (char *)single_scanline_buffer.get();
				for (size_t k = 0; k < x_size; ++k) {
					(*image)(k, j) = (double)(*charPtr);
					++charPtr;
				}
				break;
			}
			case STORAGE_TYPE_UINT16:
			{
				uint16_t *uint16tPtr = (uint16_t *)single_scanline_buffer.get();
				for (size_t k = 0; k < x_size; ++k) {
					(*image)(k, j) = (double)(*uint16tPtr);
					++uint16tPtr;
				}
				break;
			}
			case STORAGE_TYPE_UINT32:
			{
				uint32_t *uint32tPtr = (uint32_t *)single_scanline_buffer.get();
				for (size_t k = 0; k < x_size; ++k) {
					(*image)(k, j) = (double)(*uint32tPtr);
					++uint32tPtr;
				}
				break;
			}
			case STORAGE_TYPE_FP32:
			{
				float *floatPtr = (float *)single_scanline_buffer.get();
				for (size_t k = 0; k < x_size; ++k) {
					(*image)(k, j) = (double)(*floatPtr);
					++floatPtr;
				}
				break;
			}
			case STORAGE_TYPE_FP64:
			{
				double *doublePtr = (double *)single_scanline_buffer.get();
				for (size_t k = 0; k < x_size; ++k) {
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
	
	return image;
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
	
	this->checkForReasonableValues();
}

boost::shared_ptr<Eigen::MatrixXd> ImageLoaderIgor::readImage(const size_t index) {
	if (index >= total_number_of_images)
		throw IMAGE_INDEX_BEYOND_N_IMAGES(std::string("Requested more images than there are in the file"));
	
	double value[2];
	long indices[3];
	int result;
	
	boost::shared_ptr<Eigen::MatrixXd> image;
	
	// no mutex locking is required since these calls are all threadsafe
	
	indices[2] = index;
	image = boost::shared_ptr<Eigen::MatrixXd> (new Eigen::MatrixXd((int)x_size, (int)y_size));
	
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
	
	return image;
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
	
	if (overwrite == 0) {
		std::ifstream input_test;
		input_test.open(file_path.c_str(), std::ios::in | std::ios::binary);
		if (input_test.good() == 1) {
			std::string error("The output file at ");
			error += this->file_path;
			error += " already exists.";
			throw OUTPUT_FILE_ALREADY_EXISTS(error);	// escape without overwriting
		}
	}
	
	file.open(file_path.c_str(), std::ios::binary | std::ios::out | std::ios::trunc);	// DANGER: OVERWRITING THE FILE
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

void PDEImageOutputWriter::write_image(boost::shared_ptr<Eigen::MatrixXd> imageToWrite) {
	
	// determine the size of the frames
	size_t currentXSize = imageToWrite->rows();
	size_t currentYSize = imageToWrite->cols();
	size_t n_pixels = currentXSize * currentYSize;
	
	size_t offset = 0;
	
	if (this->n_images_written == 0) {
		this->x_size = currentXSize;
		this->y_size = currentYSize;
	} else {
		if ((currentXSize != this->x_size) || (currentYSize != this->y_size)) {
			throw std::runtime_error("Tried to write an image with different dimensions to an SPE file");
		}
	}
	
	switch (this->storageType) {
		case STORAGE_TYPE_UINT16:
		{
			boost::scoped_array<uint16_t> buffer(new uint16_t[n_pixels]);
			for (size_t i = 0; i < this->x_size; ++i) {
				for (size_t j = 0; j < this->y_size; ++j) {
					buffer[offset] = (uint16_t)(*imageToWrite)(i, j);
					++offset;
				}
			}
			this->file.write((char *)buffer.get(), n_pixels * sizeof(uint16_t));
			break;
		}
		case STORAGE_TYPE_UINT32:
		{
			boost::scoped_array<uint32_t> buffer(new uint32_t[n_pixels]);
			for (size_t i = 0; i < this->x_size; ++i) {
				for (size_t j = 0; j < this->y_size; ++j) {
					buffer[offset] = (uint32_t)(*imageToWrite)(i, j);
					++offset;
				}
			}
			this->file.write((char *)buffer.get(), n_pixels * sizeof(uint32_t));
			break;
		}
		case STORAGE_TYPE_FP32:
		{
			boost::scoped_array<float> buffer(new float[n_pixels]);
			for (size_t i = 0; i < this->x_size; ++i) {
				for (size_t j = 0; j < this->y_size; ++j) {
					buffer[offset] = (float)(*imageToWrite)(i, j);
					++offset;
				}
			}
			this->file.write((char *)buffer.get(), n_pixels * sizeof(float));
			break;
		}
		case STORAGE_TYPE_FP64:
		{
			boost::scoped_array<double> buffer(new double[n_pixels]);
			for (size_t i = 0; i < this->x_size; ++i) {
				for (size_t j = 0; j < this->y_size; ++j) {
					buffer[offset] = (double)(*imageToWrite)(i, j);
					++offset;
				}
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
	
	TIFFSetWarningHandler(NULL);
	this->file_path = rhs;
	this->compression = compression_rhs;
	this->storageType = storageType_rhs;
	
	if (overwrite == 0) {
		std::ifstream input_test;
		input_test.open(file_path.c_str(), std::ios::in | std::ios::binary);
		if (input_test.good() == 1) {
			std::string error("The output file at ");
			error += this->file_path;
			error += " already exists.";
			throw OUTPUT_FILE_ALREADY_EXISTS(error);	// escape without overwriting
		}
		input_test.close();
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



void TIFFImageOutputWriter::write_image(boost::shared_ptr<Eigen::MatrixXd> imageToWrite) {
	
	size_t x_size = imageToWrite->rows();
	size_t y_size = imageToWrite->cols();
	
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
	boost::scoped_array<char> scanLine(new char[x_size * bitsPerSample / 8]);
	
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
			{
				int8_t* charPtr = (int8_t *)scanLine.get();
				for (size_t i = 0; i < x_size; i++) {
					charPtr[i] = (int8_t)(*imageToWrite)(i, j);
				}
				break;
			}
			case STORAGE_TYPE_UINT8:
			{
				uint8_t* uCharPtr = (uint8_t *)scanLine.get();
				for (size_t i = 0; i < x_size; i++) {
					uCharPtr[i] = (uint8_t)(*imageToWrite)(i, j);
				}
				break;
			}
			case STORAGE_TYPE_INT16:
			{
				int16_t* int16Ptr = (int16_t *)scanLine.get();
				for (size_t i = 0; i < x_size; i++) {
					int16Ptr[i] = (int16_t)(*imageToWrite)(i, j);
				}
				break;
			}
			case STORAGE_TYPE_UINT16:
			{
				uint16_t* uInt16Ptr = (uint16_t *)scanLine.get();
				for (size_t i = 0; i < x_size; i++) {
					uInt16Ptr[i] = (uint16_t)(*imageToWrite)(i, j);
				}
				break;
			}
			case STORAGE_TYPE_INT32:
			{
				int32_t* int32Ptr = (int32_t *)scanLine.get();
				for (size_t i = 0; i < x_size; i++) {
					int32Ptr[i] = (int32_t)(*imageToWrite)(i, j);
				}
				break;
			}
			case STORAGE_TYPE_UINT32:
			{
				uint32_t* uInt32Ptr = (uint32_t *)scanLine.get();
				for (size_t i = 0; i < x_size; i++) {
					uInt32Ptr[i] = (uint32_t)(*imageToWrite)(i, j);
				}
				break;
			};
			case STORAGE_TYPE_INT64:
			{
				int64_t* int64Ptr = (int64_t *)scanLine.get();
				for (size_t i = 0; i < x_size; i++) {
					int64Ptr[i] = (int64_t)(*imageToWrite)(i, j);
				}
				break;
			}
			case STORAGE_TYPE_UINT64:
			{
				uint64_t* uInt64Ptr = (uint64_t *)scanLine.get();
				for (size_t i = 0; i < x_size; i++) {
					uInt64Ptr[i] = (uint64_t)(*imageToWrite)(i, j);
				}
				break;
			}
			case STORAGE_TYPE_FP32:
			{
				float* floatPtr = (float *)scanLine.get();
				for (size_t i = 0; i < x_size; i++) {
					floatPtr[i] = (float)(*imageToWrite)(i, j);
				}
				break;
			}
			case STORAGE_TYPE_FP64:
			{
				double* doublePtr = (double *)scanLine.get();
				for (size_t i = 0; i < x_size; i++) {
					doublePtr[i] = (double)(*imageToWrite)(i, j);
				}
				break;
			}
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
	
	this->overwrite = overwrite_rhs;
	this->outputWave = NULL;
	this->waveName = waveName_rhs;
	this->nImagesTotal = nImages_rhs;
	this->n_images_written = 0;
	this->storageType = storageType_rhs;
	
	// check if the output wave already exists
	if (this->overwrite != 1) {
		waveHndl wave = FetchWaveUsingFullPath(waveName_rhs);
		if (wave != NULL) {
			std::string error ("the output wave ");
			error += waveName_rhs;
			error += " already exists and was not overwritten (/O flag)";
		}
	}
}

void IgorImageOutputWriter::write_image(boost::shared_ptr<Eigen::MatrixXd> imageToWrite) {
	long indices[MAX_DIMENSIONS + 1];
	int result;
	double value[2];
	
	size_t x_size = imageToWrite->rows();
	size_t y_size = imageToWrite->cols();
	
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
	
	for (size_t j  = 0; j < y_size; j++) {
		for (size_t i = 0; i < x_size; i++) {
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
