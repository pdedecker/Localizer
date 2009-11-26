/*
 *  PALM_analysis_FileIO.cpp
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 30/03/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "PALM_analysis_FileIO.h"

XOPFileHandler::~XOPFileHandler() {
	if (fileRef != NULL) {
		XOPCloseFile(fileRef);
		fileRef = NULL;
	}
}

void XOPFileHandler::open(const char* path_rhs) {
	int err;
	
	err = XOPOpenFile(path_rhs, 0, &fileRef);
	if (err != 0) {
		string error;
		stringstream ss;
		ss << "Error " << err << " returned using open() on the image file at \"" << path_rhs << "\"\r";
		error = ss.str();
		throw ERROR_READING_FILE_DATA(error);
	}
	
	path = path_rhs;
}

void XOPFileHandler::close() {
	int err;
	if (fileRef != NULL) {
		err = XOPCloseFile(fileRef);
		if (err != 0) {
			string error;
			stringstream ss;
			ss << "Error " << err << " returned using close() on the image file at \"" << path << "\"\r";
			error = ss.str();
			throw ERROR_READING_FILE_DATA(error);
		}
		fileRef = NULL;
	}
}

void XOPFileHandler::get(char& c) {
	int err;
	err = XOPReadFile2(fileRef, 1, &c, NULL);
	if (err != 0) {
		string error;
		stringstream ss;
		ss << "Error " << err << " returned using get() on the image file at \"" << path << "\"\r";
		error = ss.str();
		throw ERROR_READING_FILE_DATA(error);
	}
}

void XOPFileHandler::read(char *buffer, size_t nBytes) {
	int err;
	err = XOPReadFile2(fileRef, nBytes, buffer, NULL);
	if (err != 0) {
		string error;
		stringstream ss;
		ss << "Error " << err << " returned using read() on the image file at \"" << path << "\"\r";
		error = ss.str();
		throw ERROR_READING_FILE_DATA(error);
	}
}

void XOPFileHandler::getline(char *buffer, size_t nMax) {
	int err;
	err = XOPReadLine(fileRef, buffer, nMax, NULL);
	if (err != 0) {
		string error;
		stringstream ss;
		ss << "Error " << err << " returned using getline() on the image file at \"" << path << "\"\r";
		error = ss.str();
		throw ERROR_READING_FILE_DATA(error);
	}
}

uint64_t XOPFileHandler::tellg() {
	double dPosition;
	uint64_t pos;
	int err;
	
	err = XOPGetFilePosition2(fileRef, &dPosition);
	if (err != 0) {
		string error;
		stringstream ss;
		ss << "Error " << err << " returned using tellg() on the image file at \"" << path << "\"\r";
		error = ss.str();
		throw ERROR_READING_FILE_DATA(error);
	}
	
	pos = (uint64_t) dPosition;
	return pos;
}

void XOPFileHandler::seekg(uint64_t pos) {
	double dPosition;
	
	int err;
	
	dPosition = (double)pos;
	
	err = XOPSetFilePosition2(fileRef, dPosition);
	if (err != 0) {
		string error;
		stringstream ss;
		ss << "Error " << err << " returned using seekg() on the image file at \"" << path << "\"\r";
		error = ss.str();
		throw ERROR_READING_FILE_DATA(error);
	}
}

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


ImageLoader::ImageLoader() {
	path.assign("");
	total_number_of_images = 0;
	x_size = 0;
	y_size = 0;
	header_length = 0;
}

ImageLoader::~ImageLoader() {
	if (file.is_open() == 1)
		file.close();
}

boost::shared_ptr<PALMMatrix<double> > ImageLoader::get_nth_image(const size_t n) {
	vector<boost::shared_ptr<PALMMatrix <double> > > images;
	if (n >= total_number_of_images)
		throw IMAGE_INDEX_BEYOND_N_IMAGES();
	
	size_t firstImageToLoad = n;
	size_t lastImageToLoad = n;
	
	images = ReadImagesFromDisk(firstImageToLoad, lastImageToLoad);
	
	return images.at(0);
}

ImageLoaderSPE::ImageLoaderSPE(string rhs, size_t image_cache_size_rhs) {
	path = rhs;
	
	header_length = 4100;
	
	file.open(path.c_str(), ios::binary | ios::in);
	if (file.fail() == 1) {
		throw CANNOT_OPEN_FILE();
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
			throw CANNOT_DETERMINE_SPE_STORAGE_TYPE();
			break;
	}
	current_bytes = 0;
	
	// was there an error sometime during this procedure that would have caused the reading to fail?
	if (file.fail() != 0) {
		string error;
		error = "Error parsing the header information in \"";
		error += path;
		error += "\" assuming the SPE format\r";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	file.seekg(0);
}

vector<boost::shared_ptr<PALMMatrix <double> > > ImageLoaderSPE::ReadImagesFromDisk(size_t const nStart, size_t const nEnd) {	
	uint64_t offset;
	long current_long = 0;
	float current_float = 0;
	short current_short = 0;
	unsigned short current_unsigned_short = 0;
	string error;
	boost::shared_ptr<PALMMatrix<double> > image;
	vector<boost::shared_ptr<PALMMatrix <double> > > requestedImages;
	
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
			throw CANNOT_DETERMINE_SPE_STORAGE_TYPE();
			break;
	}
	
	boost::scoped_array<float> single_image_buffer_float(new float[x_size * y_size]);
	boost::scoped_array<char> single_image_buffer(new char[n_bytes_in_single_image]);
	
	for (uint64_t i = nStart; i <= nEnd; i++) {
		
		image = boost::shared_ptr<PALMMatrix<double> >(new PALMMatrix<double>(x_size, y_size));
		
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
				throw CANNOT_DETERMINE_SPE_STORAGE_TYPE();
				break;
		}
		
		
		file.seekg(offset);
		cache_offset = 0;
		
		file.read((char *)single_image_buffer.get(), n_bytes_in_single_image);
		if (file.fail() != 0) {
			string error;
			error = "Error trying to read image data from \"";
			error += path;
			error += "\" assuming the SPE format\r";
			loadImagesMutex.unlock();
			throw ERROR_READING_FILE_DATA(error);
		}
		
		switch(storage_type) {
			case STORAGE_TYPE_FP32:	// 4-byte float
				// this is currently only safe on little-endian systems!
				for (size_t j  = 0; j < y_size; j++) {
					for (size_t i = 0; i < x_size; i++) {
						current_float = single_image_buffer_float[cache_offset];
						
						image->set(i, j, (double)current_float);
						
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
						
						image->set(i, j, (double)current_long);
						
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
						
						image->set(i, j, (double)current_short);
						
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
						
						image->set(i, j, (double)current_unsigned_short);
						
						cache_offset += 2;
						
						//	current_unsigned_short = 0;
					}
				}
				break;
				
			default:
				loadImagesMutex.unlock();
				throw CANNOT_DETERMINE_SPE_STORAGE_TYPE();
				break;
				
		}
		
		requestedImages.push_back(image);
	}
	
	file.seekg(0);
	
	loadImagesMutex.unlock();
	
	return requestedImages;
}

ImageLoaderAndor::ImageLoaderAndor(string rhs, size_t image_cache_size_rhs) {
	path = rhs;
	
	file.open(path.c_str(), ios::binary | ios::in);
	if (file.fail() == 1) {
		throw CANNOT_OPEN_FILE();
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
	boost::scoped_array<char> headerBuffer(new char[60001]);
	
	storage_type = STORAGE_TYPE_FP32;
	
	file.read(headerBuffer.get(), 60000);
	headerBuffer[60000] = '\0';	// NULL-terminate the string
	
	for (int i = 0; i < 1000; ++i) {	// make sure that there are not intermediate NULL characters
		if (headerBuffer[i] == '\0') {
			headerBuffer[i] = '1';
		}
	}
	
	string headerString(headerBuffer.get());
	stringstream ss(headerString, ios::in);
	
	// the important information on the measurement is on lines 22 and 23 (numbered from 1)
	for (int i = 0; i < 21; ++i) {
		ss.getline(headerBuffer.get(), 1024);
	}
	
	// get the first information
	temp = ss.tellg();
	ss.seekg(temp + 12);	// skip the "Pixel number"
	ss >> temp;	// extracts "65538"
	ss >> temp;	// extracts "1"
	ss >> frameYSize;
	ss >> frameXSize;
	ss >> temp;	// extracts "1"
	ss >> total_number_of_images;
	ss.getline(headerBuffer.get(), 1024);	// discard the rest of the line
	
	ss >> temp;	// extracts "65538"
	ss >> xStart;
	ss >> yEnd;
	ss >> xEnd;
	ss >> yStart;
	ss >> xBinning;
	ss >> yBinning;
	ss.getline(headerBuffer.get(), 1024);	// discard the rest of the line
	
	x_size = (xEnd - xStart + 1) / xBinning;
	y_size = (yEnd - yStart + 1) / yBinning;	// integer division
	
	// now there are some lines that may contain timestamps. There are as many lines as there are images
	for (size_t i = 0;  i < total_number_of_images; i++) {
		ss.getline(headerBuffer.get(), 1024);
	}
	
	header_length = ss.tellg();
	
	// did some error happen while reading the file?
	if (file.fail() != 0) {
		string error;
		error = "Error parsing the header information in \"";
		error += path;
		error += "\" assuming the Andor format\r";
		throw ERROR_READING_FILE_DATA(error);
	}
}

vector<boost::shared_ptr<PALMMatrix <double> > > ImageLoaderAndor::ReadImagesFromDisk(size_t const nStart, size_t const nEnd) {	
	uint64_t offset;	// off_t is the size of the file pointer used by the OS
	float current_float = 0;
	
	boost::shared_ptr<PALMMatrix<double> > image;
	vector<boost::shared_ptr<PALMMatrix <double> > > requestedImages;
	
	loadImagesMutex.lock();
	
	boost::scoped_array<float> single_image_buffer(new float[x_size * y_size]);
	uint64_t cache_offset;
	
	for (uint64_t i = nStart; i <= nEnd; i++) {
		image = boost::shared_ptr<PALMMatrix<double> >(new PALMMatrix<double>(x_size, y_size));
		
		offset = header_length + i * (x_size) * (y_size) * sizeof(float);
		file.seekg(offset);
		
		cache_offset = 0;
		
		file.read((char *)single_image_buffer.get(), (x_size * y_size * sizeof(float)));
		if (file.fail() != 0) {
			string error;
			error = "Error trying to read image data from \"";
			error += path;
			error += "\" assuming the Andor format\r";
			loadImagesMutex.unlock();
			throw ERROR_READING_FILE_DATA(error);
		}
		
		
		// this is currently only safe on little-endian systems!
		for (size_t j  = 0; j < y_size; j++) {
			for (size_t i = 0; i < x_size; i++) {
				current_float = single_image_buffer[cache_offset];
				image->set(i, j, (double)current_float);
				cache_offset++;
			}
		}
		
		requestedImages.push_back(image);
	}
	
	loadImagesMutex.unlock();
	
	return requestedImages;
}

ImageLoaderHamamatsu::ImageLoaderHamamatsu(string rhs, size_t image_cache_size_rhs) {
	path = rhs;
	
	file.open(path.c_str(), ios::binary | ios::in);
	if (file.fail() == 1) {
		throw CANNOT_OPEN_FILE();
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
		string error;
		error = "The file at \"";
		error += path;
		error += "\" specifies that it doesn't use UINT16 for storage. Please ask Peter for help.\r";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	header_length = header.commentLength + 64;
	x_size = header.xSize;
	y_size = header.ySize;
	total_number_of_images = header.nImages;
	
	// was there an error reading the file?
	if (file.fail() != 0) {
		string error;
		error = "Error parsing the header information in \"";
		error += path;
		error += "\" assuming the Hamamatsu format\r";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	file.seekg(0);
}


vector<boost::shared_ptr<PALMMatrix <double> > > ImageLoaderHamamatsu::ReadImagesFromDisk(size_t const nStart, size_t const nEnd) {	
	uint64_t offset;	// off_t is the size of the file pointer used by the OS
	
	boost::shared_ptr<PALMMatrix<double> > image;
	vector<boost::shared_ptr<PALMMatrix <double> > > requestedImages;
	
	loadImagesMutex.lock();
	
	unsigned int current_uint;
	size_t n_bytes_per_image = x_size * y_size * 2;
	boost::scoped_array<char> single_image_buffer(new char[n_bytes_per_image]);
	uint64_t cache_offset;
	
	for (uint64_t i = nStart; i <= nEnd; i++) {
		image = boost::shared_ptr<PALMMatrix<double> >(new PALMMatrix<double>(x_size, y_size));
		
		offset = (i + 1) * header_length + i * (x_size) * (y_size) * 2;	// assume a 16-bit format
		file.seekg(offset);
		
		file.read(single_image_buffer.get(), n_bytes_per_image);
		if (file.fail() != 0) {
			string error;
			error = "Error reading image data from \"";
			error += path;
			error += "\" assuming the Hamamatsu format\r";
			loadImagesMutex.unlock();
			throw ERROR_READING_FILE_DATA(error);
		}
		
		cache_offset = 0;
		
		for (size_t j  = 0; j < y_size; j++) {
			for (size_t i = 0; i < x_size; i++) {
				current_uint = 0x000000FF & single_image_buffer[cache_offset + 1];	// little endian
				current_uint *= 256;
				current_uint = current_uint | (0x000000FF & single_image_buffer[cache_offset]);
				
				image->set(i, j, (double)current_uint);
				
				cache_offset += 2;	// TWO bytes per value (UINT16)
			}
		}
		
		
		requestedImages.push_back(image);
	}
	
	file.seekg(0);
	
	loadImagesMutex.unlock();
	
	return requestedImages;
}

SimpleImageLoader::SimpleImageLoader(string rhs, size_t image_cache_size_rhs) {
	path = rhs;
	
	file.open(path.c_str(), ios::binary | ios::in);
	if (file.fail() == 1) {
		throw CANNOT_OPEN_FILE();
	}
	
	parse_header_information();
}

SimpleImageLoader::~SimpleImageLoader() {
	if (file.is_open() == 1) {
		file.close();
	}
}

vector<boost::shared_ptr<PALMMatrix <double> > > SimpleImageLoader::ReadImagesFromDisk(size_t const nStart, size_t const nEnd) {
	uint64_t offset;
	size_t array_offset;
	boost::shared_ptr<PALMMatrix<double> > image;
	vector<boost::shared_ptr<PALMMatrix <double> > > requestedImages;
	
	loadImagesMutex.lock();
	
	boost::scoped_array<float> single_image_buffer(new float[x_size * y_size]);
	
	for (uint64_t i = nStart; i <= nEnd; i++) {
		image = boost::shared_ptr<PALMMatrix<double> >(new PALMMatrix<double>(x_size, y_size));
		
		offset = 3 * sizeof(size_t) + i * (x_size) * (y_size) * sizeof(float);
		file.seekg(offset);
		
		array_offset = 0;
		
		// this is currently only safe on little-endian systems!
		
		file.read((char *)single_image_buffer.get(), (x_size * y_size * sizeof(float)));
		if (file.fail() != 0) {
			string error;
			error = "Error reading image data from \"";
			error += path;
			error += "\" assuming the simple image format\r";
			loadImagesMutex.unlock();
			throw ERROR_READING_FILE_DATA(error);
		}
		
		for (size_t j  = 0; j < y_size; j++) {
			for (size_t i = 0; i < x_size; i++) {
				image->set(i, j, single_image_buffer[array_offset]);
				array_offset++;
			}
		}
		
		requestedImages.push_back(image);
	}
	
	file.seekg(0);
	
	loadImagesMutex.unlock();
	
	return requestedImages;
}


void SimpleImageLoader::parse_header_information() {
	
	storage_type = STORAGE_TYPE_FP64;
	
	file.seekg(0);
	
	// the first 12 bytes of the file contain this information
	file.read((char *)&x_size, sizeof(size_t));
	file.read((char *)&y_size, sizeof(size_t));
	file.read((char *)&total_number_of_images, sizeof(size_t));
	
	if (file.fail() != 0) {
		string error;
		error = "Error parsing the header information in \"";
		error += path;
		error += "\" assuming the simple image format\r";
		throw ERROR_READING_FILE_DATA(error);
	}
}

ImageLoaderTIFF::ImageLoaderTIFF(string rhs, size_t image_cache_size_rhs) {
	path = rhs;
	
	tiff_file = NULL;
	
	tiff_file = TIFFOpen(path.c_str(), "r");
	if (tiff_file == NULL) {
		throw CANNOT_OPEN_FILE();
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
		string error;
		error = "The image at\"";
		error += path;
		error += "\" is not a grayscale image\r";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	if ((result_uint16 != 0) && (result_uint16 != 1)) {	// not a grayscale image
		string error;
		error = "The image at\"";
		error += path;
		error += "\" is not a grayscale image\r";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	// is it a binary image?
	result = TIFFGetField(tiff_file, TIFFTAG_BITSPERSAMPLE, &result_uint16);
	if (result != 1) {
		string error;
		error = "The image at\"";
		error += path;
		error += "\" is not a grayscale image\r";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	if (result_uint16 < 4) {	// 4 is the minimum number of bits allowed for grayscale images in the tiff specification, so this is a bilevel image
		string error;
		error = "The image at\"";
		error += path;
		error += "\" is not a grayscale image\r";
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
			string error;
			error = "The SampleFormat of the image at\"";
			error += path;
			error += "\" is unknown\r";
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
				string error;
				error = "The SampleFormat of the image at\"";
				error += path;
				error += "\" is unknown\r";
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
				string error;
				error = "The SampleFormat of the image at\"";
				error += path;
				error += "\" is unknown\r";
				throw ERROR_READING_FILE_DATA(error);
				break;
		}
	}
	
	// what is the x size?
	result = TIFFGetField(tiff_file, TIFFTAG_IMAGEWIDTH, &result_uint32);
	if (result != 1) {
		string error;
		error = "The image at\"";
		error += path;
		error += "\" does not specify a width\r";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	x_size = (size_t)(result_uint32);
	
	// what is the y size?
	result = TIFFGetField(tiff_file, TIFFTAG_IMAGELENGTH, &result_uint32);
	if (result != 1) {
		string error;
		error = "The image at\"";
		error += path;
		error += "\" does not specify a height\r";
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
		string error;
		error = "Unable to set the directory to '0' for the image at\"";
		error += path;
		error += "\"\r";
		throw ERROR_READING_FILE_DATA(error);
	}
}

vector<boost::shared_ptr<PALMMatrix <double> > > ImageLoaderTIFF::ReadImagesFromDisk(size_t const nStart, size_t const nEnd) {
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
	boost::shared_ptr<PALMMatrix<double> > image;
	vector<boost::shared_ptr<PALMMatrix <double> > > requestedImages;
	int result;
	
	loadImagesMutex.lock();
	
	single_scanline_buffer = (char *)_TIFFmalloc(TIFFScanlineSize(tiff_file));
	if (single_scanline_buffer == NULL) {
		string error;
		error = "unable to allocate scanline buffer in ImageLoaderTIFF::get_nth_image()\r";
		throw OUT_OF_MEMORY(error);
	}
	
	for (size_t i = nStart; i <= nEnd; i++) {
		try {
			image = boost::shared_ptr<PALMMatrix<double> >(new PALMMatrix<double>(x_size, y_size));
		}
		catch (std::bad_alloc) {
			_TIFFfree(single_scanline_buffer);
			string error;
			error = "unable to allocate new_image in ImageLoaderTIFF::get_nth_image()\r";
			loadImagesMutex.unlock();
			throw OUT_OF_MEMORY(error);
		}
		
		result = TIFFSetDirectory(tiff_file, directoryIndices.at(i));
		if (result != 1) {
			_TIFFfree(single_scanline_buffer);
			string error;
			error = "Unable to set the directory to '0' for the image at\"";
			error += path;
			error += "\"\r";
			loadImagesMutex.unlock();
			throw ERROR_READING_FILE_DATA(error);
		}
		
		for (size_t j = 0; j < y_size; ++j) {
			result = TIFFReadScanline(tiff_file, single_scanline_buffer, j, 0);	// sample is ignored
			if (result != 1) {
				_TIFFfree(single_scanline_buffer);
				string error;
				error = "Unable to read a scanline from the image at\"";
				error += path;
				error += "\"\r";
				loadImagesMutex.unlock();
				throw ERROR_READING_FILE_DATA(error);
			}
			
			switch (storage_type) {	// handle the different possibilities (floating, integer) and variable sizes
				case STORAGE_TYPE_UINT4:
					scanline_pointer = single_scanline_buffer;
					for (size_t k = 0; k < x_size; ++k) {
						if ((k % 2) == 0) {	// this is an even pixel, we use only the first 4 bits
							current_uint16 = 0x0000000F & (*scanline_pointer);
							image->set(k, j, (double)current_uint16);
						} else {	// this is an odd pixel, use the last 4 bits and increment the scanline_pointer
							current_uint16 = 0x000000F0 & (*scanline_pointer);
							image->set(k, j, (double)current_uint16);
							scanline_pointer += 1;
						}
					}
					break;
					
				case STORAGE_TYPE_UINT8:
					scanline_pointer = single_scanline_buffer;
					for (size_t k = 0; k < x_size; ++k) {
						current_uint16 = (uint16_t)(*scanline_pointer);
						image->set(k, j, (double)current_uint16);
						scanline_pointer += 1;
					}
					break;
					
				case STORAGE_TYPE_UINT16:
					uint16Ptr = (uint16_t*)single_scanline_buffer;
					for (size_t k = 0; k < x_size; ++k) {
						current_uint16 = (*uint16Ptr);
						image->set(k, j, (double)current_uint16);
						uint16Ptr += 1;
					}
					break;
					
				case STORAGE_TYPE_UINT32:
					uint32Ptr = (uint32_t*)single_scanline_buffer;
					for (size_t k = 0; k < x_size; ++k) {
						current_uint32 = (*uint32Ptr);
						image->set(k, j, (double)current_uint32);
						uint32Ptr += 1;
					}
					break;
					
				case STORAGE_TYPE_FP32:
					floatPtr = (float *)single_scanline_buffer;
					for (size_t k = 0; k < x_size; ++k) {
						current_float = *floatPtr;
						image->set(k, j, (double)current_float);
						floatPtr += 1;
					}
					break;
					
				case STORAGE_TYPE_FP64:
					doublePtr = (double *)single_scanline_buffer;
					for (size_t k = 0; k < x_size; ++k) {
						current_double = *doublePtr;
						image->set(k, j, current_double);
						floatPtr += 1;
					}
					break;
					
				default:
					_TIFFfree(single_scanline_buffer);
					string error;
					error = "Invalid floating point data size for the image at\"";
					error += path;
					error += "\"\r";
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
		string error;
		error = "Invalid to set the directory to '0' for the image at\"";
		error += path;
		error += "\"\r";
		loadImagesMutex.unlock();
		throw ERROR_READING_FILE_DATA(error);
	}
	
	loadImagesMutex.unlock();
	
	return requestedImages;
}

ImageLoaderIgor::ImageLoaderIgor(string waveName, size_t image_cache_size_rhs) {
	int err;
	size_t waveNameOffset;
	string dataFolderPath;
	DataFolderHandle dataFolder;
	size_t position;
	int waveType;
	// try to get images from an Igor wave
	
	// did we receive a wavename that does not reference any data folders? (no ':' in the name)
	// in that case we assume that the wave is located in the current data folder
	if (waveName.find(':') == string::npos) {
		igor_data_wave = FetchWaveFromDataFolder(NULL, waveName.c_str());
		if (igor_data_wave == NULL) {
			throw NOWAV;
		}
	} else {
		// we found a ':' in the string
		// therefore we're being passed a wave location that also includes datafolder information
		// we need to parse this string to get out the data folder handle and the wave handle
		waveNameOffset = waveName.rfind(':');
		if (waveNameOffset != 0) {	// check for the case where the user specifies something like ":waveName"
			dataFolderPath = waveName.substr(0, waveNameOffset);
		} else {
			dataFolderPath = ":";
		}
		err = GetNamedDataFolder(NULL, dataFolderPath.c_str(), &dataFolder);
		if (err != 0)
			throw err;
		
		// retain only the waveName part, discard the information about the data folder
		waveName = waveName.substr(waveNameOffset + 1, waveName.length() - waveNameOffset - 1);
		
		// if the waveName part is quoted ('') then Igor chokes on this
		// so we remove the quotes
		while ((position = waveName.find('\'')) != string::npos) {
			waveName.erase(position, 1);
		}
		
		igor_data_wave = FetchWaveFromDataFolder(dataFolder, waveName.c_str());
		if (igor_data_wave == NULL) {
			throw NOWAV;
		}
	}
	
	
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

vector<boost::shared_ptr<PALMMatrix <double> > > ImageLoaderIgor::ReadImagesFromDisk(size_t const nStart, size_t const nEnd) {
	double value[2];
	long indices[3];
	int result;
	
	boost::shared_ptr<PALMMatrix<double> > image;
	vector<boost::shared_ptr<PALMMatrix <double> > > requestedImages;
	
	// no mutex locking is required since these calls are all threadsafe
	
	for (size_t n = nStart; n <= nEnd; ++n) {
		indices[2] = n;
		image = boost::shared_ptr<PALMMatrix<double> > (new PALMMatrix<double>(x_size, y_size));
		
		for (size_t j = 0; j < y_size; j++) {
			for (size_t i = 0; i < x_size; i++) {
				indices[0] = (long)i;
				indices[1] = (long)j;
				
				result = MDGetNumericWavePointValue(igor_data_wave, indices, value);
				if (result != 0) {
					throw result;
				}
				
				image->set(i, j, value[0]);
			}
		}
		
		requestedImages.push_back(image);
		
	}
	
	return requestedImages;
}



ImageOutputWriter::ImageOutputWriter() {
	file_path.assign("");
	n_images_written = 0;
}

ImageOutputWriter::ImageOutputWriter(const string &rhs, int overwrite) {
	file_path.assign("");
	n_images_written = 0;
}



SimpleImageOutputWriter::SimpleImageOutputWriter(const string &rhs,int overwrite) {
	// if overwrite is non-zero then we overwrite any file that exists at the output path
	// if it is set to zero then we throw an error and abort instead of overwriting
	file_path = rhs;
	int header_length = 3 * sizeof(size_t);
	
	char *header = new char[header_length];
	if (header == NULL) {
		string error;
		error = "unable to allocate header in SimpleImageOutputWriter::SimpleImageOutputWriter()\r";
		throw OUT_OF_MEMORY(error);
	}
	for (int i = 0; i < header_length; i++) {
		header[i] = 0;
	}
	
	if (overwrite == 0) {
		ifstream input_test;
		input_test.open(file_path.c_str(), ios::in | ios::binary);
		input_test.close();
		if (input_test.fail() == 0) {
			throw OUTPUT_FILE_ALREADY_EXISTS();	// escape without overwriting
		} else {
			file.open(file_path.c_str(), ios::binary | ios::out);
		}
	}
	
	if (overwrite != 0) {
		file.open(file_path.c_str(), ios::binary | ios::out | ios::trunc);	// DANGER: OVERWRITING THE FILE
	}
	
	if (file.fail() != 0) {
		throw CANNOT_OPEN_OUTPUT_FILE();
	}
	
	x_size = 0;
	y_size = 0;
	n_images_written = 0;
	
	// the first 3 * 4 bytes of the file should be written in advance, we will fill them in later
	// by convention they are x_size, y_size, and n_images
	file.write(header, header_length);
	delete[] header;
}


SimpleImageOutputWriter::~SimpleImageOutputWriter() {
	
	if (file.is_open() != 0) {
		flush_and_close();
	}
}



void SimpleImageOutputWriter::write_image(boost::shared_ptr<PALMMatrix<double> > new_image) {
	
	// check whether we should write the cache to disk
	if (image_buffer.size() == N_SIMULTANEOUS_IMAGE_WRITES) {	// we need to flush the cache before accepting a new image
		flush_cache();
	}
	
	// add the new image to the queue
	image_buffer.push(new_image);
	n_images_written++;
}


void SimpleImageOutputWriter::flush_cache() {
	boost::shared_ptr<PALMMatrix<double> > current_image;
	size_t n_pixels, offset;
	
	if (image_buffer.size() == 0) {
		return;
	}
	
	// determine the size of the frames
	current_image = image_buffer.front();
	x_size = current_image->getXSize();
	y_size = current_image->getYSize();
	n_pixels = x_size * y_size;
	
	float *single_image_buffer = new float[n_pixels];	// a temporary buffer for writing a single image
	if (single_image_buffer == NULL) {
		string error;
		error = "unable to allocate buffer in SimpleImageOutputWriter::flush_cache()\r";
		throw OUT_OF_MEMORY(error);
	}
	
	while (image_buffer.size() != 0) {
		offset = 0;
		
		current_image = image_buffer.front();
		
		
		for (size_t j = 0; j < y_size; j++) {
			for (size_t i = 0; i < x_size; i++) {
				single_image_buffer[offset] = (float)current_image->get(i, j);
				offset++;
			}
		}
		
		file.write((char *)single_image_buffer, (n_pixels * sizeof(float)));
		
		// gsl_matrix_free(current_image);
		image_buffer.pop();
	}
	
	delete[] single_image_buffer;
	
}

int SimpleImageOutputWriter::flush_and_close() {
	try {
		flush_cache();
	}
	catch (OUT_OF_MEMORY) {
		file.close();
		string error;
		error = "exiting because of error in flush_cache()\r";
		throw OUT_OF_MEMORY(error);
	}
	
	file.seekp(0);
	file.write((char *)&x_size, sizeof(size_t));
	file.write((char *)&y_size, sizeof(size_t));
	file.write((char *)&n_images_written, sizeof(size_t));
	
	file.close();
	
	return 0;
}


TIFFImageOutputWriter::TIFFImageOutputWriter(const string &rhs,int overwrite, int compression_rhs) {
	// if overwrite is non-zero then we overwrite any file that exists at the output path
	// if it is set to zero then we throw an error and abort instead of overwriting
	
	file_path = rhs;
	compression = compression_rhs;
	
	if (overwrite == 0) {
		ifstream input_test;
		input_test.open(file_path.c_str(), ios::in | ios::binary);
		input_test.close();
		if (input_test.fail() == 0) {
			throw OUTPUT_FILE_ALREADY_EXISTS();	// escape without overwriting
		}
	}
	
	tiff_file = TIFFOpen(file_path.c_str(), "w");
	if (tiff_file == NULL) {
		throw CANNOT_OPEN_OUTPUT_FILE();
	}
	
	x_size = 0;
	y_size = 0;
	n_images_written = 0;
}


TIFFImageOutputWriter::~TIFFImageOutputWriter() {
	
	if (tiff_file != NULL) {
		flush_and_close();
	}
}



void TIFFImageOutputWriter::write_image(boost::shared_ptr<PALMMatrix<double> > new_image) {
	
	// check whether we should write the cache to disk
	if (image_buffer.size() == N_SIMULTANEOUS_IMAGE_WRITES) {	// we need to flush the cache before accepting a new image
		flush_cache();
	}
	
	// add the new image to the queue
	image_buffer.push(new_image);
}


void TIFFImageOutputWriter::flush_cache() {
	boost::shared_ptr<PALMMatrix<double> > current_image;
	size_t n_pixels, offset;
	
	int result;
	uint16_t current_uint16;
	uint32_t current_uint32;
	
	if (image_buffer.size() == 0) {
		return;
	}
	
	// determine the size of the frames
	current_image = image_buffer.front();
	x_size = current_image->getXSize();
	y_size = current_image->getYSize();
	n_pixels = x_size * y_size;
	
	boost::scoped_array<float> scanLine(new float[x_size]);
	
	while (image_buffer.size() != 0) {
		current_image = image_buffer.front();
		
		// make sure that all the image tags have the correct values
		result = TIFFSetField(tiff_file, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
		if (result != 1) {
			string error;
			error = "Unable to set the photometric type for the image at\"";
			error += file_path;
			error += "\"\r";
			throw ERROR_WRITING_FILE_DATA(error);
		}
		
		current_uint32 = x_size;
		result = TIFFSetField(tiff_file, TIFFTAG_IMAGEWIDTH, current_uint32);
		if (result != 1) {
			string error;
			error = "Unable to set the image width for the image at\"";
			error += file_path;
			error += "\"\r";
			throw ERROR_WRITING_FILE_DATA(error);
		}
		
		current_uint32 = y_size;
		result = TIFFSetField(tiff_file, TIFFTAG_IMAGELENGTH, current_uint32);
		if (result != 1) {
			string error;
			error = "Unable to set the image height for the image at\"";
			error += file_path;
			error += "\"\r";
			throw ERROR_WRITING_FILE_DATA(error);
		}
		
		result = TIFFSetField(tiff_file, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);	// floating point values
		if (result != 1) {
			string error;
			error = "Unable to set the SampleFormat for the image at\"";
			error += file_path;
			error += "\"\r";
			throw ERROR_WRITING_FILE_DATA(error);
		}
		
		current_uint16 = 32;
		result = TIFFSetField(tiff_file, TIFFTAG_BITSPERSAMPLE, current_uint16);	// 32 bits per float
		if (result != 1) {
			string error;
			error = "Unable to set the BitsPerSample for the image at\"";
			error += file_path;
			error += "\"\r";
			throw ERROR_WRITING_FILE_DATA(error);
		}
		
		result = TIFFSetField(tiff_file, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
		if (result != 1) {
			string error;
			error = "Unable to set the SubFileType for the image at\"";
			error += file_path;
			error += "\"\r";
			throw ERROR_WRITING_FILE_DATA(error);
		}
		
		current_uint16 = (uint16_t)n_images_written;
		result = TIFFSetField(tiff_file, TIFFTAG_PAGENUMBER, current_uint16);
		if (result != 1) {
			string error;
			error = "Unable to set the PageNumber for the image at\"";
			error += file_path;
			error += "\"\r";
			throw ERROR_WRITING_FILE_DATA(error);
		}
		
		result = TIFFSetField(tiff_file, TIFFTAG_COMPRESSION, compression);
		if (result != 1) {
			string error;
			error = "Unable to set the compression method for the image at\"";
			error += file_path;
			error += "\"\r";
			throw ERROR_WRITING_FILE_DATA(error);
		}
		
		for (size_t j = 0; j < y_size; j++) {
			offset = 0;
			for (size_t i = 0; i < x_size; i++) {
				scanLine[offset] = (float)current_image->get(i, j);
				offset++;
			}
			
			result = TIFFWriteScanline(tiff_file, (char *)scanLine.get(), j);
			if (result != 1) {
				string error;
				error = "There was an error writing a scanline for the image at\"";
				error += file_path;
				error += "\"\r";
				throw ERROR_WRITING_FILE_DATA(error);
			}
			
		}
		
		result = TIFFWriteDirectory(tiff_file);
		if (result != 1) {
			string error;
			error = "Unable to write a directory for the image at\"";
			error += file_path;
			error += "\"\r";
			throw ERROR_WRITING_FILE_DATA(error);
		}
		
		image_buffer.pop();
		++n_images_written;
	}
	
}

int TIFFImageOutputWriter::flush_and_close() {
	if (tiff_file != NULL) {
		flush_cache();
		TIFFClose(tiff_file);
		tiff_file = NULL;
	}
	
	return 0;
}


IgorImageOutputWriter::IgorImageOutputWriter(std::string waveName_rhs, size_t xSize_rhs, size_t ySize_rhs, size_t nImages_rhs, int overwrite_rhs) {
	this->waveName = waveName_rhs;
	this->x_size = xSize_rhs;
	this->y_size = ySize_rhs;
	this->nImagesTotal = nImages_rhs;
	this->n_images_written = 0;
	
	long dimensionSizes[MAX_DIMENSIONS + 1];
	dimensionSizes[0] = this->x_size;
	dimensionSizes[1] = this->y_size;
	dimensionSizes[2] = this->nImagesTotal;
	dimensionSizes[3] = 0;
	int result;
	
	if (overwrite == 0) {
		result = MDMakeWave(&(this->outputWave), this->waveName.c_str(), NULL, dimensionSizes, NT_FP64, 0);
	} else {
		result = MDMakeWave(&(this->outputWave), this->waveName.c_str(), NULL, dimensionSizes, NT_FP64, 1);
	}
	if (result != 0) {
		throw result;
	}
}

void IgorImageOutputWriter::write_image(boost::shared_ptr<PALMMatrix<double> > new_image) {
	long indices[MAX_DIMENSIONS + 1];
	int result;
	double value[2];
	
	indices[2] = n_images_written;
	
	for (size_t j = 0; j < y_size; ++j) {
		for (size_t i = 0; i < x_size; ++i) {
			indices[0] = i;
			indices[1] = j;
			
			value[0] = (*new_image)(i,j);
			result = MDSetNumericWavePointValue(outputWave, indices, value);
			if (result != 0) {
				throw result;
			}
		}
	}
	
	++n_images_written;
}



IgorResultsWriter::IgorResultsWriter(const string rhs) {
	this->waveName = rhs;
	this->totalNumberOfPositions = 0;
}

IgorResultsWriter::~IgorResultsWriter() {
	// if some positions were passed on then write the output
	if (resultsList.size() > 0) {
		WriteResults();
	}
}

void IgorResultsWriter::AppendNewResult(boost::shared_ptr<PALMResults> result) {
	resultsList.push_back(result);
	
	if (result->getFittedPositions()->size() != 0) {
		totalNumberOfPositions += result->getNumberOfFittedPositions();
	}
}

void IgorResultsWriter::WriteResults() {
	waveHndl output_wave;
	long dim_size[3];
	long indices[2];
	int status;
	double value;
	boost::shared_ptr<std::vector<LocalizedPosition> > positions;
	long number_of_positions_in_matrix, offset = 0;
	size_t frameNumber = 0;
	
	dim_size[0] = (long)this->totalNumberOfPositions;
	dim_size[1] = 12;
	dim_size[2] = 0;
	
	// the output wave will be two-dimensional
	// the first dimension contains the different individual positions
	// the second position is structured as follows:
	// n_image	amplitude	radius	x_position	y_position	offset	amp_error	radius_err	x_err	y_err	offset_err	n_iters
	
	// try to make the wave
	status = MDMakeWave(&output_wave, this->waveName.c_str(), NULL, dim_size, NT_FP64, 1);
	if (status != 0)
		throw status;
	
	while (resultsList.size() > 0) {
		
		positions = (resultsList.front())->getFittedPositions();
		
		if (positions->size() == 0) {
			// no positions were found in this frame
			resultsList.pop_front();
			continue;
		}
		
		number_of_positions_in_matrix = positions->size();
		frameNumber = (resultsList.front())->getFrameIndex();
		
		for (long j = 0; j < number_of_positions_in_matrix; j++) {
			// start by storing the index of the matrix
			// so we can tell what image the positions correspond to
			indices[0] = offset;
			indices[1] = 0;
			value = (double)frameNumber;
			MDSetNumericWavePointValue(output_wave, indices, &value);
			
			indices[1] = 1;
			value = positions->at(j).amplitude;
			MDSetNumericWavePointValue(output_wave, indices, &value);
			indices[1] = 2;
			value = positions->at(j).width;
			MDSetNumericWavePointValue(output_wave, indices, &value);
			indices[1] = 3;
			value = positions->at(j).xPos;
			MDSetNumericWavePointValue(output_wave, indices, &value);
			indices[1] = 4;
			value = positions->at(j).yPos;
			MDSetNumericWavePointValue(output_wave, indices, &value);
			indices[1] = 5;
			value = positions->at(j).offset;
			MDSetNumericWavePointValue(output_wave, indices, &value);
			indices[1] = 6;
			value = positions->at(j).amplitudeError;
			MDSetNumericWavePointValue(output_wave, indices, &value);
			indices[1] = 7;
			value = positions->at(j).widthError;
			MDSetNumericWavePointValue(output_wave, indices, &value);
			indices[1] = 8;
			value = positions->at(j).xPosError;
			MDSetNumericWavePointValue(output_wave, indices, &value);
			indices[1] = 9;
			value = positions->at(j).yPosError;
			MDSetNumericWavePointValue(output_wave, indices, &value);
			indices[1] = 10;
			value = positions->at(j).offsetError;
			MDSetNumericWavePointValue(output_wave, indices, &value);
			indices[1] = 11;
			value = (double)positions->at(j).nIterations;
			MDSetNumericWavePointValue(output_wave, indices, &value);
			
			offset++;
		}
		++frameNumber;
		resultsList.pop_front();
	}
}
