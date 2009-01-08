/*
 *  PALM_analysis_classes.cpp
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 25/01/2008.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "PALM_analysis_classes.h"
#include "PALM_analysis.h"
#include "PALM_analysis_IgorXOP.h"


encap_gsl_matrix::encap_gsl_matrix(size_t x, size_t y) {
	matrix = gsl_matrix_alloc(x, y);
	if (matrix == NULL) {
		string error;
		error = "unable to allocate matrix in encap_gsl_matrix::encap_gsl_matrix()\r";
		throw OUT_OF_MEMORY(error);
	}
	x_size = x;
	y_size = y;
}

encap_gsl_matrix::~encap_gsl_matrix() {
	if (matrix != NULL) {
		gsl_matrix_free(matrix);
	}
}

void encap_gsl_matrix::set(size_t x, size_t y, double value) {
	#ifndef ENCAP_GSL_RANGE_CHECK_OFF
	if ((x > x_size - 1) || (y > y_size - 1)) {
		error = "Out of range in encap_gsl_matrix::set()\r";
		throw std::range_error(error);
	}
	#endif
	
	gsl_matrix_set(matrix, x, y, value);
}

void encap_gsl_matrix::operator=(gsl_matrix* rhs) {
	if (rhs == matrix)
		return;
	
	if (matrix != NULL) {
		gsl_matrix_free(matrix);
	}
	
	matrix = rhs;
	
	if (rhs != NULL) {
		x_size = matrix->size1;
		y_size = matrix->size2;
	}
	
}

double encap_gsl_matrix::get(size_t x, size_t y) {
	#ifndef ENCAP_GSL_RANGE_CHECK_OFF
	if ((x > x_size - 1) || (y > y_size - 1)) {
		error = "Out of range in encap_gsl_matrix::get()\r";
		throw std::range_error(error);
	}
	#endif
	
	return gsl_matrix_get(matrix, x, y);
}

encap_gsl_matrix_uchar::encap_gsl_matrix_uchar(size_t x, size_t y) {
	matrix = gsl_matrix_uchar_alloc(x, y);
	if (matrix == NULL) {
		string error;
		error = "unable to allocate matrix in encap_gsl_matrix::encap_gsl_matrix()\r";
		throw OUT_OF_MEMORY(error);
	}
	x_size = x;
	y_size = y;
}

encap_gsl_matrix_uchar::~encap_gsl_matrix_uchar() {
	if (matrix != NULL) {
		gsl_matrix_uchar_free(matrix);
	}
}

void encap_gsl_matrix_uchar::set(size_t x, size_t y, unsigned char value) {
	#ifndef ENCAP_GSL_RANGE_CHECK_OFF
	if ((x > x_size - 1) || (y > y_size - 1)) {
		error = "Out of range in encap_gsl_matrix::set()\r";
		throw std::range_error(error);
	}
	#endif
	
	gsl_matrix_uchar_set(matrix, x, y, value);
}

void encap_gsl_matrix_uchar::operator=(gsl_matrix_uchar * rhs) {
	if (rhs == matrix)
		return;
	
	if (matrix != NULL) {
		gsl_matrix_uchar_free(matrix);
	}
	
	matrix = rhs;
	
	if (rhs != NULL) {
		x_size = matrix->size1;
		y_size = matrix->size2;
	}
	
}

unsigned char encap_gsl_matrix_uchar::get(size_t x, size_t y) {
	#ifndef ENCAP_GSL_RANGE_CHECK_OFF
	if ((x > x_size - 1) || (y > y_size - 1)) {
		error = "Out of range in encap_gsl_matrix::get()\r";
		throw std::range_error(error);
	}
	#endif
	
	return gsl_matrix_uchar_get(matrix, x, y);
}

encap_gsl_matrix_long::encap_gsl_matrix_long(size_t x, size_t y) {
	matrix = gsl_matrix_long_alloc(x, y);
	if (matrix == NULL) {
		string error;
		error = "unable to allocate matrix in encap_gsl_matrix::encap_gsl_matrix()\r";
		throw OUT_OF_MEMORY(error);
	}
	x_size = x;
	y_size = y;
}

encap_gsl_matrix_long::~encap_gsl_matrix_long() {
	if (matrix != NULL) {
		gsl_matrix_long_free(matrix);
	}
}

void encap_gsl_matrix_long::set(size_t x, size_t y, long value) {
	#ifndef ENCAP_GSL_RANGE_CHECK_OFF
	if ((x > x_size - 1) || (y > y_size - 1)) {
		error = "Out of range in encap_gsl_matrix::set()\r";
		throw std::range_error(error);
	}
	#endif
	
	gsl_matrix_long_set(matrix, x, y, value);
}

void encap_gsl_matrix_long::operator=(gsl_matrix_long * rhs) {
	if (rhs == matrix)
		return;
	
	if (matrix != NULL) {
		gsl_matrix_long_free(matrix);
	}
	
	matrix = rhs;
	
	if (rhs != NULL) {
		x_size = matrix->size1;
		y_size = matrix->size2;
	}
	
}

long encap_gsl_matrix_long::get(size_t x, size_t y) {
	#ifndef ENCAP_GSL_RANGE_CHECK_OFF
	if ((x > x_size - 1) || (y > y_size - 1)) {
		error = "Out of range in encap_gsl_matrix::get()\r";
		throw std::range_error(error);
	}
	#endif
	
	return gsl_matrix_long_get(matrix, x, y);
}

ImageLoader::ImageLoader() {
	path.assign("");
	total_number_of_images = 0;
	x_size = 0;
	y_size = 0;
	header_length = 0;
	image_cache.reserve(N_SIMULTANEOUS_IMAGE_LOADS);
	images_in_cache.reserve(N_SIMULTANEOUS_IMAGE_LOADS);
}

ImageLoader::~ImageLoader() {
	if (file.is_open() == 1)
		file.close();
	
	/*// free any images that may remain in the cache
	for (unsigned long i = 0; i < N_SIMULTANEOUS_IMAGE_LOADS; i++) {
		if (image_cache[i] != NULL) {
			gsl_matrix_free(image_cache[i]);
			image_cache[i] = NULL;
		}
	}*/
}


ImageLoaderSPE::ImageLoaderSPE(string rhs) {
	path = rhs;
	image_cache_size = N_SIMULTANEOUS_IMAGE_LOADS;
	
	header_length = 4100;
	
	file.open(path.c_str(), ios::binary | ios::in);
	if (file.fail() == 1) {
		throw CANNOT_OPEN_FILE();
	}
	
	parse_header_information();
	
	image_cache.resize(image_cache_size, boost::shared_ptr<encap_gsl_matrix>());
	images_in_cache.resize(image_cache_size, (unsigned long)-1);
}

ImageLoaderSPE::ImageLoaderSPE(string rhs, unsigned long image_cache_size_rhs) {
	path = rhs;
	image_cache_size = image_cache_size_rhs;
	if (image_cache_size > N_SIMULTANEOUS_IMAGE_LOADS) {
		image_cache.reserve(image_cache_size);
	}
	
	header_length = 4100;
	
	file.open(path.c_str(), ios::binary | ios::in);
	if (file.fail() == 1) {
		throw CANNOT_OPEN_FILE();
	}
	
	parse_header_information();
	
	image_cache.resize(image_cache_size, boost::shared_ptr<encap_gsl_matrix>());
	images_in_cache.resize(image_cache_size, (unsigned long)-1);
}

ImageLoaderSPE::~ImageLoaderSPE() {
	if (file.is_open() == 1)
		file.close();
}

int ImageLoaderSPE::parse_header_information() {
	long current_bytes = 0;	// assume that this is a four-byte variable
	char byte_reader1, byte_reader2, byte_reader3, byte_reader4;	// we want to only read 1 byte at a time
	
	// is there a file open?
	if (file.is_open() == 0) {
		return -1;
	}
	
	file.seekg(42);
	file.get(byte_reader1);	// we know that the values are little-endian
	file.get(byte_reader2);
	current_bytes = 0x000000FF & byte_reader2;
	current_bytes *= 256;
	current_bytes = current_bytes | (0x000000FF & byte_reader1);
	
	x_size = current_bytes;
	current_bytes = 0;
	
	file.seekg(656);
	file.get(byte_reader1);	// we know that the values are little-endian
	file.get(byte_reader2);
	current_bytes = 0x000000FF & byte_reader2;
	current_bytes *= 256;
	current_bytes = current_bytes | (0x000000FF & byte_reader1);
	
	y_size = current_bytes;
	current_bytes = 0;
	
	file.seekg(1446);
	file.get(byte_reader1);	// we know that the values are little-endian
	file.get(byte_reader2);
	file.get(byte_reader3);	// we know that the values are little-endian
	file.get(byte_reader4);
	current_bytes = 0x000000FF  & byte_reader4;
	current_bytes *= 256;
	current_bytes = current_bytes | (0x000000FF & byte_reader3);
	current_bytes *= 256;
	current_bytes = current_bytes | (0x000000FF & byte_reader2);
	current_bytes *= 256;
	current_bytes = current_bytes | (0x000000FF & byte_reader1);
	
	total_number_of_images = current_bytes;
	current_bytes = 0;
	
	file.seekg(108);
	file.get(byte_reader1);	// we know that the values are little-endian
	file.get(byte_reader2);
	current_bytes = 0x000000FF & byte_reader2;
	current_bytes *= 256;
	current_bytes = current_bytes | (0x000000FF & byte_reader1);
	
	storage_type = (int)current_bytes;
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
	
	return 0;
}

boost::shared_ptr<encap_gsl_matrix> ImageLoaderSPE::get_nth_image(const unsigned long n) {	
	unsigned long offset;
	long current_long = 0;
	float current_float = 0;
	short current_short = 0;
	unsigned short current_unsigned_short = 0;
	string error;
	boost::shared_ptr<encap_gsl_matrix> image;
	// gsl_matrix *image_copy;
	
	if (n >= total_number_of_images)
		throw IMAGE_INDEX_BEYOND_N_IMAGES();
	
	// is there a file open?
	if (file.is_open() != 1)
		throw GET_NTH_IMAGE_FILE_NOT_OPEN();
	
	// is the requested image in the cache?
	// if it is then we return a copy
	for (unsigned long i = 0; i < image_cache_size; i++) {
		if (images_in_cache[i] == n) {
			return image_cache[i];
	//		image_copy = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(x_size, y_size));
			
	//		gsl_matrix_memcpy(image_copy, image_cache[i]);
			
	//		return image_copy;
		}
	}
	
	// if we are here then the image wasn't in the cache
	// we will load a new sequence of images in the cache, starting with the one that is requested
	
	// first we free the images that are already in the cache
	// this is safe because we only pass copies to other functions
	for (unsigned long i = 0; i < image_cache_size; i++) {
		images_in_cache[i] = -1;
	}
	
	// now load the new set of images
	
	unsigned long n_bytes_in_single_image;
	unsigned long cache_offset;
	
	// determine how big we have to make the single image buffer
	switch(storage_type) {
		case 0:	// 4 byte float
			n_bytes_in_single_image = x_size * y_size * 4;
			break;
		case 1:	// 4-byte long
			n_bytes_in_single_image = x_size * y_size * 4;
			break;
		case 2:	// 2 byte signed short
			n_bytes_in_single_image = x_size * y_size * 2;
			break;
		case 3:	// 2 byte unsigned short
			n_bytes_in_single_image = x_size * y_size * 2;
			break;
		default:
			throw CANNOT_DETERMINE_SPE_STORAGE_TYPE();
			break;
	}
	
		boost::scoped_array<float> single_image_buffer_float(new float[x_size * y_size]);
		boost::scoped_array<char> single_image_buffer(new char[n_bytes_in_single_image]);
	
	for (unsigned long i = n; i < n + image_cache_size; i++) {
		
		if (i >= total_number_of_images) {
			continue;	// there are no more images to load into this cache location
			// by freeing the locations above we have already made sure that the arrays are set to a NULL pointer and -1
		}
		
		image = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(x_size, y_size));
		
		switch(storage_type) {
			case 0:	// 4 byte float
				offset = header_length + i * (x_size) * (y_size) * 4;
				break;
			case 1:	// 4-byte long
				offset = header_length + i * (x_size) * (y_size) * 4;
				break;
			case 2:	// 2 byte signed short
				offset = header_length + i * (x_size) * (y_size) * 2;
				break;
			case 3:	// 2 byte unsigned short
				offset = header_length + i * (x_size) * (y_size) * 2;
				break;
			default:
				throw CANNOT_DETERMINE_SPE_STORAGE_TYPE();	// if we have already loaded some images in the cache then this causes a memory leak
															// fortunately that's not a real problem because this would have been caught already while
															// we were determining the number of bytes in each image
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
			throw ERROR_READING_FILE_DATA(error);
		}
		
		switch(storage_type) {
			case 0:	// 4-byte float
					// this is currently only safe on little-endian systems!
					for (unsigned long j  = 0; j < y_size; j++) {
						for (unsigned long i = 0; i < x_size; i++) {
							//	file.read((char *)&current_float, sizeof(float));
							
							current_float = single_image_buffer_float[cache_offset];
							
							image->set(i, j, (double)current_float);
							
							cache_offset++;
						}
					}
					break;
				
				case 1:	// 4-byte long
					for (unsigned long j  = 0; j < y_size; j++) {
						for (unsigned long i = 0; i < x_size; i++) {
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
					
				case 2:	// 2-byte signed short
					for (unsigned long j  = 0; j < y_size; j++) {
						for (unsigned long i = 0; i < x_size; i++) {
							current_short = 0x000000FF & single_image_buffer[cache_offset + 1];	// little endian
							current_short *= 256;
							current_short = current_short | (0x000000FF & single_image_buffer[cache_offset + 0]);
							
							image->set(i, j, (double)current_short);
							
							cache_offset += 2;
						}
					}
					break;
					
				case 3: // 2-byte unsigned short
					for (unsigned long j  = 0; j < y_size; j++) {
						for (unsigned long i = 0; i < x_size; i++) {
						//	file.read((char *)&current_unsigned_short, 2);
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
					throw CANNOT_DETERMINE_SPE_STORAGE_TYPE();
					break;
				
		}
		
		// store the image at the correct location in the cache
		image_cache[i - n] = image;
		images_in_cache[i - n] = i;
	}
	
	file.seekg(0);
	
	return image_cache[0];
	
}

ImageLoaderAndor::ImageLoaderAndor(string rhs) {
	path = rhs;
	image_cache_size = N_SIMULTANEOUS_IMAGE_LOADS;
	
	file.open(path.c_str(), ios::binary | ios::in);
	if (file.fail() == 1) {
		throw CANNOT_OPEN_FILE();
	}
	
	parse_header_information();
	
	image_cache.resize(image_cache_size, boost::shared_ptr<encap_gsl_matrix>());
	images_in_cache.resize(image_cache_size, (unsigned long)-1);
}

ImageLoaderAndor::ImageLoaderAndor(string rhs, unsigned long image_cache_size_rhs) {
	path = rhs;
	image_cache_size = image_cache_size_rhs;
	if (image_cache_size > N_SIMULTANEOUS_IMAGE_LOADS) {
		image_cache.reserve(image_cache_size);
	}
	
	file.open(path.c_str(), ios::binary | ios::in);
	if (file.fail() == 1) {
		throw CANNOT_OPEN_FILE();
	}
	
	parse_header_information();
	
	image_cache.resize(image_cache_size, boost::shared_ptr<encap_gsl_matrix>());
	images_in_cache.resize(image_cache_size, (unsigned long)-1);
}

ImageLoaderAndor::~ImageLoaderAndor() {
	if (file.is_open() == 1)
		file.close();
	// free up any images that may remain in the cache
	/*for (unsigned long i = 0; i < N_SIMULTANEOUS_IMAGE_LOADS; i++) {
		if (image_cache[i] != NULL) {
			gsl_matrix_free(image_cache[i]);
			image_cache[i] = NULL;		// we need to do this as otherwise the ImageLoader destructor tries to free the same memory again
		}
	}*/
}

int ImageLoaderAndor::parse_header_information() {
	char header_buffer_char[1024];
	string header_buffer;
	string header_subset;
	size_t start_of_info, end_of_info;
	// size_t start_of_data;
	int temp;
	// unsigned long file_position;
	// unsigned long length;
	
	// we need to get storage type out of the way as Andor doesn't use this
	storage_type = 0;
	
	// the information that we are looking for is in the 22nd line (starting from 1) of the file
	int image_params_line = 21;
	
	for (int i = 0; i < image_params_line; i++) {
		file.getline(header_buffer_char, 1024);
	}
	
	// now we read a few hundred bytes at the beginning of the file into the buffer
	
	// file_position = file.tellg();
	// file.read(header_buffer_char, 24000);
	// header_buffer_char[24000] = '\0';
	
	// we have to make sure that there are no terminating \0 characters
	// for (int i = 0; i < 24000; i++) {
	//	if (header_buffer_char[i] == '\0') {
	//		header_buffer_char[i] = '0';
	//	}
	//}
	
	// now we read the line containing the first part of the data
	file.getline(header_buffer_char, 1024);
	header_buffer.assign(header_buffer_char);
	
	// read the next part
	file.getline(header_buffer_char, 1024);
	header_buffer.append(header_buffer_char);
	
	// length = header_buffer.length();
	
	start_of_info = header_buffer.find("65538");
	end_of_info = header_buffer.find("65538", start_of_info + 1);
	
	header_subset = header_buffer.substr(start_of_info, end_of_info);
	stringstream header_data(header_subset, (ios::in | ios::out));
	
	header_data >> temp;	// this extracts 65538
	header_data >> temp;
	
	header_data >> y_size;
	header_data >> x_size;
	
	header_data >> temp;
	
	header_data >> total_number_of_images;
	
	// now we have one additional line and a bunch of zeroes (equal to the number of images) on separate lines until the file data begins.
	// we need to find the position of the start of the data
	
	// the number of zeroes is equal to the number of frames in the image
	for (unsigned long i = 0;  i < total_number_of_images; i++) {
		file.getline(header_buffer_char, 1024);
	}
	
	header_length = file.tellg();
	
//	start_of_data = header_buffer.rfind("0\n0\n") + 4;
//	header_length = start_of_data + file_position + 1;
	
	// did some error happen while reading the file?
	if (file.fail() != 0) {
		string error;
		error = "Error parsing the header information in \"";
		error += path;
		error += "\" assuming the Andor format\r";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	file.seekg(0);
	
	return 0;
}
	
boost::shared_ptr<encap_gsl_matrix> ImageLoaderAndor::get_nth_image(const unsigned long n) {	
	unsigned long offset;
	float current_float = 0;
	
	
	boost::shared_ptr<encap_gsl_matrix> image;
	// gsl_matrix *image_copy;
	
	if (n >= total_number_of_images)
		throw IMAGE_INDEX_BEYOND_N_IMAGES();
	
	// is there a file open?
	if (file.is_open() != 1)
		throw GET_NTH_IMAGE_FILE_NOT_OPEN();
	
	// is the requested image in the cache?
	// if it is then we return a copy
	for (unsigned long i = 0; i < image_cache_size; i++) {
		if (images_in_cache[i] == n) {
			
	//		image_copy = gsl_matrix_alloc(x_size, y_size);
			
	//		if (image_copy == NULL)
	//			throw OUT_OF_MEMORY();
			
	//		gsl_matrix_memcpy(image_copy, image_cache[i]);
			
	//		return image_copy;
			return image_cache[i];
		}
	}
	
	// if we are here then the image wasn't in the cache
	// we will load a new sequence of images in the cache, starting with the one that is requested
	
	// first we free the images that are already in the cache
	// this is safe because we only pass copies to other functions
	for (unsigned long i = 0; i < image_cache_size; i++) {
		/*if (image_cache[i] != NULL) {
			gsl_matrix_free(image_cache[i]);
			image_cache[i] = NULL;
		}*/
		images_in_cache[i] = -1;
	}
	
	// now load the new set of images
	
	boost::scoped_array<float> single_image_buffer(new float[x_size * y_size]);
	unsigned long cache_offset;
	
	for (unsigned long i = n; i < n + image_cache_size; i++) {
		
		if (i >= total_number_of_images) {
			continue;	// there are no more images to load into this cache location
			// by freeing the locations above we have already made sure that the arrays are set to a NULL pointer and -1
		}
		
		image = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(x_size, y_size));
		
		offset = header_length + i * (x_size) * (y_size) * sizeof(float);
		file.seekg(offset);
		
		cache_offset = 0;
		
		file.read((char *)single_image_buffer.get(), (x_size * y_size * sizeof(float)));
		if (file.fail() != 0) {
			string error;
			error = "Error trying to read image data from \"";
			error += path;
			error += "\" assuming the Andor format\r";
			throw ERROR_READING_FILE_DATA(error);
		}
			
		
		// this is currently only safe on little-endian systems!
		for (unsigned long j  = 0; j < y_size; j++) {
			for (unsigned long i = 0; i < x_size; i++) {
				current_float = single_image_buffer[cache_offset];
				image->set(i, j, (double)current_float);
				cache_offset++;
			}
		}
		
		// store the image at the correct location in the cache
		image_cache[i - n] = image;
		images_in_cache[i - n] = i;
	}
	
	file.seekg(0);
	
	// the originally requested image is at the start of the cache
	// however, we will return a copy of the matrix, so that no synchronization problems can occur
	// if another thread were to access a new image at a different time
	/*image_copy = gsl_matrix_alloc(x_size, y_size);
	if (image_copy == NULL) {
		throw OUT_OF_MEMORY();
	}
	gsl_matrix_memcpy(image_copy, image_cache[0]);
	return image_copy;*/
	
	return image_cache[0];
	
}


ImageLoaderHamamatsu::ImageLoaderHamamatsu(string rhs) {
	path = rhs;
	image_cache_size = N_SIMULTANEOUS_IMAGE_LOADS;
	
	file.open(path.c_str(), ios::binary | ios::in);
	if (file.fail() == 1) {
		throw CANNOT_OPEN_FILE();
	}
	
	parse_header_information();
	
	image_cache.resize(image_cache_size, boost::shared_ptr<encap_gsl_matrix>());
	images_in_cache.resize(image_cache_size, (unsigned long)-1);
}

ImageLoaderHamamatsu::ImageLoaderHamamatsu(string rhs, unsigned long image_cache_size_rhs) {
	path = rhs;
	image_cache_size = image_cache_size_rhs;
	if (image_cache_size > N_SIMULTANEOUS_IMAGE_LOADS) {
		image_cache.reserve(image_cache_size);
	}
	
	file.open(path.c_str(), ios::binary | ios::in);
	if (file.fail() == 1) {
		throw CANNOT_OPEN_FILE();
	}
	
	parse_header_information();
	
	image_cache.resize(image_cache_size, boost::shared_ptr<encap_gsl_matrix>());
	images_in_cache.resize(image_cache_size, (unsigned long)-1);
}

ImageLoaderHamamatsu::~ImageLoaderHamamatsu() {
	if (file.is_open() == 1)
		file.close();
}


int ImageLoaderHamamatsu::parse_header_information() {
	char header_buffer_char[1024];	// too large, but we'd better be safe
//	string header_buffer;
//	string header_subset;
//	size_t start_of_info, end_of_info;
//	size_t start_of_data;
//	int temp;
	string header_string;
	unsigned long wasabi_position;
	
	char byte_reader1, byte_reader2;
	total_number_of_images = 0;
	x_size = 0;
	y_size = 0;
	storage_type = 0;
	
	
	// bytes 4-7 contain the x and y size as UINT16
	file.seekg(4);
	file.get(byte_reader1);
	file.get(byte_reader2);
	
	x_size = 0x00000000FF & byte_reader2;	// little-endian
	x_size *= 256;
	x_size = x_size | (0x00000000FF & byte_reader1);
	
	file.get(byte_reader1);
	file.get(byte_reader2);
	y_size = 0x000000FF & byte_reader2;	// little-endian
	y_size *= 256;
	y_size = y_size | (0x00000000FF & byte_reader1);
	
	// get the number of images, stored at bytes 14 and 15
	file.seekg(14);
	file.get(byte_reader1);
	file.get(byte_reader2);
	
	total_number_of_images = 0x00000000FF & byte_reader2;
	total_number_of_images *= 256;
	total_number_of_images = total_number_of_images | (0x00000000FF & byte_reader1);
	
	
	// now we need to determine the length of the header
	// we will do this by looking for the string "~WASABI~" in the file
	// first we load a set of data
	file.seekg(0);
	file.read(header_buffer_char, 1023);
	header_buffer_char[1023] = '\0';
	// we have to remove extraneous NULL characters that may be present in the string we read from the file
	// this will cause the string to be interpreted as a very short array
	for (int i = 0; i < 1023; i++) {
		if (header_buffer_char[i] == '\0')
			header_buffer_char[i] = 1;
	}
	
	header_string.assign(header_buffer_char);
	
	// there seem to be different versions of the files, but they all follow the same stupid principle
	// first we check if it's one of the old files
	
	wasabi_position = header_string.find("~WASABI~", 0);
	if (wasabi_position == (unsigned long)(-1)) {	// we didn't the "~WASABI~", it's probably a new file
		wasabi_position = header_string.find("~Hokawo~", 0);
	}
	
	header_length = wasabi_position + 8;	// we found the length of the header
	
	// was there an error reading the file?
	if (file.fail() != 0) {
		string error;
		error = "Error parsing the header information in \"";
		error += path;
		error += "\" assuming the Hamamatsu format\r";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	file.seekg(0);
	
	return 0;
}


boost::shared_ptr<encap_gsl_matrix> ImageLoaderHamamatsu::get_nth_image(const unsigned long n) {	
	unsigned long offset;
//	int current_int = 0;
//	char byte_reader1, byte_reader2;
	boost::shared_ptr<encap_gsl_matrix> image;
	// gsl_matrix *image_copy;
	
	if (n >= total_number_of_images)
		throw IMAGE_INDEX_BEYOND_N_IMAGES();
	
	// is there a file open?
	if (file.is_open() != 1)
		throw GET_NTH_IMAGE_FILE_NOT_OPEN();
	
	
	// is the requested image in the cache?
	// if it is then we return a copy
	for (unsigned long i = 0; i < image_cache_size; i++) {
		if (images_in_cache[i] == n) {
			
		//	image_copy = gsl_matrix_alloc(x_size, y_size);
			
		//	if (image_copy == NULL)
		//		throw OUT_OF_MEMORY();
			
		//	gsl_matrix_memcpy(image_copy, image_cache[i]);
			
		//	return image_copy;
			return image_cache[i];
		}
	}
	
	// if we are here then the image wasn't in the cache
	// we will load a new sequence of images in the cache, starting with the one that is requested
	
	// first we free the images that are already in the cache
	// this is safe because we only pass copies to other functions
	for (unsigned long i = 0; i < image_cache_size; i++) {
		/*if (image_cache[i] != NULL) {
			gsl_matrix_free(image_cache[i]);
			image_cache[i] = NULL;
		}*/
		images_in_cache[i] = -1;
	}
	
	// now load the new set of images
	
	// we are going to load the images first into a char buffer
	// because we assume that a char s equal to one byte
	
	unsigned int current_uint;
	unsigned long n_bytes_per_image = x_size * y_size * 2;
	boost::scoped_array<char> single_image_buffer(new char[n_bytes_per_image]);
	unsigned long cache_offset;
	
	for (unsigned long i = n; i < n + image_cache_size; i++) {
		
		if (i >= total_number_of_images) {
			continue;	// there are no more images to load into this cache location
			// by freeing the locations above we have already made sure that the arrays are set to a NULL pointer and -1
		}
		
		image = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(x_size, y_size));
		
		offset = (i + 1) * header_length + i * (x_size) * (y_size) * 2;	// assume a 16-bit format
		file.seekg(offset);
		
		file.read(single_image_buffer.get(), n_bytes_per_image);
		if (file.fail() != 0) {
			string error;
			error = "Error reading image data from \"";
			error += path;
			error += "\" assuming the Hamamatsu format\r";
			throw ERROR_READING_FILE_DATA(error);
		}
		
		cache_offset = 0;
		
		for (unsigned long j  = 0; j < y_size; j++) {
			for (unsigned long i = 0; i < x_size; i++) {
				current_uint = 0x000000FF & single_image_buffer[cache_offset + 1];	// little endian
				current_uint *= 256;
				current_uint = current_uint | (0x000000FF & single_image_buffer[cache_offset]);
				
				image->set(i, j, (double)current_uint);
				
				cache_offset += 2;	// TWO bytes per value (UINT16)
			}
		}
		
		
		// store the image at the correct location in the cache
		image_cache[i - n] = image;
		images_in_cache[i - n] = i;
	}
	
	file.seekg(0);
	return image_cache[0];
}


SimpleImageLoader::SimpleImageLoader(string rhs) {
	path = rhs;
	image_cache_size = N_SIMULTANEOUS_IMAGE_LOADS;
	
	file.open(path.c_str(), ios::binary | ios::in);
	if (file.fail() == 1) {
		throw CANNOT_OPEN_FILE();
	}
	
	parse_header_information();
	
	image_cache.resize(image_cache_size, boost::shared_ptr<encap_gsl_matrix>());
	images_in_cache.resize(image_cache_size, (unsigned long)-1);
}

SimpleImageLoader::SimpleImageLoader(string rhs, unsigned long image_cache_size_rhs) {
	path = rhs;
	image_cache_size = image_cache_size_rhs;
	if (image_cache_size > N_SIMULTANEOUS_IMAGE_LOADS) {
		image_cache.reserve(image_cache_size);
	}
	
	file.open(path.c_str(), ios::binary | ios::in);
	if (file.fail() == 1) {
		throw CANNOT_OPEN_FILE();
	}
	
	parse_header_information();
	
	image_cache.resize(image_cache_size, boost::shared_ptr<encap_gsl_matrix>());
	images_in_cache.resize(image_cache_size, (unsigned long)-1);
}

SimpleImageLoader::~SimpleImageLoader() {
	if (file.is_open() == 1) {
		file.close();
	}
}

boost::shared_ptr<encap_gsl_matrix> SimpleImageLoader::get_nth_image(const unsigned long n) {

	// gsl_matrix *image_copy;
	
	if (n >= total_number_of_images)
		throw IMAGE_INDEX_BEYOND_N_IMAGES();
	
	// is there a file open?
	if (file.is_open() != 1)
		throw GET_NTH_IMAGE_FILE_NOT_OPEN();
	
	
	// is the requested image in the cache?
	// if it is then we return a copy
	for (unsigned long i = 0; i < image_cache_size; i++) {
		if (images_in_cache[i] == n) {
			return image_cache[i];
		}
	}
	
	// if we are here then the image wasn't in the cache
	// we will load a new sequence of images in the cache, starting with the one that is requested
	
	// first we free the images that are already in the cache
	// this is safe because we only pass copies to other functions
	for (unsigned long i = 0; i < image_cache_size; i++) {
		images_in_cache[i] = -1;
	}
	
	// now load the new set of images
	unsigned long offset, array_offset;
	boost::shared_ptr<encap_gsl_matrix> new_image;
	boost::scoped_array<float> single_image_buffer(new float[x_size * y_size]);
	
	for (unsigned long i = n; i < n + image_cache_size; i++) {
		
		if (i >= total_number_of_images) {
			continue;	// there are no more images to load into this cache location
			// by freeing the locations above we have already made sure that the arrays are set to a NULL pointer and -1
		}
		
		new_image = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(x_size, y_size));
		
		offset = 3 * sizeof(unsigned long) + i * (x_size) * (y_size) * sizeof(float);
		file.seekg(offset);
		
		array_offset = 0;
		
		// this is currently only safe on little-endian systems!
		
		file.read((char *)single_image_buffer.get(), (x_size * y_size * sizeof(float)));
		if (file.fail() != 0) {
			string error;
			error = "Error reading image data from \"";
			error += path;
			error += "\" assuming the simple image format\r";
			throw ERROR_READING_FILE_DATA(error);
		}
		
		for (unsigned long j  = 0; j < y_size; j++) {
			for (unsigned long i = 0; i < x_size; i++) {
				new_image->set(i, j, single_image_buffer[array_offset]);
				array_offset++;
			}
		}
		
		// store the image at the correct location in the cache
		image_cache[i - n] = new_image;
		images_in_cache[i - n] = i;
	}
	
	file.seekg(0);
	return image_cache[0];
}


int SimpleImageLoader::parse_header_information() {
	
	file.seekg(0);
	
	// the first 12 bytes of the file contain this information
	file.read((char *)&x_size, sizeof(unsigned long));
	file.read((char *)&y_size, sizeof(unsigned long));
	file.read((char *)&total_number_of_images, sizeof(unsigned long));
	
	if (file.fail() != 0) {
		string error;
		error = "Error parsing the header information in \"";
		error += path;
		error += "\" assuming the simple image format\r";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	return 0;
	
}
	

ImageLoaderTIFF_Igor::ImageLoaderTIFF_Igor(string rhs) {
	int result;
	
	path = rhs;
	image_cache_size = N_SIMULTANEOUS_IMAGE_LOADS;
	
	file.open(path.c_str(), ios::binary | ios::in);
	if (file.fail() == 1) {
		throw CANNOT_OPEN_FILE();
	}
	
	file.close();
	
	native_Igor_FilePath = path;
	
	#ifdef _MACINTOSH_	// if we are on a mac then convert the POSIX path back to a native HFS path
						// converting the paths from HFS to POSIX is unfortunately quite messy
	
		CFStringRef POSIX_Path_ref;
		CFURLRef urlRef;
		char buffer[4096];
		
		POSIX_Path_ref = CFStringCreateWithCString(kCFAllocatorDefault, path.c_str(), kCFStringEncodingASCII);
		
		urlRef = CFURLCreateWithFileSystemPath(kCFAllocatorDefault, POSIX_Path_ref, kCFURLPOSIXPathStyle, false);
		
		POSIX_Path_ref = CFURLCopyFileSystemPath(urlRef, kCFURLHFSPathStyle);
		
		result = CFStringGetCString(POSIX_Path_ref, buffer, 4096, kCFStringEncodingASCII);
		if (result != true) {
			throw CANNOT_OPEN_FILE();
		}
		
		native_Igor_FilePath = buffer;
		
	#endif
	
	#ifdef _WINDOWS_
		// we need to escape the backslashes in the filepath so that it can be passed directly to Igor
		for (unsigned long i = 0; i < native_Igor_FilePath.length(); i++) {
			if (native_Igor_FilePath[i] == '\\') {
				native_Igor_FilePath.insert(i, "\\");
				i++;
			}
		}
	#endif
	
	parse_header_information();
	
	tiff_cache_wave = NULL;
	images_in_cache.resize(image_cache_size, (unsigned long)-1);
	
}

ImageLoaderTIFF_Igor::ImageLoaderTIFF_Igor(string rhs, unsigned long image_cache_size_rhs) {
	int result;
	image_cache_size = image_cache_size_rhs;
	if (image_cache_size > N_SIMULTANEOUS_IMAGE_LOADS) {
		image_cache.reserve(image_cache_size);
	}
	
	path = rhs;
	
	file.open(path.c_str(), ios::binary | ios::in);
	if (file.fail() == 1) {
		throw CANNOT_OPEN_FILE();
	}
	
	file.close();
	
	native_Igor_FilePath = path;
	
#ifdef _MACINTOSH_	// if we are on a mac then convert the POSIX path back to a native HFS path
	// converting the paths from HFS to POSIX is unfortunately quite messy
	
	CFStringRef POSIX_Path_ref;
	CFURLRef urlRef;
	char buffer[4096];
	
	POSIX_Path_ref = CFStringCreateWithCString(kCFAllocatorDefault, path.c_str(), kCFStringEncodingASCII);
	
	urlRef = CFURLCreateWithFileSystemPath(kCFAllocatorDefault, POSIX_Path_ref, kCFURLPOSIXPathStyle, false);
	
	POSIX_Path_ref = CFURLCopyFileSystemPath(urlRef, kCFURLHFSPathStyle);
	
	result = CFStringGetCString(POSIX_Path_ref, buffer, 4096, kCFStringEncodingASCII);
	if (result != true) {
		throw CANNOT_OPEN_FILE();
	}
	
	native_Igor_FilePath = buffer;
	
#endif
	
#ifdef _WINDOWS_
	// we need to escape the backslashes in the filepath so that it can be passed directly to Igor
	for (unsigned long i = 0; i < native_Igor_FilePath.length(); i++) {
		if (native_Igor_FilePath[i] == '\\') {
			native_Igor_FilePath.insert(i, "\\");
			i++;
		}
	}
#endif
	
	parse_header_information();
	
	tiff_cache_wave = NULL;
	images_in_cache.resize(image_cache_size, (unsigned long)-1);
	
}


ImageLoaderTIFF_Igor::~ImageLoaderTIFF_Igor() {
	int result;
	if (tiff_cache_wave != NULL) {
		result = KillWave(tiff_cache_wave);
		if (result != 0) {
			throw result;
		}
	}
}
		

int ImageLoaderTIFF_Igor::parse_header_information() {
	// we let Igor do the parsing for us
	string IgorCommand;
	double value[2];
	int result;
	
	IgorCommand = "ImageFileInfo \"";
	IgorCommand += native_Igor_FilePath;
	IgorCommand += "\"";
	// XOPNotice(IgorCommand.c_str());
	
	result = XOPSilentCommand(IgorCommand.c_str());
	if (result != 0) {
		string error;
		error = "Error processing ImageFileInfo while loading the TIFF file at\"";
		error += path;
		error += "\" using Igor\r";
		throw result;
	}
	
	result = FetchNumVar("V_Flag", &value[0], &value[1]);
	if (result == -1) {
		string error;
		error = "Error processing V_Flag from ImageFileInfo while loading the TIFF file at\"";
		error += path;
		error += "\" using Igor\r";
		throw ERROR_READING_FILE_DATA(error);
	}
	result = (int)(value[0] + 0.5);
	if (result != 1) {
		string error;
		error = "V_Flag from ImageFileInfo != 1 while loading the TIFF file at\"";
		error += path;
		error += "\" using Igor\r";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	result = FetchNumVar("V_numRows", &value[0], &value[1]);
	if (result == -1) {
		string error;
		error = "Error processing V_numRows from ImageFileInfo while loading the TIFF file at\"";
		error += path;
		error += "\" using Igor\r";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	x_size = (unsigned long)(value[0] + 0.5);
	
	result = FetchNumVar("V_numCols", &value[0], &value[1]);
	if (result == -1) {
		string error;
		error = "Error processing V_numCols from ImageFileInfo while loading the TIFF file at\"";
		error += path;
		error += "\" using Igor\r";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	y_size = (unsigned long)(value[0] + 0.5);
	
	result = FetchNumVar("V_numImages", &value[0], &value[1]);
	if (result == -1) {
		string error;
		error = "Error processing V_numImages from ImageFileInfo while loading the TIFF file at\"";
		error += path;
		error += "\" using Igor\r";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	total_number_of_images = (unsigned long)(value[0] + 0.5);
	
	return 0;
	
}
	
boost::shared_ptr<encap_gsl_matrix> ImageLoaderTIFF_Igor::get_nth_image(const unsigned long n) {
	
	boost::shared_ptr<encap_gsl_matrix> image_copy;
	int result;
	int number_of_images_loaded;
	double value[2];
	string IgorCommand;
	string part;
	stringstream sstream;
	long indices[MAX_DIMENSIONS];
	
	if (n >= total_number_of_images)
		throw IMAGE_INDEX_BEYOND_N_IMAGES();
	
	
	// is the requested image in the cache?
	// if it is then we return a copy
	for (unsigned long i = 0; i < image_cache_size; i++) {
		if (images_in_cache[i] == n) {
			
			image_copy = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(x_size, y_size));
			/*if (image_copy == NULL) {
				throw OUT_OF_MEMORY();
			}*/
			
			indices[2] = i;
			for (long j = 0; j < (long)y_size; j++) {
				for (long i = 0; i < (long)x_size; i++) {
					indices[0] = i;
					indices[1] = j;
					result = MDGetNumericWavePointValue(tiff_cache_wave, indices, value);
					if (result != 0) {
					//	gsl_matrix_free(image_copy);
						throw result;
					}
					image_copy->set(i, j, value[0]);
				}
			}
			
			return image_copy;
		}
	}
	
	// if we are here then the image wasn't in the cache
	// we will load a new sequence of images in the cache, starting with the one that is requested
	
	// if the cache wave already exists then we will delete it first
	if (tiff_cache_wave != NULL) {
		result = KillWave(tiff_cache_wave);
		if (result != 0) {
			throw result;
		}
		tiff_cache_wave = NULL;
		for (unsigned long i = 0; i < image_cache_size; i++) {
			images_in_cache[i] = -1;
		}
	}
	
	// now load the new set of images
	sstream << n;
	sstream << ' ';
	sstream << N_SIMULTANEOUS_IMAGE_LOADS;
	
	IgorCommand = "ImageLoad /Q /LR3D /N=tmp_TIFF_cache";
	
	sstream >> part;
	IgorCommand += " /S=";
	IgorCommand += part;
	IgorCommand += " /C=";
	sstream >> part;
	IgorCommand += part;
	IgorCommand += " \"";
	IgorCommand += native_Igor_FilePath;
	IgorCommand += "\"";
	
	// XOPNotice(IgorCommand.c_str());
	
	result = XOPSilentCommand(IgorCommand.c_str());
	if (result != 0) {
		XOPNotice("Error processing ImageLoad\r");
		throw result;
	}
	result = FetchNumVar("V_Flag", &value[0], &value[1]);
	if (result == -1) {
		string error;
		error = "Error processing V_Flag from ImageLoad while loading the TIFF file at\"";
		error += path;
		error += "\" using Igor\r";
		throw ERROR_READING_FILE_DATA(error);
	}
	result = (int)(value[0] + 0.5);
	if (result != 1) {
		string error;
		error = "V_Flag from ImageLoad != 1 while loading the TIFF file at\"";
		error += path;
		error += "\" using Igor\r";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	tiff_cache_wave = FetchWave("tmp_TIFF_cache");
	if (tiff_cache_wave == NULL) {
		string error;
		error = "Unable to retrieve the cache wave for the image at\"";
		error += path;
		error += "\"\r";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	result = FetchNumVar("V_numImages", &value[0], &value[1]);
	if (result == -1) {
		string error;
		error = "Error processing V_numImages from ImageLoad while loading the TIFF file at\"";
		error += path;
		error += "\" using Igor\r";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	number_of_images_loaded = (int)(value[0] + 0.5);
	
	// now update the array that holds the indices of the loaded images so that it matches
	for (int i = 0; i < number_of_images_loaded; i++) {
		images_in_cache[i] = n + i;
	}
	
	// now return a gsl matrix corresponding to the image that we want
	image_copy = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(x_size, y_size));
	/*if (image_copy == NULL) {
		throw OUT_OF_MEMORY();
	}*/
	
	// we want the first image if we just loaded a new set
	indices[2] = 0;
	for (long j = 0; j < (long)y_size; j++) {
		for (long i = 0; i < (long)x_size; i++) {
			indices[0] = i;
			indices[1] = j;
			result = MDGetNumericWavePointValue(tiff_cache_wave, indices, value);
			if (result != 0) {
				throw result;
			}
			image_copy->set(i, j, value[0]);
		}
	}
	
	return image_copy;
	
}
	

ImageLoaderTIFF::ImageLoaderTIFF(string rhs) {
	path = rhs;
	tiff_file = NULL;
	image_cache_size = N_SIMULTANEOUS_IMAGE_LOADS;
	
	tiff_file = TIFFOpen(path.c_str(), "r");
	if (tiff_file == NULL) {
		throw CANNOT_OPEN_FILE();
	}
	
	parse_header_information();
	
	image_cache.resize(image_cache_size, boost::shared_ptr<encap_gsl_matrix>());
	images_in_cache.resize(image_cache_size, (unsigned long)-1);
	
}

ImageLoaderTIFF::ImageLoaderTIFF(string rhs, unsigned long image_cache_size_rhs) {
	path = rhs;
	image_cache_size = image_cache_size_rhs;
	if (image_cache_size > N_SIMULTANEOUS_IMAGE_LOADS) {
		image_cache.reserve(image_cache_size);
	}
	tiff_file = NULL;
	
	tiff_file = TIFFOpen(path.c_str(), "r");
	if (tiff_file == NULL) {
		throw CANNOT_OPEN_FILE();
	}
	
	parse_header_information();
	
	image_cache.resize(image_cache_size, boost::shared_ptr<encap_gsl_matrix>());
	images_in_cache.resize(image_cache_size, (unsigned long)-1);
	
}


ImageLoaderTIFF::~ImageLoaderTIFF() {
	
	if (tiff_file != NULL) {
		TIFFClose(tiff_file);
	}
	
}


int ImageLoaderTIFF::parse_header_information() {
	int result;
	uint16_t result_uint16;
	uint32_t result_uint32;
	unsigned long index = 0;
	
	
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
			sampleFormat = 1;
			break;
		case 3:
			sampleFormat = 3;
			break;
		default:
			string error;
			error = "The SampleFormat of the image at\"";
			error += path;
			error += "\" is unknown\r";
			throw ERROR_READING_FILE_DATA(error);
			break;
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
	
	x_size = (unsigned long)(result_uint32);
	
	// what is the y size?
	result = TIFFGetField(tiff_file, TIFFTAG_IMAGELENGTH, &result_uint32);
	if (result != 1) {
		string error;
		error = "The image at\"";
		error += path;
		error += "\" does not specify a height\r";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	y_size = (unsigned long)(result_uint32);
	
	
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
	
	return 0;
	
}

boost::shared_ptr<encap_gsl_matrix> ImageLoaderTIFF::get_nth_image(const unsigned long n) {
	
	if (n >= total_number_of_images)
		throw IMAGE_INDEX_BEYOND_N_IMAGES();
	
	// is there a file open?
	if (tiff_file == NULL)
		throw GET_NTH_IMAGE_FILE_NOT_OPEN();
	
	
	// is the requested image in the cache?
	// if it is then we return a copy
	for (unsigned long i = 0; i < image_cache_size; i++) {
		if (images_in_cache[i] == n) {
			return image_cache[i];
		}
	}
	
	// if we are here then the image wasn't in the cache
	// we will load a new sequence of images in the cache, starting with the one that is requested
	
	// first we free the images that are already in the cache
	// this is safe because we only pass copies to other functions
	for (unsigned long i = 0; i < image_cache_size; i++) {
		images_in_cache[i] = -1;
	}
	
	// now load the new set of images
	
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
	boost::shared_ptr<encap_gsl_matrix> new_image;
	int result;
	
	single_scanline_buffer = (char *)_TIFFmalloc(TIFFScanlineSize(tiff_file));
	if (single_scanline_buffer == NULL) {
		string error;
		error = "unable to allocate scanline buffer in ImageLoaderTIFF::get_nth_image()\r";
		throw OUT_OF_MEMORY(error);
	}
	
	for (unsigned long i = n; i < n + image_cache_size; i++) {
		
		if (i >= total_number_of_images) {
			continue;	// there are no more images to load into this cache location
			// by freeing the locations above we have already made sure that the arrays are set to a NULL pointer and -1
		}
		
		try {
			new_image = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(x_size, y_size));
		}
		catch (std::bad_alloc) {
			for (unsigned long j = 0; j < i - n; j++) {
				images_in_cache[j] = -1;
			}
			_TIFFfree(single_scanline_buffer);
			string error;
			error = "unable to allocate new_image in ImageLoaderTIFF::get_nth_image()\r";
			throw OUT_OF_MEMORY(error);
		}
		
		result = TIFFSetDirectory(tiff_file, directoryIndices.at(i));
		if (result != 1) {
			string error;
			error = "Unable to set the directory to '0' for the image at\"";
			error += path;
			error += "\"\r";
			throw ERROR_READING_FILE_DATA(error);
		}
		
		for (unsigned long j = 0; j < y_size; ++j) {
			result = TIFFReadScanline(tiff_file, single_scanline_buffer, j, 0);	// sample is ignored
			if (result != 1) {
				string error;
				error = "Unable to read a scanline from the image at\"";
				error += path;
				error += "\"\r";
				throw ERROR_READING_FILE_DATA(error);
			}
			
			switch (sampleFormat) {	// handle the different possibilities (floating, integer) and variable sizes
				case 1:	// the data is present as an integer
					switch (bitsPerPixel) {	// handle the different sampling sizes
						case 4:
							scanline_pointer = single_scanline_buffer;
							for (unsigned long k = 0; k < x_size; ++k) {
								if ((k % 2) == 0) {	// this is an even pixel, we use only the first 4 bits
									current_uint16 = 0x0000000F & (*scanline_pointer);
									new_image->set(k, j, (double)current_uint16);
								} else {	// this is an odd pixel, use the last 4 bits and increment the scanline_pointer
									current_uint16 = 0x000000F0 & (*scanline_pointer);
									new_image->set(k, j, (double)current_uint16);
									scanline_pointer += 1;
								}
							}
							break;
						case 8:
							scanline_pointer = single_scanline_buffer;
							for (unsigned long k = 0; k < x_size; ++k) {
								current_uint16 = (uint16_t)(*scanline_pointer);
								new_image->set(k, j, (double)current_uint16);
								scanline_pointer += 1;
							}
							break;
						case 16:
							uint16Ptr = (uint16_t*)single_scanline_buffer;
							for (unsigned long k = 0; k < x_size; ++k) {
								current_uint16 = (*uint16Ptr);
								new_image->set(k, j, (double)current_uint16);
								uint16Ptr += 1;
							}
							break;
						case 32:
							uint32Ptr = (uint32_t*)single_scanline_buffer;
							for (unsigned long k = 0; k < x_size; ++k) {
								current_uint32 = (*uint32Ptr);
								new_image->set(k, j, (double)current_uint32);
								uint32Ptr += 1;
							}
							break;
						default:
							_TIFFfree(single_scanline_buffer);
							string error;
							error = "Invalid integer data size for the image at\"";
							error += path;
							error += "\"\r";
							throw ERROR_READING_FILE_DATA(error);
							break;
					}
					break;
				case 3:	// the data are stored as floating point values
					switch (bitsPerPixel) {
						case 32:
							floatPtr = (float *)single_scanline_buffer;
							for (unsigned long k = 0; k < x_size; ++k) {
								current_float = *floatPtr;
								new_image->set(k, j, (double)current_float);
								floatPtr += 1;
							}
							break;
						case 64:
							doublePtr = (double *)single_scanline_buffer;
							for (unsigned long k = 0; k < x_size; ++k) {
								current_double = *doublePtr;
								new_image->set(k, j, current_double);
								floatPtr += 1;
							}
							break;
						default:
							_TIFFfree(single_scanline_buffer);
							_TIFFfree(single_scanline_buffer);
							string error;
							error = "Invalid floating point data size for the image at\"";
							error += path;
							error += "\"\r";
							throw ERROR_READING_FILE_DATA(error);
							break;
					}
					break;
				default:
					_TIFFfree(single_scanline_buffer);
					_TIFFfree(single_scanline_buffer);
					string error;
					error = "Unknown SampleFormat for the image at\"";
					error += path;
					error += "\"\r";
					throw ERROR_READING_FILE_DATA(error);
					break;
			}
		}
		
		// store the image at the correct location in the cache
		image_cache[i - n] = new_image;
		images_in_cache[i - n] = i;
	}
	
	_TIFFfree(single_scanline_buffer);
	result = TIFFSetDirectory(tiff_file, 0);
	if (result != 1) {
		_TIFFfree(single_scanline_buffer);
		string error;
		error = "Invalid to set the directory to '0' for the image at\"";
		error += path;
		error += "\"\r";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	return image_cache[0];
	
}



ImageLoaderIgor::ImageLoaderIgor(string waveName) {
	// try to get images from an Igor wave
	// the string that is passed in has a leading slash added by the code that converts between Macintosh and Windows paths
	// so we need to correct for that
	#ifdef _MACINTOSH_
		waveName.erase(0, 1);
	#endif
	
	igor_data_wave = FetchWave(waveName.c_str());
	if (igor_data_wave == NULL) {
		throw NOWAV;
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
	
	x_size = (unsigned long)DimensionSizes[0];
	y_size = (unsigned long)DimensionSizes[1];
	total_number_of_images = (unsigned long)DimensionSizes[2];
	
	// special case: if the wave contains only a single image then it is usually two-dimensional, that is, DimensionSizes[2] == 0
	// in that case total_number_of_images is still one
	if (DimensionSizes[2] == 0) {
		total_number_of_images = 1;
	}
}

boost::shared_ptr<encap_gsl_matrix> ImageLoaderIgor::get_nth_image(const unsigned long n) {
	boost::shared_ptr<encap_gsl_matrix> image;
	double value[2];
	long indices[3];
	int result;
	
	if (n >= total_number_of_images)
		throw IMAGE_INDEX_BEYOND_N_IMAGES();
	
	image = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(x_size, y_size));
	
	indices[2] = n;
	
	for (unsigned long j = 0; j < y_size; j++) {
		for (unsigned long i = 0; i < x_size; i++) {
			indices[0] = (long)i;
			indices[1] = (long)j;
			
			result = MDGetNumericWavePointValue(igor_data_wave, indices, value);
			if (result != 0) {
				throw result;
			}
			
			image->set(i, j, value[0]);
		}
	}
	
	return image;
}



CCDImagesProcessorAverageSubtraction::CCDImagesProcessorAverageSubtraction(ImageLoader *i_loader, OutputWriter *o_writer) {
	image_loader = i_loader;
	output_writer = o_writer;
	
	total_number_of_images = image_loader->get_total_number_of_images();
	x_size = image_loader->get_x_size();
	y_size = image_loader->get_y_size();
	
	n_frames_averaging = 0;
}


void CCDImagesProcessorAverageSubtraction::set_n_parameter(double n) {
	// check is we are indeed averaging over an odd number
	// if not then we throw an error
	// the exception is that a value of '0' means that we have to average over the entire sequence
	
	// by convention the calling function checks that n is positive
	n_frames_averaging = (unsigned long)(n + 0.5);
	
	if (((n_frames_averaging % 2) != 1) && (n_frames_averaging != 0)) {
		throw NUMBER_OF_AVERAGING_FRAMES_SHOULD_BE_ODD();
	}
	
	if (n_frames_averaging > total_number_of_images) {
		throw NUMBER_OF_AVERAGING_LARGER_THAN_N_FRAMES_IN_FILE();
	}
}


int CCDImagesProcessorAverageSubtraction::convert_images() {
	if (n_frames_averaging == 0) {	// we want to average over the entire trace
		subtract_average_of_entire_trace();
	} else {
		subtract_partial_average();
	}
	
	return 0;
}

void CCDImagesProcessorAverageSubtraction::subtract_average_of_entire_trace() {
	unsigned long n;
	boost::shared_ptr<encap_gsl_matrix> average_image;
	boost::shared_ptr<encap_gsl_matrix> loaded_image;
	boost::shared_ptr<encap_gsl_matrix> subtracted_image;
	double current_double;
	double value;
	
	// we pass through the images two times:
	// the first pass calculates the average,
	// the second pass subtracts it from the image
	// fortunately out intermediate format uses doubles to store the data!
	
	average_image = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(x_size, y_size));
	
	average_image->set_all(0);	// zero the matrix
	
	for (n = 0; n < total_number_of_images; n++) {
		loaded_image = image_loader->get_nth_image(n);
		
		// average_image->add(*loaded_image);
		for (unsigned long l = 0; l < y_size; ++l) {
			for (unsigned long k = 0; k < x_size; ++k) {
				value = average_image->get(k, l);
				value += loaded_image->get(k, l);
				average_image->set(k, l, value);
			}
		}
	}
	
	// now divide each point so that we get the average
	// gsl_matrix_scale(average_image, (1.0 / (double)n_frames_averaging));
	
	// the approach using the gsl functions seems off so we use a different one instead
	
	for (unsigned long i = 0; i < x_size; i++) {
		for (unsigned long j = 0; j < y_size; j++) {
			current_double = average_image->get(i, j);
			current_double /= (double)total_number_of_images;
			average_image->set(i, j, current_double);
		}
	}
	
	// now subtract the average for each frame
	for (n = 0; n < total_number_of_images; n++) {
		loaded_image = image_loader->get_nth_image(n);
		
		// loaded_image->sub(*average_image);
		for (unsigned long k = 0; k < x_size; ++k) {
			for (unsigned long l = 0; l < y_size; ++l) {
				value = loaded_image->get(k, l);
				value -= average_image->get(k, l);
				loaded_image->set(k, l, value);
			}
		}
		
		subtracted_image = loaded_image;
		
		output_writer->write_image(subtracted_image);	// the output writer will take care of freeing the memory
	}
}



void CCDImagesProcessorAverageSubtraction::subtract_partial_average() {
	boost::shared_ptr<encap_gsl_matrix> current_image;
	boost::shared_ptr<encap_gsl_matrix> average_image;
	boost::shared_ptr<encap_gsl_matrix> subtracted_image;
	// gsl_matrix *averaging_buffer[n_frames_averaging];	// strictly speaking this isn't legal C++ code because declarations on the stack should be
															// a const size
															// g++ accepts this, but the Microsoft compiler throws an error on this
															// so we will convert it to an assignment on the heap instead
	vector<boost::shared_ptr<encap_gsl_matrix> > averaging_buffer;
	averaging_buffer.resize(n_frames_averaging, boost::shared_ptr<encap_gsl_matrix> ());
	
	long average_starting_index, average_ending_index;
	unsigned long cache_loading_offset = 0;
	
	double current_double;
	double value;
	
	average_image = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(x_size, y_size));
	
	/*// the buffer for the averaging will be allocated by the get_nth_image() routine
	for (unsigned long i = 0; i < n_frames_averaging; i++) {
		averaging_buffer[i] = NULL;
	}*/
	
	// we normally try to subtract the average of the frames surrounding the frame that we are interested in
	// this is why we demand that n_frames_averaging is odd
	
	// however, if we have a frame at the beginning or end of the image stack then we cannot do this
	// instead we construct an average of n_frames_averaging as close as possible to the frame we are interested in
	
	// loop over all the images
	for (unsigned long n = 0; n < total_number_of_images; n++) {
		average_starting_index = n - floor((double)total_number_of_images / 2.0);
		average_ending_index = n + floor((double)total_number_of_images / 2.0);
		
		// SPECIAL CASE: check if we are too close to the START of the image stack to calculate a normal average
		
		if (average_starting_index <= 0) {
			// all the points that fulfill this condition calculate their average using the same set of frames
			// so we only load these frames once, for the first image
			if (n == 0) {
				for (unsigned long i = 0; i < n_frames_averaging; i++) {
					
					averaging_buffer.at(i) = image_loader->get_nth_image(i);
					
				}
				
				// now calculate the average once
				average_image->set_all(0);
				for (unsigned long i = 0; i < n_frames_averaging; i++) {
					for (unsigned long k = 0; k < x_size; ++k) {
						for (unsigned long l = 0; l < y_size; ++l) {
							value = average_image->get(k, l);
							value += averaging_buffer.at(i)->get(k, l);
							average_image->set(k, l, value);
						}
					}
				}
				// gsl_matrix_scale(average_image, (1.0 / (double)n_frames_averaging));
				// the calculation using the gsl seems to be off so we use a direct one instead
				
				for (unsigned long i = 0; i < x_size; i++) {
					for (unsigned long j = 0; j < y_size; j++) {
						current_double = average_image->get(i, j);
						current_double /= n_frames_averaging;
						average_image->set(i, j, current_double);
					}
				}
				
			}
			// the frame that we are interested in (index n) can now also be found at
			// averaging_buffer[n]
			current_image = averaging_buffer.at(n);
			
			// now allocate the subtracted image
			subtracted_image = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(x_size, y_size));
			
			// subtract the average from the current frame
			// subtracted_image->copy(*current_image);
			for (unsigned long j = 0; j < y_size; ++j) {
				for (unsigned long i = 0; i < x_size; ++i) {
					subtracted_image->set(i, j, current_image->get(i, j));
				}
			}
			// subtracted_image->sub(*average_image);
			for (unsigned long k = 0; k < x_size; ++k) {
				for (unsigned long l = 0; l < y_size; ++l) {
					value = subtracted_image->get(k, l);
					value -= average_image->get(k, l);
					subtracted_image->set(k, l, value);
				}
			}
			
			// store the subtracted image
			output_writer->write_image(subtracted_image);	// output_writer will take care of freeing the memory
			
			continue;	// don't go through any of the other stuff in the main for loop
						// while we are too close to the edge of the image stack
		}
		
		
		// SPECIAL CASE: check if we are too close to the END of the image stack to calculate a normal average
		if ((unsigned long)average_ending_index >= total_number_of_images) {
			
			// we don't need to update the image buffer anymore, all the images we need to calculate the average are stored in memory
			// also, we don't need to calculate this average as it will have been calculated by the last image
			// that could calculate a normal average
			
			current_image = averaging_buffer[(unsigned long)floor((double)(n_frames_averaging) / 2.0) + cache_loading_offset];
			
			// now allocate the subtracted image
			subtracted_image = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(x_size, y_size));
			
			// subtracted_image->copy(*current_image);
			for (unsigned long j = 0; j < y_size; ++j) {
				for (unsigned long i = 0; i < x_size; ++i) {
					subtracted_image->set(i, j, current_image->get(i, j));
				}
			}
			
			
			// subtracted_image->sub(*average_image);
			for (unsigned long k = 0; k < x_size; ++k) {
				for (unsigned long l = 0; l < y_size; ++l) {
					value = subtracted_image->get(k, l);
					value -= average_image->get(k, l);
					subtracted_image->set(k, l, value);
				}
			}
			
			output_writer->write_image(subtracted_image);
			
			cache_loading_offset++;
			continue;
		}
				
		
		// NORMAL CASE: we are somewhere in the middle of the image stack
		
		averaging_buffer.erase(averaging_buffer.begin());
		// we shift the images to the right by one position
		for (unsigned long i = 0; i < (n_frames_averaging - 1); i++) {
			
			averaging_buffer.at(i) = averaging_buffer.at(i + 1);
		}
		
		averaging_buffer.at(n_frames_averaging - 1) = image_loader->get_nth_image(n + floor((double)(n_frames_averaging) / 2.0));
		
		// if we are here then we can calculate a new average
		average_image->set_all(0);
		
		for (unsigned long i = 0; i < n_frames_averaging; i++) {
			for (unsigned long i = 0; i < x_size; ++i) {
				for (unsigned long j = 0; j < y_size; ++j) {
					value = average_image->get(i, j);
					value += current_image->get(i, j);
					average_image->set(i, j, value);
				}
			}
		}
		
		// gsl_matrix_scale(average_image, (1.0 / (double)n_frames_averaging));
		// the calculation using the gsl seems to be off so we use a direct one instead
		
		for (unsigned long i = 0; i < x_size; i++) {
			for (unsigned long j = 0; j < y_size; j++) {
				current_double = average_image->get(i, j);
				current_double /= n_frames_averaging;
				average_image->set(i, j, current_double);
			}
		}
		
		// the image that we are interested in is at the center of the imaging buffer
		current_image = averaging_buffer.at((unsigned long)floor((double)(n_frames_averaging) / 2.0));
		
		// now allocate the subtracted image
		subtracted_image = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(x_size, y_size));
		
		// subtracted_image->copy(*current_image);
		for (unsigned long j = 0; j < y_size; ++j) {
			for (unsigned long i = 0; i < x_size; ++i) {
				subtracted_image->set(i, j, current_image->get(i, j));
			}
		}
		
		// subtracted_image->sub(*average_image);
		for (unsigned long k = 0; k < x_size; ++k) {
			for (unsigned long l = 0; l < y_size; ++l) {
				value = subtracted_image->get(k, l);
				value -= average_image->get(k, l);
				subtracted_image->set(k, l, value);
			}
		}
		
		output_writer->write_image(subtracted_image);
	}
	
	// phew! we're done
	// now we just have to make sure that we don't leave a mess behind when we close
	
}
	
	
CCDImagesProcessorDifferenceImage::CCDImagesProcessorDifferenceImage(ImageLoader *i_loader, OutputWriter *o_writer) {
	image_loader = i_loader;
	output_writer = o_writer;
	
	total_number_of_images = image_loader->get_total_number_of_images();
	x_size = image_loader->get_x_size();
	y_size = image_loader->get_y_size();
}

int CCDImagesProcessorDifferenceImage::convert_images() {
	
	boost::shared_ptr<encap_gsl_matrix> current_image;
	boost::shared_ptr<encap_gsl_matrix> next_image;
	
	// we start by loading the first image in next_image
	// this is required so the loop that follows can start properly
	
	next_image = image_loader->get_nth_image(0);
	
	
	for (unsigned long n = 0; n < (total_number_of_images - 1); n++) {
		
		// the previous image for this run of the loop is the image that was previously in current_image
		// so we have to shift it down
		current_image = next_image;
		
		next_image = image_loader->get_nth_image(n + 1);
		
		// now do the actual subtraction
		gsl_matrix_sub(current_image->get_ptr(), next_image->get_ptr());
		
		// current_image now contains the subtracted image, we should write it to disk
		output_writer->write_image(current_image);
		
		// the output_writer also takes care of freeing current_image
	}
	
	// before we exit the function we need to free next_image
	// gsl_matrix_free(next_image);
	
	return 0;
}
		

CCDImagesProcessorConvertToSimpleFileFormat::CCDImagesProcessorConvertToSimpleFileFormat(ImageLoader *i_loader, OutputWriter *o_writer) {
	image_loader = i_loader;
	output_writer = o_writer;
	
	total_number_of_images = image_loader->get_total_number_of_images();
	x_size = image_loader->get_x_size();
	y_size = image_loader->get_y_size();
}

int CCDImagesProcessorConvertToSimpleFileFormat::convert_images() {
	
	boost::shared_ptr<encap_gsl_matrix> current_image;
	
	for (unsigned long n = 0; n < total_number_of_images; ++n) {
		current_image = image_loader->get_nth_image(n);
		output_writer->write_image(current_image);
	}
	return 0;
}


boost::shared_ptr<encap_gsl_matrix> ParticleFinder_radius::findPositions(boost::shared_ptr<encap_gsl_matrix> image, boost::shared_ptr<encap_gsl_matrix_uchar> threshold_image) {
	vector<position> positions;
	position position;
	// we store the pixels above the treshold as a vector containing x,y,intensity
	unsigned long x_size = image->get_x_size(), y_size = image->get_y_size();
	unsigned long number_of_positions;
	double current_intensity, previous_intensity, current_x, current_y;
	double distance_squared;
	double radius_squared = (double)(radius * radius);
	int skip;
	double x, y;
	boost::shared_ptr<encap_gsl_matrix> output_positions;
	
	// we run over all the points in the image to see if they are above the treshold
	for (unsigned long j = minDistanceFromEdge; j < y_size - minDistanceFromEdge; j++) {
		for (unsigned long i = minDistanceFromEdge; i < x_size - minDistanceFromEdge; i++) {
			
			if (threshold_image->get(i, j) < 128)
				continue;	// we don't care about this point, it's not included in the thresholded image
			
			current_intensity = image->get(i, j);
			
			current_x = (double)i;
			current_y = (double)j;
			
			// is this point too close to the edge of the image?
			if ((current_x < minDistanceFromEdge) || (current_x > (x_size - minDistanceFromEdge)))
				continue;
			if ((current_y < minDistanceFromEdge) || (current_y > (y_size - minDistanceFromEdge)))
				continue;
			
			// if we are still here then we need to take a closer look at this point
			// check if the current point overlaps with a previous point
			skip = 0;
			number_of_positions = positions.size();
			for (unsigned long k = 0; k < number_of_positions; k++) {
				x = positions[k].get_x();
				y = positions[k].get_y();
				distance_squared = (current_x - x) * (current_x - x) + (current_y - y) * (current_y - y);
				
				if (distance_squared < radius_squared) {
					// we need to skip one of the two pixels that we are comparing
					// we will keep the pixel with the largest intensity
					previous_intensity = positions[k].get_intensity();
					if (current_intensity > previous_intensity) {
						positions[k].set_intensity(current_intensity);
						positions[k].set_x(current_x);
						positions[k].set_y(current_y);
					}
					skip = 1;
					break;
				}
			}
			
			if (skip == 0) {	// we should store this point
				position.set_intensity(current_intensity);
				position.set_x(current_x);
				position.set_y(current_y);
				positions.push_back(position);
			}
		}
	}
	
	// now we need to store the data in the standard matrix format
	// this means that the columns are oriented as intensity, x, y
	
	number_of_positions = positions.size();
	
	if (number_of_positions == 0) {	// no positions were found
		return output_positions;	// is equal to NULL due to its initialization as a shared_ptr
	}
	
	output_positions = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(number_of_positions, 3));
	
	for (unsigned long i = 0; i < number_of_positions; i++) {
		output_positions->set(i, 0, positions[i].get_intensity());
		output_positions->set(i, 1, positions[i].get_x());
		output_positions->set(i, 2, positions[i].get_y());
	}
	
	return output_positions;
}


boost::shared_ptr<encap_gsl_matrix> ParticleFinder_adjacent4::findPositions(boost::shared_ptr<encap_gsl_matrix> image, boost::shared_ptr<encap_gsl_matrix_uchar> threshold_image) {
	
	boost::shared_ptr<encap_gsl_matrix> output_positions;
	list<position> positionsInCurrentParticleList;
	vector<position> positionsInCurrentParticle;
	position currentPosition;
	vector<position> particles;
	boost::shared_ptr<encap_gsl_matrix_long> mapped_image;	// keeps track of which pixels have already been mapped to a particle
	unsigned long x_size = image->get_x_size();
	unsigned long y_size = image->get_y_size();
	unsigned long x, y;
	long particleIndex = 0;
	double average_x, average_y;
	double maxIntensity;
	
	mapped_image = boost::shared_ptr<encap_gsl_matrix_long>(new encap_gsl_matrix_long(x_size, y_size));
	
	mapped_image->set_all(-1);
	
	for (unsigned long j = minDistanceFromEdge; j < y_size - minDistanceFromEdge; ++j) {
		for (unsigned long i = minDistanceFromEdge; i < x_size - minDistanceFromEdge; ++i) {	// loop over the entire image
			
			if (threshold_image->get(i, j) < 128) {
				continue;
			}
			
			if (mapped_image->get(i, j) != -1) {	// this point is already assigned to another particle
				continue;
			}
			
			// if we are still here then we found an active pixel that is not yet assigned to a new particle
			// we create a new particle at this position
			positionsInCurrentParticle.clear();
			positionsInCurrentParticleList.clear();
			
			mapped_image->set(i, j, particleIndex);
			
			// store this position
			currentPosition.set_x((double)i);
			currentPosition.set_y((double)j);
			currentPosition.set_intensity(image->get(i, j));
			
			positionsInCurrentParticleList.push_back(currentPosition);
			
			// growParticle(currentPosition, i, j, threshold_image, mapped_image);
			
			
			while (positionsInCurrentParticleList.size() > 0) {
				currentPosition = positionsInCurrentParticleList.front();
				x = (unsigned long)(currentPosition.get_x() + 0.5);
				y = (unsigned long)(currentPosition.get_y() + 0.5);
				
				growParticle(currentPosition, positionsInCurrentParticleList, image, threshold_image, mapped_image);
				// growParticle will update the list with new positions
				// this position has been checked, so we don't need to include it in future searches
				positionsInCurrentParticleList.pop_front();
				positionsInCurrentParticle.push_back(currentPosition);
			}
			
			// we have found all positions in this particle
			++particleIndex;
			
			// store the output positions
			// first calculate an average in x and y
			// and get an estimate for the intensity of the particle
			maxIntensity = 0;
			average_x = 0;
			average_y = 0;
			for (unsigned long k = 0; k < positionsInCurrentParticle.size(); ++k) {
				average_x += positionsInCurrentParticle[k].get_x();
				average_y += positionsInCurrentParticle[k].get_y();
				if (positionsInCurrentParticle[k].get_intensity() > maxIntensity) {
					maxIntensity = positionsInCurrentParticle[k].get_intensity();
				}
			}
			average_x /= positionsInCurrentParticle.size();
			average_y /= positionsInCurrentParticle.size();
			currentPosition.set_intensity(maxIntensity);
			currentPosition.set_x(average_x);
			currentPosition.set_y(average_y);
			
			// we store the particle only if it is not too close to the edge of the frame
			if ((average_x < minDistanceFromEdge) || (average_x > ((double)x_size - minDistanceFromEdge)))
				continue;
			if ((average_y < minDistanceFromEdge) || (average_y > ((double)y_size - minDistanceFromEdge)))
				continue;
			
			particles.push_back(currentPosition);
		}
	}
	
	// if we have found some particles then return them, else return a null pointer
	if (particles.size() == 0) {
		return output_positions;
	}
	
	output_positions = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(particles.size(), 3));
	
	// now copy the output to a gsl matrix
	for (unsigned long k = 0; k < particles.size(); ++k) {
		output_positions->set(k, 0, particles[k].get_intensity());
		output_positions->set(k, 1, particles[k].get_x());
		output_positions->set(k, 2, particles[k].get_y());
	}
	
	return output_positions;
	
}

void ParticleFinder_adjacent4::growParticle(position centerPosition, list<position> &positionsInCurrentParticle, boost::shared_ptr<encap_gsl_matrix> image, boost::shared_ptr<encap_gsl_matrix_uchar> threshold_image, boost::shared_ptr<encap_gsl_matrix_long> mapped_image) {
	// the pixel at position (x,y) belongs to a particle
	// do the surrounding pixels belong to the same particle?
	
	// the function checks which of the pixels surrounding pos are active
	// if one or more of these is active, then it checks if they are already assigned to the particle by checking mapped_image
	// if they are not known then they are added to to the list with positions of the current particle
	// and also added to mapped_image
	
	unsigned long x_size = image->get_x_size();
	unsigned long y_size = image->get_y_size();
	position currentPosition;
	
	unsigned long x = (unsigned long)(centerPosition.get_x() + 0.5);
	unsigned long y = (unsigned long)(centerPosition.get_y() + 0.5);
	
	long particleIndex = mapped_image->get(x, y);
	
	if ((x < minDistanceFromEdge) || (x > x_size - minDistanceFromEdge - 1))
		return;
	if ((y < minDistanceFromEdge) || (y > y_size - minDistanceFromEdge - 1))
		return;
	
	// is the pixel to the left of the current one active?
	if (x > 0) {
		if (threshold_image->get(x - 1, y) > 128) {
			// it's active
			// did we already include this pixel?
			if (mapped_image->get(x - 1, y) == -1) {
				mapped_image->set(x - 1, y, particleIndex);
				// add the point to the vector
				currentPosition.set_x((double)x - 1);
				currentPosition.set_y((double)y);
				currentPosition.set_intensity(image->get(x - 1, y));
				positionsInCurrentParticle.push_back(currentPosition);
				// now call this function itself (recursion) on the new pixel
				// growParticle(positions, x - 1, y, threshold_image, mapped_image);
			}
		}
	}
	// is the pixel to the right of the current one active?
	if (x < x_size - 1) {
		if (threshold_image->get(x + 1, y) > 128) {
			// it's active
			// did we already include this pixel?
			if (mapped_image->get(x + 1, y) == -1) {
				mapped_image->set(x + 1, y, particleIndex);
				// add the point to the vector
				currentPosition.set_x((double)x + 1);
				currentPosition.set_y((double)y);
				currentPosition.set_intensity(image->get(x + 1, y));
				positionsInCurrentParticle.push_back(currentPosition);
				// now call this function itself (recursion) on the new pixel
				// growParticle(positions, x + 1, y, threshold_image, mapped_image);
			}
		}
	}
	// is the pixel above the current one active?
	if (y > 0) {
		if (threshold_image->get(x, y - 1) > 128) {
			// it's active
			// did we already include this pixel?
			if (mapped_image->get(x, y - 1) == -1) {
				mapped_image->set(x, y - 1, particleIndex);
				// add the point to the vector
				currentPosition.set_x((double)x);
				currentPosition.set_y((double)y - 1);
				currentPosition.set_intensity(image->get(x, y - 1));
				positionsInCurrentParticle.push_back(currentPosition);
				// now call this function itself (recursion) on the new pixel
				// growParticle(positions, x, y - 1, threshold_image, mapped_image);
			}
		}
	}
	// is the pixel below the current one active?
	if (y < y_size - 1) {
		if (threshold_image->get(x, y + 1) > 128) {
			// it's active
			// did we already include this pixel?
			if (mapped_image->get(x, y + 1) == -1) {
				mapped_image->set(x, y + 1, particleIndex);
				// add the point to the vector
				currentPosition.set_x((double)x);
				currentPosition.set_y((double)y + 1);
				currentPosition.set_intensity(image->get(x, y + 1));
				positionsInCurrentParticle.push_back(currentPosition);
				// now call this function itself (recursion) on the new pixel
				// growParticle(positions, x, y + 1, threshold_image, mapped_image);
			}
		}
	}
}


boost::shared_ptr<encap_gsl_matrix> ParticleFinder_adjacent8::findPositions(boost::shared_ptr<encap_gsl_matrix> image, boost::shared_ptr<encap_gsl_matrix_uchar> threshold_image) {
	
	boost::shared_ptr<encap_gsl_matrix> output_positions;
	list<position> positionsInCurrentParticleList;
	vector<position> positionsInCurrentParticle;
	position currentPosition;
	vector<position> particles;
	boost::shared_ptr<encap_gsl_matrix_long> mapped_image;	// keeps track of which pixels have already been mapped to a particle
	unsigned long x_size = image->get_x_size();
	unsigned long y_size = image->get_y_size();
	unsigned long x, y;
	long particleIndex = 0;
	double average_x, average_y;
	double maxIntensity;
	
	mapped_image = boost::shared_ptr<encap_gsl_matrix_long>(new encap_gsl_matrix_long(x_size, y_size));
	
	mapped_image->set_all(-1);
	
	for (unsigned long j = minDistanceFromEdge; j < y_size - minDistanceFromEdge; ++j) {
		for (unsigned long i = minDistanceFromEdge; i < x_size - minDistanceFromEdge; ++i) {	// loop over the entire image
			
			if (threshold_image->get(i, j) < 128) {
				continue;
			}
			
			if (mapped_image->get(i, j) != -1) {	// this point is already assigned to another particle
				continue;
			}
			
			// if we are still here then we found an active pixel that is not yet assigned to a new particle
			// we create a new particle at this position
			positionsInCurrentParticle.clear();
			positionsInCurrentParticleList.clear();
			
			mapped_image->set(i, j, particleIndex);
			
			// store this position
			currentPosition.set_x((double)i);
			currentPosition.set_y((double)j);
			currentPosition.set_intensity(image->get(i, j));
			
			positionsInCurrentParticleList.push_back(currentPosition);
			
			// growParticle(currentPosition, i, j, threshold_image, mapped_image);
			
			
			while (positionsInCurrentParticleList.size() > 0) {
				currentPosition = positionsInCurrentParticleList.front();
				x = (unsigned long)(currentPosition.get_x() + 0.5);
				y = (unsigned long)(currentPosition.get_y() + 0.5);
				
				growParticle(currentPosition, positionsInCurrentParticleList, image, threshold_image, mapped_image);
				// growParticle will update the list with new positions
				// this position has been checked, so we don't need to include it in future searches
				positionsInCurrentParticleList.pop_front();
				positionsInCurrentParticle.push_back(currentPosition);
			}
			
			
			// we're done with this particle, time for the next one
			++particleIndex;
			
			// store the output positions
			// first calculate an average in x and y
			// and get an estimate for the intensity of the particle
			maxIntensity = 0;
			average_x = 0;
			average_y = 0;
			for (unsigned long k = 0; k < positionsInCurrentParticle.size(); ++k) {
				average_x += positionsInCurrentParticle[k].get_x();
				average_y += positionsInCurrentParticle[k].get_y();
				if (positionsInCurrentParticle[k].get_intensity() > maxIntensity) {
					maxIntensity = positionsInCurrentParticle[k].get_intensity();
				}
			}
			average_x /= positionsInCurrentParticle.size();
			average_y /= positionsInCurrentParticle.size();
			currentPosition.set_intensity(maxIntensity);
			currentPosition.set_x(average_x);
			currentPosition.set_y(average_y);
			
			// we store the particle only if it is not too close to the edge of the frame
			if ((average_x < minDistanceFromEdge) || (average_x > ((double)x_size - minDistanceFromEdge)))
				continue;
			if ((average_y < minDistanceFromEdge) || (average_y > ((double)y_size - minDistanceFromEdge)))
				continue;
			
			particles.push_back(currentPosition);
		}
	}
	
	// if we have found some particles then return them, else return a null pointer
	if (particles.size() == 0) {
		return output_positions;
	}
	
	output_positions = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(particles.size(), 3));
	
	// now copy the output to a gsl matrix
	for (unsigned long k = 0; k < particles.size(); ++k) {
		output_positions->set(k, 0, particles[k].get_intensity());
		output_positions->set(k, 1, particles[k].get_x());
		output_positions->set(k, 2, particles[k].get_y());
	}
	
	return output_positions;
	
}

void ParticleFinder_adjacent8::growParticle(position centerPosition, list<position> &positionsInCurrentParticle, boost::shared_ptr<encap_gsl_matrix> image, boost::shared_ptr<encap_gsl_matrix_uchar> threshold_image, boost::shared_ptr<encap_gsl_matrix_long> mapped_image) {
	// the pixel at position (x,y) belongs to a particle
	// do the surrounding pixels belong to the same particle?
	
	// the function checks which of the pixels surrounding pos are active
	// if one or more of these is active, then it checks if they are already assigned to the particle by checking mapped_image
	// if they are not known then they are added to to the list with positions of the current particle
	// and also added to mapped_image
	
	unsigned long x_size = image->get_x_size();
	unsigned long y_size = image->get_y_size();
	position currentPosition;
	
	unsigned long x = (unsigned long)(centerPosition.get_x() + 0.5);
	unsigned long y = (unsigned long)(centerPosition.get_y() + 0.5);
	
	long particleIndex = mapped_image->get(x, y);
	
	if ((x < minDistanceFromEdge) || (x > x_size - minDistanceFromEdge - 1))
		return;
	if ((y < minDistanceFromEdge) || (y > y_size - minDistanceFromEdge - 1))
		return;
	
	for (long j = y - 1; j <= y + 1; ++j) {
		for (long i = x - 1; i <= x + 1; ++i) {
			
			if ((i == 0) && (j == 0))
				continue;
			
			if (threshold_image->get(i, j) < 128)
				continue;
			
			if (mapped_image->get(i, j) != -1)
				continue;
			
			// add the current position
			mapped_image->set(i, j, particleIndex);
			currentPosition.set_x((double)i);
			currentPosition.set_y((double)j);
			currentPosition.set_intensity(image->get(i, j));
			positionsInCurrentParticle.push_back(currentPosition);
		}
	}
}


boost::shared_ptr<encap_gsl_matrix_uchar> ThresholdImage_Direct::do_thresholding() {
	unsigned long x_size, y_size;
	double current_value;
	double threshold = parameter;
	
	x_size = CCD_Frame->get_x_size();
	y_size = CCD_Frame->get_y_size();
	
	thresholded_image = boost::shared_ptr<encap_gsl_matrix_uchar>(new encap_gsl_matrix_uchar(x_size, y_size));
	
	for (unsigned long i = 0; i < x_size; i++) {
		for (unsigned long j = 0; j < y_size; j++) {
			current_value = CCD_Frame->get(i, j);
			if (current_value >= threshold) {
				thresholded_image->set(i, j, 255);
			} else {
				thresholded_image->set(i, j, 0);
			}
		}
	}
	
	return thresholded_image;
}

boost::shared_ptr<encap_gsl_matrix_uchar> ThresholdImage_Direct::do_thresholding(boost::shared_ptr<encap_gsl_matrix> image) {
	
	boost::shared_ptr<encap_gsl_matrix_uchar> result_image;
	
	CCD_Frame = image;
	
	result_image = do_thresholding();
	
	return result_image;

}
	

boost::shared_ptr<encap_gsl_matrix_uchar> ThresholdImage_Igor_Iterative::do_thresholding() {
	// we use a built-in Igor operation to handle this case
	// the approach is quite clumsy: we create a temporary wave in Igor that has the same values as the frame we are interested in
	// the igor routine then runs on that frame
	// and next we copy the data back out into a new image
	
	// this is wasteful, but in this way the interface is the same as the other classes
	unsigned long x_size, y_size;
	long dimensionSizes[MAX_DIMENSIONS+1];
	long indices[MAX_DIMENSIONS];
	double value[2];
	int result;
	unsigned char threshold_result;
	
	waveHndl tmp_storage_wave;
	
	x_size = CCD_Frame->get_x_size();
	y_size = CCD_Frame->get_y_size();
	
	dimensionSizes[0] = (long)x_size;
	dimensionSizes[1] = (long)y_size;
	dimensionSizes[2] = 0;
	
	// try to make the wave for the temporary storage
	result = MDMakeWave(&tmp_storage_wave, "tmp_thresh_storage", NULL, dimensionSizes, NT_FP64, 1);
	if (result != 0) {
		throw result;
	}
	
	thresholded_image = boost::shared_ptr<encap_gsl_matrix_uchar>(new encap_gsl_matrix_uchar(x_size, y_size));
	
	// now copy the data
	for (unsigned long i = 0; i < x_size; i++) {
		for (unsigned long j = 0; j < y_size; j++) {
			indices[0] = i;
			indices[1] = j;
			
			value[0] = CCD_Frame->get(i, j);
			
			result = MDSetNumericWavePointValue(tmp_storage_wave, indices, value) ;
			if (result != 0) {
				throw result;
			}
		}
	}
	
	// now let Igor handle the calculation of the threshold
	result = XOPSilentCommand("ImageThreshold /Q /O /M=1 tmp_thresh_storage");
	if (result != 0) {
		throw result;
	}
	
	// now copy the data back to the output matrix
	for (unsigned long i = 0; i < x_size; i++) {
		for (unsigned long j = 0; j < y_size; j++) {
			indices[0] = i;
			indices[1] = j;
			
			result = MDGetNumericWavePointValue(tmp_storage_wave, indices, value) ;
			if (result != 0) {
				throw result;
			}
			
			threshold_result = (unsigned char)(value[0] + 0.5);
			
			thresholded_image->set(i, j, threshold_result);
		}
	}
	
	// now delete the temporary storage wave
	result = KillWave(tmp_storage_wave);
	if (result != 0) {
		throw result;
	}
	
	return thresholded_image;
}

boost::shared_ptr<encap_gsl_matrix_uchar> ThresholdImage_Igor_Iterative::do_thresholding(boost::shared_ptr<encap_gsl_matrix> image) {
	boost::shared_ptr<encap_gsl_matrix_uchar> result_image;
	
	CCD_Frame = image;
	
	result_image = do_thresholding();
	
	return result_image;
}


boost::shared_ptr<encap_gsl_matrix_uchar> ThresholdImage_Igor_Bimodal::do_thresholding() {
	// we use a built-in Igor operation to handle this case
	// the approach is quite clumsy: we create a temporary wave in Igor that has the same values as the frame we are interested in
	// the igor routine then runs on that frame
	// and next we copy the data back out into a new image
	
	// this is wasteful, but the interface is the same as the other classes
	unsigned long x_size, y_size;
	long dimensionSizes[MAX_DIMENSIONS+1];
	long indices[MAX_DIMENSIONS];
	double value[2];
	int result;
	unsigned char threshold_result;
	
	waveHndl tmp_storage_wave;
	
	x_size = CCD_Frame->get_x_size();
	y_size = CCD_Frame->get_y_size();
	
	dimensionSizes[0] = (long)x_size;
	dimensionSizes[1] = (long)y_size;
	dimensionSizes[2] = 0;
	
	// try to make the wave for the temporary storage
	result = MDMakeWave(&tmp_storage_wave, "tmp_thresh_storage", NULL, dimensionSizes, NT_FP64, 1);
	if (result != 0) {
		throw result;
	}
	
	thresholded_image = boost::shared_ptr<encap_gsl_matrix_uchar>(new encap_gsl_matrix_uchar(x_size, y_size));
	
	// now copy the data
	for (unsigned long i = 0; i < x_size; i++) {
		for (unsigned long j = 0; j < y_size; j++) {
			indices[0] = i;
			indices[1] = j;
			
			value[0] = CCD_Frame->get(i, j);
			
			result = MDSetNumericWavePointValue(tmp_storage_wave, indices, value) ;
			if (result != 0) {
				throw result;
			}
		}
	}
	
	// now let Igor handle the calculation of the threshold
	result = XOPSilentCommand("ImageThreshold /Q /O /M=2 tmp_thresh_storage");
	if (result != 0) {
		throw result;
	}
	
	// now copy the data back to the output matrix
	for (unsigned long i = 0; i < x_size; ++i) {
		for (unsigned long j = 0; j < y_size; ++j) {
			indices[0] = i;
			indices[1] = j;
			
			result = MDGetNumericWavePointValue(tmp_storage_wave, indices, value) ;
			if (result != 0) {
				throw result;
			}
			
			threshold_result = (unsigned char)(value[0] + 0.5);
			
			thresholded_image->set(i, j, threshold_result);
		}
	}
	
	// now delete the temporary storage wave
	result = KillWave(tmp_storage_wave);
	if (result != 0) {
		throw result;
	}
	
	return thresholded_image;
}

boost::shared_ptr<encap_gsl_matrix_uchar> ThresholdImage_Igor_Bimodal::do_thresholding(boost::shared_ptr<encap_gsl_matrix> image) {
	boost::shared_ptr<encap_gsl_matrix_uchar> result_image;
	
	CCD_Frame = image;
	
	result_image = do_thresholding();
	
	return result_image;
}


boost::shared_ptr<encap_gsl_matrix_uchar> ThresholdImage_Igor_Adaptive::do_thresholding() {
	// we use a built-in Igor operation to handle this case
	// the approach is quite clumsy: we create a temporary wave in Igor that has the same values as the frame we are interested in
	// the igor routine then runs on that frame
	// and next we copy the data back out into a new image
	
	// this is wasteful, but the interface is the same as the other classes
	unsigned long x_size, y_size;
	long dimensionSizes[MAX_DIMENSIONS+1];
	long indices[MAX_DIMENSIONS];
	double value[2];
	int result;
	unsigned char threshold_result;
	boost::shared_ptr<encap_gsl_matrix_uchar> original_thresholded;
	boost::shared_ptr<encap_gsl_matrix_uchar> transposed_tresholded;
	
	waveHndl tmp_storage_wave;
	waveHndl thresholded_wave;
	
	x_size = CCD_Frame->get_x_size();
	y_size = CCD_Frame->get_y_size();
	
	// we make two images for the original and the transposed threshold
	original_thresholded = boost::shared_ptr<encap_gsl_matrix_uchar>(new encap_gsl_matrix_uchar(x_size, y_size));
	transposed_tresholded = boost::shared_ptr<encap_gsl_matrix_uchar>(new encap_gsl_matrix_uchar(x_size, y_size));
	
	dimensionSizes[0] = (long)x_size;
	dimensionSizes[1] = (long)y_size;
	dimensionSizes[2] = 0;
	
	// try to make the wave for the temporary storage
	result = MDMakeWave(&tmp_storage_wave, "tmp_thresh_storage", NULL, dimensionSizes, NT_FP64, 1);
	if (result != 0) {
		throw result;
	}
	
	// now copy the data
	for (unsigned long i = 0; i < x_size; i++) {
		for (unsigned long j = 0; j < y_size; j++) {
			indices[0] = i;
			indices[1] = j;
			
			value[0] = CCD_Frame->get(i, j);
			
			result = MDSetNumericWavePointValue(tmp_storage_wave, indices, value) ;
			if (result != 0) {
				throw result;
			}
		}
	}
	
	// now let Igor handle the calculation of the threshold
	result = XOPSilentCommand("ImageThreshold /Q /M=3 tmp_thresh_storage");
	if (result != 0) {
		throw result;
	}
	
	
	thresholded_wave = FetchWave("M_ImageThresh");
	if (thresholded_wave == NULL) {
		throw NOWAV;
	}
	
	// now copy the thresholded image back to the first output wave
	for (unsigned long i = 0; i < x_size; i++) {
		for (unsigned long j = 0; j < y_size; j++) {
			indices[0] = i;
			indices[1] = j;
			
			result = MDGetNumericWavePointValue(thresholded_wave, indices, value) ;
			if (result != 0) {
				throw result;
			}
			
			threshold_result = (unsigned char)(value[0] + 0.5);
			
			original_thresholded->set(i, j, threshold_result);
		}
	}
	
	// now transpose the image
	result = XOPSilentCommand("MatrixTranspose tmp_thresh_storage");
	if (result != 0) {
		throw result;
	}
	
	// now calculate the threshold again
	result = XOPSilentCommand("ImageThreshold /Q /M=3 tmp_thresh_storage");
	if (result != 0) {
		throw result;
	}
	
	// we transpose the thresholded image back to the original orientation
	result = XOPSilentCommand("MatrixTranspose M_ImageThresh");
	if (result != 0) {
		throw result;
	}
	
	thresholded_wave = FetchWave("M_ImageThresh");
	if (thresholded_wave == NULL) {
		throw NOWAV;
	}
	
	// now copy the thresholded image back to the second output wave
	for (unsigned long i = 0; i < x_size; i++) {
		for (unsigned long j = 0; j < y_size; j++) {
			indices[0] = i;
			indices[1] = j;
			
			result = MDGetNumericWavePointValue(thresholded_wave, indices, value) ;
			if (result != 0) {
				throw result;
			}
			
			threshold_result = (unsigned char)(value[0] + 0.5);
			
			transposed_tresholded->set(i, j, threshold_result);
		}
	}
	
	// now construct the combined threshold image by AND'ing the original and transpose together
	for (unsigned long i = 0; i < x_size; i++) {
		for (unsigned long j = 0; j < y_size; j++) {
			if (original_thresholded->get(i, j) < 128) {	// below the threshold, we skip it
				continue;
			} else {
				// is the transposed matrix also above the threshold?
				if (transposed_tresholded->get(i, j) < 128) {	// below the threshold. We should not include this point
					original_thresholded->set(i, j, 0);
				} else {
					continue;
				}
			}
		}
	}
	
	
	
	// now delete the temporary storage wave
	result = KillWave(tmp_storage_wave);
	if (result != 0) {
		throw result;
	}
	
	return original_thresholded;
}

boost::shared_ptr<encap_gsl_matrix_uchar> ThresholdImage_Igor_Adaptive::do_thresholding(boost::shared_ptr<encap_gsl_matrix> image) {
	boost::shared_ptr<encap_gsl_matrix_uchar> result_image;
	
	CCD_Frame = image;
	
	result_image = do_thresholding();
	
	return result_image;
}


boost::shared_ptr<encap_gsl_matrix_uchar> ThresholdImage_Igor_Fuzzy1::do_thresholding() {
	// we use a built-in Igor operation to handle this case
	// the approach is quite clumsy: we create a temporary wave in Igor that has the same values as the frame we are interested in
	// the igor routine then runs on that frame
	// and next we copy the data back out into a new image
	
	// this is wasteful, but the interface is the same as the other classes
	unsigned long x_size, y_size;
	long dimensionSizes[MAX_DIMENSIONS+1];
	long indices[MAX_DIMENSIONS];
	double value[2];
	int result;
	unsigned char threshold_result;
	
	waveHndl tmp_storage_wave;
	
	x_size = CCD_Frame->get_x_size();
	y_size = CCD_Frame->get_y_size();
	
	dimensionSizes[0] = (long)x_size;
	dimensionSizes[1] = (long)y_size;
	dimensionSizes[2] = 0;
	
	// try to make the wave for the temporary storage
	result = MDMakeWave(&tmp_storage_wave, "tmp_thresh_storage", NULL, dimensionSizes, NT_FP64, 1);
	if (result != 0) {
		throw result;
	}
	
	thresholded_image = boost::shared_ptr<encap_gsl_matrix_uchar>(new encap_gsl_matrix_uchar(x_size, y_size));
	
	// now copy the data
	for (unsigned long i = 0; i < x_size; i++) {
		for (unsigned long j = 0; j < y_size; j++) {
			indices[0] = i;
			indices[1] = j;
			
			value[0] = CCD_Frame->get(i, j);
			
			result = MDSetNumericWavePointValue(tmp_storage_wave, indices, value) ;
			if (result != 0) {
				throw result;
			}
		}
	}
	
	// now let Igor handle the calculation of the threshold
	result = XOPSilentCommand("ImageThreshold /Q /O /M=4 tmp_thresh_storage");
	if (result != 0) {
		throw result;
	}
	
	// now copy the data back to the output matrix
	for (unsigned long i = 0; i < x_size; i++) {
		for (unsigned long j = 0; j < y_size; j++) {
			indices[0] = i;
			indices[1] = j;
			
			result = MDGetNumericWavePointValue(tmp_storage_wave, indices, value) ;
			if (result != 0) {
				throw result;
			}
			
			threshold_result = (unsigned char)(value[0] + 0.5);
			
			thresholded_image->set(i, j, threshold_result);
		}
	}
	
	// now delete the temporary storage wave
	result = KillWave(tmp_storage_wave);
	if (result != 0) {
		throw result;
	}
	
	return thresholded_image;
}

boost::shared_ptr<encap_gsl_matrix_uchar> ThresholdImage_Igor_Fuzzy1::do_thresholding(boost::shared_ptr<encap_gsl_matrix> image) {
	boost::shared_ptr<encap_gsl_matrix_uchar> result_image;
	
	CCD_Frame = image;
	
	result_image = do_thresholding();
	
	return result_image;
}


boost::shared_ptr<encap_gsl_matrix_uchar> ThresholdImage_Igor_Fuzzy2::do_thresholding() {
	// we use a built-in Igor operation to handle this case
	// the approach is quite clumsy: we create a temporary wave in Igor that has the same values as the frame we are interested in
	// the igor routine then runs on that frame
	// and next we copy the data back out into a new image
	
	// this is wasteful, but the interface is the same as the other classes
	unsigned long x_size, y_size;
	long dimensionSizes[MAX_DIMENSIONS+1];
	long indices[MAX_DIMENSIONS];
	double value[2];
	int result;
	unsigned char threshold_result;
	
	waveHndl tmp_storage_wave;
	
	x_size = CCD_Frame->get_x_size();
	y_size = CCD_Frame->get_y_size();
	
	dimensionSizes[0] = (long)x_size;
	dimensionSizes[1] = (long)y_size;
	dimensionSizes[2] = 0;
	
	// try to make the wave for the temporary storage
	result = MDMakeWave(&tmp_storage_wave, "tmp_thresh_storage", NULL, dimensionSizes, NT_FP64, 1);
	if (result != 0) {
		throw result;
	}
	
	thresholded_image = boost::shared_ptr<encap_gsl_matrix_uchar>(new encap_gsl_matrix_uchar(x_size, y_size));
	
	// now copy the data
	for (unsigned long i = 0; i < x_size; i++) {
		for (unsigned long j = 0; j < y_size; j++) {
			indices[0] = i;
			indices[1] = j;
			
			value[0] = CCD_Frame->get(i, j);
			
			result = MDSetNumericWavePointValue(tmp_storage_wave, indices, value) ;
			if (result != 0) {
				throw result;
			}
		}
	}
	
	// now let Igor handle the calculation of the threshold
	result = XOPSilentCommand("ImageThreshold /Q /O /M=5 tmp_thresh_storage");
	if (result != 0) {
		throw result;
	}
	
	// now copy the data back to the output matrix
	for (unsigned long i = 0; i < x_size; i++) {
		for (unsigned long j = 0; j < y_size; j++) {
			indices[0] = i;
			indices[1] = j;
			
			result = MDGetNumericWavePointValue(tmp_storage_wave, indices, value) ;
			if (result != 0) {
				throw result;
			}
			
			threshold_result = (unsigned char)(value[0] + 0.5);
			
			thresholded_image->set(i, j, threshold_result);
		}
	}
	
	// now delete the temporary storage wave
	result = KillWave(tmp_storage_wave);
	if (result != 0) {
		throw result;
	}
	
	return thresholded_image;
}

boost::shared_ptr<encap_gsl_matrix_uchar> ThresholdImage_Igor_Fuzzy2::do_thresholding(boost::shared_ptr<encap_gsl_matrix> image) {
	boost::shared_ptr<encap_gsl_matrix_uchar> result_image;
	
	CCD_Frame = image;
	
	result_image = do_thresholding();
	
	return result_image;
}

boost::shared_ptr<encap_gsl_matrix_uchar> ThresholdImage_Isodata::do_thresholding() {
	
	gsl_histogram *hist;
	boost::shared_ptr<encap_gsl_matrix_uchar> threshold_image;
	
	unsigned long x_size = CCD_Frame->get_x_size();
	unsigned long y_size = CCD_Frame->get_y_size();
	
	unsigned long number_of_bins = 256;
	unsigned long current_threshold_bin = 127;
	int n_iterations = 0;
	int max_iters = 50;
	double lower_mean, upper_mean;
	double sum, denominator_sum;
	double current_threshold = -1, previous_threshold;
	double lower_bin_limit, upper_bin_limit;
	unsigned long bin_threshold;
	double intensity_threshold;
	int converged = 0;
	
	// since this is a histogram-based approach we start by constructing the histogram
	hist = make_histogram_from_matrix(CCD_Frame, number_of_bins);
	
	// because this approach is based on thresholding it makes sense to only express the threshold in terms of bins, stored in "current_threshold_bin".
	// a value of 0 for "current_threshold_bin" means that all bins at index 0 and higher are considered to be 'signal', not 'background'.
	
	while ((converged == 0) && (n_iterations < max_iters)) {
		
		previous_threshold = current_threshold;
		
		// calculate the lower mean
		sum = 0;
		denominator_sum = 0;
		for (unsigned long i = 0; i < current_threshold_bin; i++) {
			sum += (double)i * gsl_histogram_get(hist, i);
			denominator_sum += gsl_histogram_get(hist, i);
		}
	
		lower_mean = sum / denominator_sum;
		
		// calculate the upper mean
		sum = 0;
		denominator_sum = 0;
		for (unsigned long i = current_threshold_bin; i < number_of_bins; i++) {
			sum += (double)i * gsl_histogram_get(hist, i);
			denominator_sum += gsl_histogram_get(hist, i);
		}
		
		upper_mean = sum / denominator_sum;
		
		current_threshold = (lower_mean + upper_mean) / 2.0;
		current_threshold_bin = floor(current_threshold);
		
		if (floor(current_threshold + 0.5) == floor(previous_threshold + 0.5)) {
			bin_threshold = floor(current_threshold + 0.5);
			converged = 1;
		}
		
	}
	
	if (converged == 0) {	// the iterations did not converge, there is no clear threshold
							// to indicate this we set everything to 'off' (0)
		threshold_image = boost::shared_ptr<encap_gsl_matrix_uchar>(new encap_gsl_matrix_uchar(x_size, y_size));
		threshold_image->set_all(0);
		return threshold_image;
	}
	
	// now translate the threshold value to an intensity instead of being in bins
	gsl_histogram_get_range(hist, bin_threshold, &lower_bin_limit, &upper_bin_limit);
	intensity_threshold = lower_bin_limit;
	
	// get another threshold class to do the work for us
	ThresholdImage_Direct thresholder(CCD_Frame, intensity_threshold);
	
	try {
		threshold_image = thresholder.do_thresholding();
	}
	catch (OUT_OF_MEMORY) {
		gsl_histogram_free(hist);
		string error;
		error = "unable to allocate threshold_image in ThresholdImage_Isodata::do_thresholding()\r";
		throw OUT_OF_MEMORY(error);
	}
	
	gsl_histogram_free(hist);
	
	return threshold_image;
}

boost::shared_ptr<encap_gsl_matrix_uchar> ThresholdImage_Isodata::do_thresholding(boost::shared_ptr<encap_gsl_matrix> image) {
	boost::shared_ptr<encap_gsl_matrix_uchar> result_image;
	
	CCD_Frame = image;
	
	result_image = do_thresholding();
	
	return result_image;
}

boost::shared_ptr<encap_gsl_matrix_uchar> ThresholdImage_Triangle::do_thresholding() {
	gsl_histogram *hist;
	unsigned long number_of_bins = 256;
	unsigned long maximum_bin;
	double max_val, max_bin_double;
	double end_val, end_bin_double;
	double slope, intercept, perpendicular_slope, perpendicular_intercept;
	double current_bin_value, double_i;
	double intercept_x, intercept_y;
	double distance;
	double max_distance = -1;
	unsigned long max_index;
	double lower_bin_limit, upper_bin_limit, intensity_threshold;
	
	boost::shared_ptr<encap_gsl_matrix_uchar> threshold_image;
	unsigned long x_size = CCD_Frame->get_x_size();
	unsigned long y_size = CCD_Frame->get_y_size();
	
	// since this is a histogram-based approach we start by constructing the histogram
	
	try {
		hist = make_histogram_from_matrix(CCD_Frame, number_of_bins);
	}
	catch (OUT_OF_MEMORY) {
		string error;
		error = "unable to allocate buffer in ThresholdImage_Triangle::do_thresholding()\r";
		throw OUT_OF_MEMORY(error);
	}
	
	maximum_bin = gsl_histogram_max_bin(hist);
	max_bin_double = (double)maximum_bin;
	max_val = gsl_histogram_max_val(hist);
	end_val = gsl_histogram_get(hist, number_of_bins - 1);	// the bin that contains the largest intensity is the last bin in the histogram
	end_bin_double = (double)(number_of_bins - 1);
	
	// catch an unlikely case where the maximum corresponds to the last bin
	if (maximum_bin == (number_of_bins - 1)) {
		gsl_histogram_free(hist);
		threshold_image = boost::shared_ptr<encap_gsl_matrix_uchar>(new encap_gsl_matrix_uchar(x_size, y_size));
		threshold_image->set_all(0);
		return threshold_image;
	}
	
	// calculate the line that connects the maximum and highest-intensity value
	slope = (end_val - max_val) / (end_bin_double - max_bin_double);
	intercept = max_val / (slope * max_bin_double);
	
	// catch an unlikely case where the connecting line is flat (the histogram is apparently uniform)
	if (slope == 0) {
		gsl_histogram_free(hist);
		threshold_image = boost::shared_ptr<encap_gsl_matrix_uchar>(new encap_gsl_matrix_uchar(x_size, y_size));
		threshold_image->set_all(0);
		return threshold_image;
	}
	
	
	// calculate the slope of a line perpendicular to the connecting line
	perpendicular_slope = - 1.0 / slope;
	
	for (unsigned long i = maximum_bin + 1; i < number_of_bins; i++) {	// determine the minimum distance in the triangle
		
		// what is the intercept for the perpendicular line if it has to go through the bin that we're currently looking at?
		current_bin_value = gsl_histogram_get(hist, i);
		double_i = (double)i;
		
		perpendicular_intercept = current_bin_value / (perpendicular_slope * double_i);
		
		// where does the perpendicular line intercept the connecting line?
		// x = (b1 - b2) / (a1 - a2) and y = a1 * x + b1
		intercept_x = (intercept - perpendicular_intercept) / (slope - perpendicular_slope);
		intercept_y = slope * intercept_x + intercept;
		
		// what is the distance to the connecting line?
		distance = sqrt((intercept_x - double_i) * (intercept_x - double_i) + (intercept_y - current_bin_value) * (intercept_y - current_bin_value));
		
		if (distance > max_distance) {
			max_index = i;
			max_distance = distance;
		}
	}
	
	// translate the maximal index to a threshold value
	gsl_histogram_get_range(hist, max_index, &lower_bin_limit, &upper_bin_limit);
	intensity_threshold = lower_bin_limit;
	
	// get another threshold class to do the work for us
	ThresholdImage_Direct thresholder(CCD_Frame, intensity_threshold);
	
	try {
		threshold_image = thresholder.do_thresholding();
	}
	catch (OUT_OF_MEMORY) {
		string error;
		error = "unable to do direct thresholding in ThresholdImage_Triangle::do_thresholding()\r";
		gsl_histogram_free(hist);
		throw OUT_OF_MEMORY(error);
	}
	
	gsl_histogram_free(hist);
	
	return threshold_image;
}

boost::shared_ptr<encap_gsl_matrix_uchar> ThresholdImage_Triangle::do_thresholding(boost::shared_ptr<encap_gsl_matrix> image) {
	boost::shared_ptr<encap_gsl_matrix_uchar> result_image;
	
	CCD_Frame = image;
	
	result_image = do_thresholding();
	
	return result_image;
}


boost::shared_ptr<encap_gsl_matrix_uchar> ThresholdImage_MTT::do_thresholding() {
	
	// the code is based on a series of matlab files sent by Didier Marguet, corresponding author of the original paper
	boost::shared_ptr<encap_gsl_matrix_uchar> threshold_image;
	unsigned long x_size = CCD_Frame->get_x_size();
	unsigned long y_size = CCD_Frame->get_y_size();
	
	boost::shared_ptr<encap_gsl_matrix> averages(new encap_gsl_matrix (x_size, y_size));	// contains an estimation of the average value at every position, calculation over the number of pixels in the window
	boost::shared_ptr<encap_gsl_matrix> CCD_Frame_squared(new encap_gsl_matrix (x_size, y_size));
	boost::shared_ptr<encap_gsl_matrix> summed_squares(new encap_gsl_matrix (x_size, y_size));
	boost::shared_ptr<encap_gsl_matrix> null_hypothesis(new encap_gsl_matrix (x_size, y_size));
	boost::shared_ptr<encap_gsl_matrix> Gaussian_window(new encap_gsl_matrix (x_size, y_size));
	boost::shared_ptr<encap_gsl_matrix> image_Gaussian_convolved(new encap_gsl_matrix (x_size, y_size));	// this is 'alpha' in the original matlab code
	boost::shared_ptr<encap_gsl_matrix> hypothesis_test(new encap_gsl_matrix (x_size, y_size));	// this is 'test' in the original matlab code
	
	parameter = 13;	// TODO: test different values for the window size, and allow a way for this to be set from the command line
	
	double average = 0;
	unsigned long window_size = parameter;	// this is the size of the window over which we calculate the hypotheses
	unsigned long half_window_size = window_size / 2;	// integer division takes care of the floor() aspect
	unsigned long window_pixels = window_size * window_size;
	double double_window_pixels = (double)(window_pixels);
	double current_value;
	double distance_x, distance_y;
	double sum;
	double sum_squared_Gaussian;
	
	
	threshold_image = boost::shared_ptr<encap_gsl_matrix_uchar>(new encap_gsl_matrix_uchar(x_size, y_size));
	threshold_image->set_all(0);
	
	averages->set_all(0);
	summed_squares->set_all(0);
	hypothesis_test->set_all(0);
	
	// calculate the square of the pixel values
	// we'll use this later
	for (unsigned long i = 0; i < x_size; i++) {
		for (unsigned long j = 0; j < y_size; j++) {
			current_value = CCD_Frame->get(i, j);
			current_value = current_value * current_value;
			CCD_Frame_squared->set(i, j, current_value);
		}
	}
	
	// NULL HYPOTHESIS: there is no emitter at a certain position
	
	// start by estimating the mean at every position
	// in the original code this done by a convolution of a unity matrix with the window size-> This is done using an FFT in the matlab code, but we'll do it directly
	
	for (unsigned long k = half_window_size; k < x_size - half_window_size; k++) {
		for (unsigned long l = half_window_size; l < y_size - half_window_size; l++) {
			// these two loops run over the entire image-> We exclude the edges where the window doesn't fit
			
			// now we loop over the size of the window, to determine the average at every pixel
			average = 0;
			
			for (unsigned long j = l - half_window_size; j <= l + half_window_size; j++) {
				for (unsigned long i = k - half_window_size; i <= k + half_window_size; i++) {
					average += CCD_Frame->get(i, j);
				}
			}
			average /= double_window_pixels;
			
			averages->set(k, l, average);
		}
	}
	
	// now we do the same, but for the square of the pixel values, and we don't divide (no averaging)
	for (unsigned long k = half_window_size; k < x_size - half_window_size; k++) {
		for (unsigned long l = half_window_size; l < y_size - half_window_size; l++) {
			// these two loops run over the entire image-> We exclude the edges where the window doesn't fit
			current_value = 0;
			for (unsigned long i = k - half_window_size; i <= k + half_window_size; i++) {
				for (unsigned long j = l - half_window_size; j <= l + half_window_size; j++) {
					current_value += CCD_Frame_squared->get(i, j);
				}
			}
			summed_squares->set(k, l, current_value);	// summed_squares is now equal to "Sim2" in the orignal matlab code
		}
	}
	
	// now calculate the null hypothesis image-> This is T_sig0_2 in the original matlab source
	for (unsigned long l = half_window_size; l < y_size - half_window_size; l++) {
		for (unsigned long k = half_window_size; k < x_size - half_window_size; k++) {
			current_value = summed_squares->get(k, l) - double_window_pixels * averages->get(k, l) * averages->get(k, l);
			null_hypothesis->set(k, l, current_value);
		}
	}
	
	// calculate the hypothesis H1 that there is an emitter
	
	// store the values of a Gaussian in a matrix with the size of a window
	// this is so we can cache the values, and do not to calculate it every time
	sum = 0;
	for (unsigned long j = 0; j < window_size; j++) {
		for (unsigned long i = 0; i < window_size; i++) {
			// the Gaussian is assumed to be in the center of the window
			distance_x = (double)half_window_size - (double)i;
			distance_y = (double)half_window_size - (double)j;
			current_value = 1.0 / (1.77245385 * gaussianWidth) * exp(- 1.0 / (2.0 * gaussianWidth * gaussianWidth) * (distance_x * distance_x + distance_y * distance_y));
			
			Gaussian_window->set(i, j, current_value);
			
			sum += current_value;	// we will use this below
		}
	}
	
	// now we re-normalize this Gaussian matrix
	// at this point Gaussian_window becomes equal to 'gc' in the original matlab code
	sum /= double_window_pixels;
	sum_squared_Gaussian = 0;
	for (unsigned long j = 0; j < window_size; j++) {
		for (unsigned long i = 0; i < window_size; i++) {
			current_value = Gaussian_window->get(i, j);
			current_value = current_value - sum;
			Gaussian_window->set(i, j, current_value);
			sum_squared_Gaussian += current_value * current_value;	// this is 'Sgc2' in the original code
		}
	}
	
	// now we need to again convolve this Gaussian_window ('gc') with the original image-> As before, we'll do it directly
	// these two loops run over the entire image-> We exclude the edges where the window doesn't fit
	for (unsigned long l = half_window_size; l < y_size - half_window_size; l++) {
		for (unsigned long k = half_window_size; k < x_size - half_window_size; k++) {
			
			sum = 0;
			for (unsigned long j = l - half_window_size; j <= l + half_window_size; j++) {
				for (unsigned long i = k - half_window_size; i <= k + half_window_size; i++) {
					current_value = Gaussian_window->get(i - k + half_window_size, j - l + half_window_size);
					current_value *= CCD_Frame->get(i, j);
					sum += current_value;
				}
			}
			sum /= sum_squared_Gaussian;
			image_Gaussian_convolved->set(k, l, sum);
		}
	}
	
	// "image_Gaussian_convolved" is now equal to 'alpha' in the original matlab code
	
	// calculate the image that will determine whether to accept or reject the null hypothesis
	for (unsigned long l = half_window_size; l < y_size - half_window_size; l++) {
		for (unsigned long k = half_window_size; k < x_size - half_window_size; k++) {
			current_value = 1 - (sum_squared_Gaussian * image_Gaussian_convolved->get(k, l) * image_Gaussian_convolved->get(k, l)) / null_hypothesis->get(k , l);
			current_value = (current_value > 0) * current_value + (current_value <= 0);	// the equivalent of test = (test > 0) ->* test + (test <= 0) in the original code
			current_value = - double_window_pixels * log(current_value);
			hypothesis_test->set(k, l, current_value);
		}
	}
	
	// at this point 'hypothesis_test' is equal to 'carte_MV' in the original image
	// check where we have to reject the hypothesis test
	for (unsigned long l = half_window_size; l < y_size - half_window_size; l++) {
		for (unsigned long k = half_window_size; k < x_size - half_window_size; k++) {
			if (hypothesis_test->get(k, l) > PFA) {
				threshold_image->set(k, l, 255);
			}
		}
	}
	
	return threshold_image;
}


boost::shared_ptr<encap_gsl_matrix_uchar> ThresholdImage_MTT::do_thresholding(boost::shared_ptr<encap_gsl_matrix> image) {
	boost::shared_ptr<encap_gsl_matrix_uchar> result_image;
	
	CCD_Frame = image;
	
	result_image = do_thresholding();
	
	return result_image;
}


boost::shared_ptr<encap_gsl_matrix_uchar> ThresholdImage_MTT_FFT::do_thresholding() {
	
	// the code is based on a series of matlab files sent by Didier Marguet, corresponding author of the original paper
	boost::shared_ptr<encap_gsl_matrix_uchar> threshold_image;
	unsigned long x_size = CCD_Frame->get_x_size();
	unsigned long y_size = CCD_Frame->get_y_size();
	
	boost::shared_ptr<encap_gsl_matrix> averages;
	boost::shared_ptr<encap_gsl_matrix> CCD_Frame_squared;
	boost::shared_ptr<encap_gsl_matrix> summed_squares;
	boost::shared_ptr<encap_gsl_matrix> null_hypothesis;
	boost::shared_ptr<encap_gsl_matrix> image_Gaussian_convolved;
	boost::shared_ptr<encap_gsl_matrix> hypothesis_test;
	boost::shared_ptr<encap_gsl_matrix> Gaussian_window;
	
	unsigned long window_size = 13;
	unsigned long half_window_size = window_size / 2;	// integer division takes care of the floor() aspect
	unsigned long center_x = x_size / 2;
	unsigned long center_y = y_size / 2;
	unsigned long window_pixels = window_size * window_size;
	double double_window_pixels = (double)(window_pixels);
	double current_value;
	double sum_squared_Gaussian;
	double distance_x, distance_y;
	double sum;
	
	threshold_image = boost::shared_ptr<encap_gsl_matrix_uchar>(new encap_gsl_matrix_uchar(x_size, y_size));
	threshold_image->set_all(0);
	
	averages = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(x_size, y_size));
	averages->set_all(0);
	
	CCD_Frame_squared = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(x_size, y_size));
	
	summed_squares = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(x_size, y_size));
	summed_squares->set_all(0);
	
	null_hypothesis = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(x_size, y_size));
	
	Gaussian_window = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(x_size, y_size));
	
	image_Gaussian_convolved = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(x_size, y_size));	// this is 'alpha' in the original matlab code
	
	hypothesis_test = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(x_size, y_size));	// this is 'test' in the original matlab code
	hypothesis_test->set_all(0);
	
	// calculate the square of the pixel values
	// we'll use this later
	for (unsigned long j = 0; j < y_size; j++) {
		for (unsigned long i = 0; i < x_size; i++) {
			current_value = CCD_Frame->get(i, j);
			current_value = current_value * current_value;
			CCD_Frame_squared->set(i, j, current_value);
		}
	}
	
	// NULL HYPOTHESIS: there is no emitter at a certain position
	
	// start by estimating the mean at every position
	// in the original code this done by a convolution of a unity matrix with the window size. We now do this using an FFT-based approach
	
	
	// convolve the image with a "box function", that will get us the average
	// if we have a lot of images then we only need to make this kernel once
	if ((average_kernel.get() == NULL) || (kernel_x_size != x_size) || (kernel_y_size != y_size)) {
		average_kernel = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(x_size, y_size));
		
		average_kernel->set_all(0);
		
		for (unsigned long j = center_y - half_window_size; j <= center_y + half_window_size; j++) {
			for (unsigned long i = center_x - half_window_size; i <= center_x + half_window_size; i++) {
				average_kernel->set(i, j, 1);
			}
		}
	}
	
	averages = convolve_matrices_using_fft(CCD_Frame, average_kernel);
	
	// normalize the result, so that we get averages
	for (unsigned long i = 0; i < x_size; i++) {
		for (unsigned long j = 0; j < y_size; j++) {
			current_value = averages->get(i, j);
			current_value /= double_window_pixels;
			averages->set(i, j, current_value);
		}
	}
	
	// do the same for the squared CCD_Frame
	summed_squares = convolve_matrices_using_fft(CCD_Frame_squared, average_kernel);
	
	// now calculate the null hypothesis image. This is T_sig0_2 in the original matlab source
	for (unsigned long l = half_window_size; l < y_size - half_window_size; l++) {
		for (unsigned long k = half_window_size; k < x_size - half_window_size; k++) {
			current_value = summed_squares->get(k, l) - double_window_pixels * averages->get(k, l) * averages->get(k, l);
			null_hypothesis->set(k, l, current_value);
		}
	}
	
	// calculate the hypothesis H1 that there is an emitter
	
	 // create a Gaussian kernel for the convolution
	// we only need to do this once if we are looking at a series of images
	if ((Gaussian_kernel.get() == NULL) || (kernel_x_size != x_size) || (kernel_y_size != y_size)) {		
		Gaussian_kernel = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(x_size, y_size));
		
		Gaussian_window = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(window_size, window_size));
		
		sum = 0;
		for (unsigned long j = 0; j < window_size; j++) {
			for (unsigned long i = 0; i < window_size; i++) {
				// the Gaussian is assumed to be in the center of the window
				distance_x = (double)half_window_size - (double)i;
				distance_y = (double)half_window_size - (double)j;
				current_value = 1.0 / (1.77245385 * gaussianWidth) * exp(- 1.0 / (2.0 * gaussianWidth * gaussianWidth) * (distance_x * distance_x + distance_y * distance_y));
				
				Gaussian_window->set(i, j, current_value);
				
				sum += current_value;	// we will use this below
			}
		}
		
		// now we re-normalize this Gaussian matrix
		// at this point Gaussian_window becomes equal to 'gc' in the original matlab code
		sum /= double_window_pixels;
		sum_squared_Gaussian = 0;
		for (unsigned long j = 0; j < window_size; j++) {
			for (unsigned long i = 0; i < window_size; i++) {
				current_value = Gaussian_window->get(i, j);
				current_value = current_value - sum;
				Gaussian_window->set(i, j, current_value);
				sum_squared_Gaussian += current_value * current_value;	// this is 'Sgc2' in the original code
			}
		}
		
		// now introduce this small kernel into a larger one that is the same size as the image
		Gaussian_kernel->set_all(0);
		
		for (unsigned long j = center_y - half_window_size; j <= center_y + half_window_size; j++) {
			for (unsigned long i = center_x - half_window_size; i <= center_x + half_window_size; i++) {
				Gaussian_kernel->set(i, j, Gaussian_window->get(i - center_x + half_window_size, j - center_y + half_window_size));
			}
		}
	}
	
	// now we need to again convolve this Gaussian_window ('gc') with the original image. 
	// we now do this using the FFT
	
	image_Gaussian_convolved = convolve_matrices_using_fft(CCD_Frame, Gaussian_kernel);
	
	// now normalize this convolved image so that it becomes equal to 'alpha' in the original matlab code
	for (unsigned long i = 0; i < x_size; i++) {
		for (unsigned long j = 0; j < y_size; j++) {
			current_value = image_Gaussian_convolved->get(i, j);
			current_value /= sum_squared_Gaussian;
			image_Gaussian_convolved->set(i, j, current_value);
		}
	}
	
	// calculate the image that will determine whether to accept or reject the null hypothesis
	for (unsigned long l = half_window_size; l < y_size - half_window_size; l++) {
		for (unsigned long k = half_window_size; k < x_size - half_window_size; k++) {
			current_value = 1 - (sum_squared_Gaussian * image_Gaussian_convolved->get(k, l) * image_Gaussian_convolved->get(k, l)) / null_hypothesis->get(k , l);
			current_value = (current_value > 0) * current_value + (current_value <= 0);	// the equivalent of test = (test > 0) .* test + (test <= 0) in the original code
			current_value = - double_window_pixels * log(current_value);
			hypothesis_test->set(k, l, current_value);
		}
	}
	
	// at this point 'hypothesis_test' is equal to 'carte_MV' in the original image
	// check where we have to reject the hypothesis test
	for (unsigned long l = half_window_size; l < y_size - half_window_size; l++) {
		for (unsigned long k = half_window_size; k < x_size - half_window_size; k++) {
			if (hypothesis_test->get(k, l) > PFA) {
				threshold_image->set(k, l, 255);
			}
		}
	}
	
	return threshold_image;
}


boost::shared_ptr<encap_gsl_matrix_uchar> ThresholdImage_MTT_FFT::do_thresholding(boost::shared_ptr<encap_gsl_matrix> image) {
	boost::shared_ptr<encap_gsl_matrix_uchar> result_image;
	
	CCD_Frame = image;
	
	result_image = do_thresholding();
	
	return result_image;
}


boost::shared_ptr<encap_gsl_matrix> ThresholdImage_Preprocessor_MedianFilter::do_preprocessing(boost::shared_ptr<encap_gsl_matrix> image) {
	
	unsigned long kernel_size = kernel_x_size * kernel_y_size;
	unsigned long half_kernel_size_x = kernel_x_size / 2;
	unsigned long half_kernel_size_y = kernel_y_size / 2;
	unsigned long x_size = image->get_x_size();
	unsigned long y_size = image->get_y_size();
	unsigned long offset;
	double value, median;
	
	gsl_vector *median_environment;
	boost::shared_ptr<encap_gsl_matrix> filtered_image;
	
	// allocate a gsl_vector with the correct size
	median_environment = gsl_vector_alloc(kernel_size);
	unsigned long sorted_center = kernel_size / 2;
	
	// make a copy of the image
	// this copy will be median-filtered
	// close to the edges (where the kernel doesn't fit we will not modify the image)
	filtered_image = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(x_size, y_size));
	
	for (unsigned long j = 0; j < y_size; ++j) {
		for (unsigned long i = 0; i < x_size; ++i) {
			filtered_image->set(i, j, image->get(i, j));
		}
	}
	
	
	// the main loop
	// for now we handle this using a direct, slow calculation
	
	for (unsigned long j = half_kernel_size_y; j < y_size - half_kernel_size_y; j++) {
		for (unsigned long i = half_kernel_size_x; i < x_size - half_kernel_size_x; i++) {
			
			offset = 0;
			for (unsigned long l = j - half_kernel_size_y; l <= j + half_kernel_size_y; l++) {
				for (unsigned long k = i - half_kernel_size_x; k <= i + half_kernel_size_x; k++) {
					value = image->get(k, l);
					gsl_vector_set(median_environment, offset, value);
					offset++;
				}
			}
			gsl_sort_vector(median_environment);
			median = gsl_vector_get(median_environment, sorted_center);
			filtered_image->set(i, j, median);
		}
	}
	
	gsl_vector_free(median_environment);
	
	return filtered_image;
}


void ThresholdImage_Preprocessor_GaussianSmoothing::generate_Gaussian_kernel(unsigned long x_size, unsigned long y_size) {
	
	unsigned long window_size = 31;
	unsigned long half_window_size = window_size / 2;
	unsigned long center_x = x_size / 2;
	unsigned long center_y = y_size / 2;
	double current_value, distance_x, distance_y;
	
	boost::shared_ptr<encap_gsl_matrix> Gaussian_window(new encap_gsl_matrix(window_size, window_size));
	
	Gaussian_kernel = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(x_size, y_size));
	
	
	// calculate the values of a Gaussian with the correct width in a smaller window
	for (unsigned long j = 0; j < window_size; j++) {
		for (unsigned long i = 0; i < window_size; i++) {
			// the Gaussian is assumed to be in the center of the window
			distance_x = (double)half_window_size - (double)i;
			distance_y = (double)half_window_size - (double)j;
			current_value = 1.0 / (6.28318531 * width * width) * exp(- 1.0 / (2.0 * width * width) * (distance_x * distance_x + distance_y * distance_y));
			// normalized Gaussian in two dimensions
			
			Gaussian_window->set(i, j, current_value);
		}
	}
	
	// now introduce this small kernel into a larger one that is the same size as the image
	Gaussian_kernel->set_all(0);
	
	for (unsigned long j = center_y - half_window_size; j <= center_y + half_window_size; j++) {
		for (unsigned long i = center_x - half_window_size; i <= center_x + half_window_size; i++) {
			Gaussian_kernel->set(i, j, Gaussian_window->get(i - center_x + half_window_size, j - center_y + half_window_size));
		}
	}
}
	
	


boost::shared_ptr<encap_gsl_matrix> ThresholdImage_Preprocessor_GaussianSmoothing::do_preprocessing(boost::shared_ptr<encap_gsl_matrix> image) {
	
	boost::shared_ptr<encap_gsl_matrix> filtered_image;
	unsigned long x_size = image->get_x_size();
	unsigned long y_size = image->get_y_size();
	
	// do we already have a Gaussian kernel stored, or is this the first run?
	if (Gaussian_kernel.get() == NULL) {	// we don't have a kernel, we need to generate it
		
		generate_Gaussian_kernel(x_size, y_size);
		
	} else {	// we already have a kernel stored, is it the correct size?
				// if not we will calculate a new one
		if ((x_size != Gaussian_kernel->get_x_size()) || (y_size != Gaussian_kernel->get_y_size())) {
			generate_Gaussian_kernel(x_size, y_size);
		}
	}
	
	filtered_image = convolve_matrices_using_fft(image, Gaussian_kernel);
	
	return filtered_image;
}


boost::shared_ptr<encap_gsl_matrix> ThresholdImage_Preprocessor_MeanFilter::do_preprocessing(boost::shared_ptr<encap_gsl_matrix> image) {
	
	unsigned long kernel_size = kernel_x_size * kernel_y_size;
	double double_kernel_pixels = (double)kernel_size;
	unsigned long half_kernel_size_x = kernel_x_size / 2;
	unsigned long half_kernel_size_y = kernel_y_size / 2;
	unsigned long x_size = image->get_x_size();
	unsigned long y_size = image->get_y_size();
	double mean;
	
	boost::shared_ptr<encap_gsl_matrix> filtered_image;
	
	// make a copy of the image
	// this copy will be mean-filtered
	// close to the edges, where the kernel doesn't fit we will not modify the image
	filtered_image = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(x_size, y_size));
	
	for (unsigned long i = 0; i < x_size; ++i) {
		for (unsigned long j = 0; j < y_size; ++j) {
			filtered_image->set(i, j, image->get(i, j));
		}
	}
	
	
	// the main loop
	// for now we handle this using a direct, slow calculation
	
	for (unsigned long j = half_kernel_size_y; j < y_size - half_kernel_size_y; j++) {
		for (unsigned long i = half_kernel_size_x; i < x_size - half_kernel_size_x; i++) {
			
			mean = 0;
			for (unsigned long l = j - half_kernel_size_y; l <= j + half_kernel_size_y; l++) {
				for (unsigned long k = i - half_kernel_size_x; k <= i + half_kernel_size_x; k++) {
					mean += image->get(k, l);
				}
			}
			mean /= double_kernel_pixels;
			filtered_image->set(i, j, mean);
		}
	}
	
	return filtered_image;
}
	
			
boost::shared_ptr<encap_gsl_matrix_uchar> ThresholdImage_Postprocessor_RemoveIsolatedPixels::do_postprocessing(boost::shared_ptr<encap_gsl_matrix_uchar> thresholded_image, boost::shared_ptr<encap_gsl_matrix> image) {
	// we don't care about the edges, they are ignored anyway in the fitting
	unsigned long x_size = thresholded_image->get_x_size();
	unsigned long y_size = thresholded_image->get_y_size();
	unsigned char value;
	double meanIntensity = 0;
	
	boost::shared_ptr<encap_gsl_matrix_uchar> processed_thresholded_image;
	
	processed_thresholded_image = boost::shared_ptr<encap_gsl_matrix_uchar>(new encap_gsl_matrix_uchar(x_size, y_size));
	
	processed_thresholded_image->set_all(0);
	
	// calculate the mean intensity
	for (unsigned long i = 0; i < x_size; i++) {
		for (unsigned long j = 0; j < y_size; j++) {
			meanIntensity += image->get(i, j);
		}
	}
	meanIntensity /= x_size * y_size;
	
	for (unsigned long i = 0; i < x_size; i++) {
		for (unsigned long j = 0; j < y_size; j++) {
			
			value = thresholded_image->get(i, j);
			if (value < 128) {	// this is an 'off' pixel
				continue;
			}
			
			// if we are here then the current pixel is active
			// we only mark the pixel as active in the processed image if it is above the mean
			if (image->get(i, j) > meanIntensity) {
				processed_thresholded_image->set(i, j, 255);
			}
		}
	}
	
	return processed_thresholded_image;
}


boost::shared_ptr<encap_gsl_matrix_uchar> ThresholdImage_Postprocessor_RemovePixelsBelowMean::do_postprocessing(boost::shared_ptr<encap_gsl_matrix_uchar> thresholded_image, boost::shared_ptr<encap_gsl_matrix> image) {
	// we don't care about the edges, they are ignored anyway in the fitting
	unsigned long x_size = thresholded_image->get_x_size();
	unsigned long y_size = thresholded_image->get_y_size();
	unsigned char value;
	bool neighbour_found;
	
	boost::shared_ptr<encap_gsl_matrix_uchar> processed_thresholded_image;
	
	processed_thresholded_image = boost::shared_ptr<encap_gsl_matrix_uchar>(new encap_gsl_matrix_uchar(x_size, y_size));
	
	for (unsigned long i = 0; i < x_size; ++i) {
		for (unsigned long j = 0; j < y_size; ++j) {
			processed_thresholded_image->set(i, j, thresholded_image->get(i, j));
		}
	}
	
	// we will return a copy
	
	for (unsigned long j = 1; j < y_size - 1; j++) {
		for (unsigned long i = 1; i < x_size - 1; i++) {
			
			value = thresholded_image->get(i, j);
			if (value < 128) {	// this is an 'off' pixel
				continue;
			}
			
			neighbour_found = 0;
			for (unsigned long l = j - 1; l <= j + 1; l++) {
				for (unsigned long k = i - 1; k <= i + 1; k++) {
					if ((k == i) && (l == j)) {
						continue;	// this is the pixel that we are considering itself, not the environment
					}
					
					if (thresholded_image->get(k, l) > 128) {
						neighbour_found = 1;
						break;
					}
				}
			}
			
			if (neighbour_found == 0) {
				// we didn't find an active point in the neighborhood, it was an isolated pixel
				processed_thresholded_image->set(i, j, 0);
			}
		}
	}
	
	return processed_thresholded_image;
}
				


boost::shared_ptr<encap_gsl_matrix> convolve_matrices_using_fft(boost::shared_ptr<encap_gsl_matrix> image1, boost::shared_ptr<encap_gsl_matrix> image2) {
	size_t x_size1, y_size1, x_size2, y_size2;
	size_t half_index;
	
	x_size1 = image1->get_x_size();
	y_size1 = image1->get_y_size();
	x_size2 = image2->get_x_size();
	y_size2 = image2->get_y_size();
	
	size_t n_pixels, offset;
	size_t n_FFT_values;
	
	double value;
	
	double *array1;
	double *array2;
	fftw_complex *array1_FFT;
	fftw_complex *array2_FFT;
	fftw_complex complex_value;
	fftw_plan transform_plan;
	boost::shared_ptr<encap_gsl_matrix> convolved_image;
	boost::shared_ptr<encap_gsl_matrix> convolved_image_aligned;
	
	// are the dimensions equal?
	if ((x_size1 != x_size2) || (y_size1 != y_size2)) {
		XOPNotice("Error in 'convolve_matrices_using_fft': dimensions of the matrices are not equal.");
		throw INCOMPATIBLE_DIMENSIONING;
	}
	
	n_pixels = x_size1 * y_size1;
	double normalization_factor = (double)(n_pixels);
	
	// convert both of the matrices to an array of doubles
	// we allocate the memory using fftw_malloc()
	// as this is recommended by the library
	array1 = (double *)fftw_malloc(sizeof(double) * n_pixels);
	if (array1 == NULL) {
		string error;
		error = "unable to allocate array1 in convolve_matrices_using_fft()\r";
		throw OUT_OF_MEMORY(error);
	}
	
	array2 = (double *)fftw_malloc(sizeof(double) * n_pixels);
	if (array2 == NULL) {
		fftw_free(array1);
		string error;
		error = "unable to allocate array2 in convolve_matrices_using_fft()\r";
		throw OUT_OF_MEMORY(error);
	}
	
	offset = 0;
	// now copy the matrices to the arrays
	
	for (size_t i = 0; i < x_size1; i++) {
		for (size_t j = 0; j < y_size1; j++) {
			// IMPORTANT: the data in the array is assumed to be in ROW-MAJOR order, so we loop over y first
			array1[offset] = image1->get(i, j);
			array2[offset] = image2->get(i, j);
			
			offset++;
		}
	}
	
	// now allocate the arrays that will hold the transformed result
	// the dimensions of these arrays are a bit unusual, and are x_size * (y_size / 2 + 1)
	n_FFT_values = x_size1 * (y_size1 / 2 + 1);
	
	array1_FFT = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * n_FFT_values);
	if (array1_FFT == NULL) {
		fftw_free(array1);
		fftw_free(array2);
		string error;
		error = "unable to allocate array1_FFT in convolve_matrices_using_fft()\r";
		throw OUT_OF_MEMORY(error);
	}
	array2_FFT = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * n_FFT_values);
	if (array2_FFT == NULL) {
		fftw_free(array1);
		fftw_free(array2);
		fftw_free(array1_FFT);
		string error;
		error = "unable to allocate array2_FFT in convolve_matrices_using_fft()\r";
		throw OUT_OF_MEMORY(error);
	}
	
	// prepare the transform and execute it on the first array
	transform_plan = fftw_plan_dft_r2c_2d((int)(x_size1), (int)(y_size1), array1, array1_FFT, FFTW_ESTIMATE);
	fftw_execute(transform_plan);
	
	fftw_destroy_plan(transform_plan);
	
	// do the same on the second array
	transform_plan = fftw_plan_dft_r2c_2d((int)(x_size1), (int)(y_size1), array2, array2_FFT, FFTW_ESTIMATE);
	fftw_execute(transform_plan);
	
	fftw_destroy_plan(transform_plan);
	
	// now do the convolution
	for (size_t i = 0; i < n_FFT_values; i++) {
		complex_value[0] = array1_FFT[i][0] * array2_FFT[i][0] - array1_FFT[i][1] * array2_FFT[i][1];
		complex_value[1] = array1_FFT[i][0] * array2_FFT[i][1] + array1_FFT[i][1] * array2_FFT[i][0];
		
		// store the result in the first array
		array1_FFT[i][0] = complex_value[0];
		array1_FFT[i][1] = complex_value[1];
	}
	
	// now do the reverse transform
	// we overwrite the original array
	transform_plan = fftw_plan_dft_c2r_2d((int)(x_size1), (int)(y_size1), array1_FFT, array1, FFTW_ESTIMATE);
	fftw_execute(transform_plan);
	fftw_destroy_plan(transform_plan);
	
	// and store the result back in a new gsl_matrix (we don't overwrite the input arguments)
	try {
		convolved_image = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(x_size1, y_size1));
	}
	catch(std::bad_alloc) {
		fftw_free(array1);
		fftw_free(array2);
		fftw_free(array1_FFT);
		fftw_free(array2_FFT);
		string error;
		error = "unable to allocate convolved_image in convolve_matrices_using_fft()\r";
		throw OUT_OF_MEMORY(error);
	}
	
	offset = 0;
	for (size_t i = 0; i < x_size1; i++) {
		for (size_t j = 0; j < y_size1; j++) {
			// the data in the array is assumed to be in ROW-MAJOR order, so we loop over x first
			// we also normalize the result
			convolved_image->set(i, j, (array1[offset] / normalization_factor));
			
			offset++;
		}
	}
	
	// cleanup
	fftw_free(array1);
	fftw_free(array2);
	fftw_free(array1_FFT);
	fftw_free(array2_FFT);
	
	// now shift the output matrix so that the image aligns properly, and is not shifted into quadrants
	convolved_image_aligned = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(x_size1, y_size1));
	
	// start by shifting the top and bottom half
	half_index = y_size1 / 2;
	for (unsigned long j = 0; j < y_size1; j++) {
		for (unsigned long i = 0; i < x_size1; i++) {
			value = convolved_image->get(i, ((j + half_index) % y_size1));
			convolved_image_aligned->set(i, j, value);
		}
	}
	
	// now copy the values
	// convolved_image->copy(*convolved_image_aligned);
	for (unsigned long j = 0; j < y_size1; j++) {
		for (unsigned long i = 0; i < x_size1; i++) {
			convolved_image->set(i, j, convolved_image_aligned->get(i,j));
		}
	}
	
	// now shift the left and right part
	half_index = x_size1 / 2;
	for (unsigned long j = 0; j < y_size1; j++) {
		for (unsigned long i = 0; i < x_size1; i++) {
			value = convolved_image->get(((i + half_index) % x_size1), j);
			convolved_image_aligned->set(i, j, value);
		}
	}
	
	return convolved_image_aligned;
}
	
	
	

int Gauss_2D_fit_function(const gsl_vector *params, void *measured_intensities_struct, gsl_vector *model_values) {
	
	measured_data_Gauss_fits *intensities_local = (measured_data_Gauss_fits *)measured_intensities_struct;
	
	unsigned long number_of_intensities = intensities_local->get_number_of_intensities();
	double *measured_intensities = intensities_local->get_intensities();
	double *sigma = intensities_local->get_sigma();
	unsigned long x_size = intensities_local->get_x_size();
	unsigned long y_size = intensities_local->get_y_size();
	double x_offset = intensities_local->get_x_offset();
	double y_offset = intensities_local->get_y_offset();
	
	double amplitude = gsl_vector_get(params, 0);
	double r = gsl_vector_get(params, 1);
	double x0 = gsl_vector_get(params, 2);
	double y0 = gsl_vector_get(params, 3);
	double offset = gsl_vector_get(params, 4);
	
	double x,y, function_value, square_deviation;
	
	if (r == 0) {
		return GSL_FAILURE;
	}
	
	for (unsigned long i = 0; i < number_of_intensities; i++) {
		return_xy_from_array(i, x_size, y_size, x_offset, y_offset, x, y);
		
		function_value = offset + amplitude * exp(- (((x0 - x)/ r) * ((x0 - x)/ r) + ((y0 - y) / r) * ((y0 - y) / r)));
		square_deviation = (function_value - measured_intensities[i]) / sigma[i];
		gsl_vector_set(model_values, i, square_deviation);
	}
	return GSL_SUCCESS;
}

int Gauss_2D_fit_function_Jacobian(const gsl_vector *params, void *measured_intensities, gsl_matrix *jacobian) {
	
	measured_data_Gauss_fits *intensities_local = (measured_data_Gauss_fits *)measured_intensities;
	
	unsigned long number_of_intensities = intensities_local->get_number_of_intensities();
	//	double *measured_intensities = intensities_local->get_intensities();
	double *sigma = intensities_local->get_sigma();
	unsigned long x_size = intensities_local->get_x_size();
	unsigned long y_size = intensities_local->get_y_size();
	double x_offset = intensities_local->get_x_offset();
	double y_offset = intensities_local->get_y_offset();
	
	double amplitude = gsl_vector_get(params, 0);
	double r = gsl_vector_get(params, 1);
	double x0 = gsl_vector_get(params, 2);
	double y0 = gsl_vector_get(params, 3);
	//	double offset = gsl_vector_get(params, 4);
	
	double x,y, exp_factor;
	double dfdA, dfdr, dfdx0, dfdy0,dfdoffset;
	
	if (r == 0) {
		return GSL_FAILURE;
	}
	
	for (unsigned long i = 0; i < number_of_intensities; i++) {
		return_xy_from_array(i, x_size, y_size, x_offset, y_offset, x, y);
		
		exp_factor = exp(- ((x0 - x)/ r) * ((x0 - x)/ r) - ((y0 - y) / r) * ((y0 - y) / r));
		
		dfdA = exp_factor / sigma[i];
		dfdr = (2 * (y - y0) * (y - y0) / r / r / r + 2 * (x - x0) * (x - x0) / r / r / r) * exp_factor * amplitude / sigma[i];
		dfdx0 = (2 * (x - x0) * exp_factor * amplitude) / (r * r *sigma[i]);
		dfdy0 = (2 * (y - y0) * exp_factor * amplitude) / (r * r *sigma[i]);
		dfdoffset = 1/sigma[i];
		
		gsl_matrix_set(jacobian, i, 0, dfdA);
		gsl_matrix_set(jacobian, i, 1, dfdr);
		gsl_matrix_set(jacobian, i, 2, dfdx0);
		gsl_matrix_set(jacobian, i, 3, dfdy0);
		gsl_matrix_set(jacobian, i, 4, dfdoffset);
	}
	
	return GSL_SUCCESS;
	
}

int Gauss_2D_fit_function_and_Jacobian(const gsl_vector *params, void *measured_intensities_struct, gsl_vector *model_values, gsl_matrix *jacobian) {
	
	Gauss_2D_fit_function(params, measured_intensities_struct, model_values);
	Gauss_2D_fit_function_Jacobian(params, measured_intensities_struct, jacobian);
	
	return GSL_SUCCESS;
}

int Gauss_2D_Poissonian_fit_function(const gsl_vector *params, void *measured_intensities_struct, gsl_vector *model_values) {
	measured_data_Gauss_fits *intensities_local = (measured_data_Gauss_fits *)measured_intensities_struct;
	
	unsigned long number_of_intensities = intensities_local->get_number_of_intensities();
	double *measured_intensities = intensities_local->get_intensities();
	// double *sigma = intensities_local->get_sigma();
	unsigned long x_size = intensities_local->get_x_size();
	unsigned long y_size = intensities_local->get_y_size();
	double x_offset = intensities_local->get_x_offset();
	double y_offset = intensities_local->get_y_offset();
	
	double amplitude = gsl_vector_get(params, 0);
	double r = gsl_vector_get(params, 1);
	double x0 = gsl_vector_get(params, 2);
	double y0 = gsl_vector_get(params, 3);
	double offset = gsl_vector_get(params, 4);
	
	double x,y, function_value, square_deviation;
	
	if (r == 0) {
		return GSL_FAILURE;
	}
	
	for (unsigned long i = 0; i < number_of_intensities; i++) {
		return_xy_from_array(i, x_size, y_size, x_offset, y_offset, x, y);
		
		function_value = offset + amplitude * exp(- (((x0 - x)/ r) * ((x0 - x)/ r) + ((y0 - y) / r) * ((y0 - y) / r)));
		// square_deviation = (function_value - measured_intensities[i]) / sigma[i];
		square_deviation = sqrt(2 * measured_intensities[i] * log(measured_intensities[i] / function_value));
		gsl_vector_set(model_values, i, square_deviation);
	}
	return GSL_SUCCESS;
}

int Gauss_2D_Poissonian_fit_function_Jacobian(const gsl_vector *params, void *measured_intensities, gsl_matrix *jacobian) {
	
	measured_data_Gauss_fits *intensities_local = (measured_data_Gauss_fits *)measured_intensities;
	
	unsigned long number_of_intensities = intensities_local->get_number_of_intensities();
	// double *measured_intensities = intensities_local->get_intensities();
	// double *sigma = intensities_local->get_sigma();
	unsigned long x_size = intensities_local->get_x_size();
	unsigned long y_size = intensities_local->get_y_size();
	double x_offset = intensities_local->get_x_offset();
	double y_offset = intensities_local->get_y_offset();
	double *measured_intensities_array = intensities_local->get_intensities();
	
	double amplitude = gsl_vector_get(params, 0);
	double r = gsl_vector_get(params, 1);
	double x0 = gsl_vector_get(params, 2);
	double y0 = gsl_vector_get(params, 3);
	double offset = gsl_vector_get(params, 4);
	
	double x,y, exp_factor;
	double dfdA, dfdr, dfdx0, dfdy0,dfdoffset;
	
	if (r == 0) {
		return GSL_FAILURE;
	}
	
	// the maxima code to get the expression to derive from:
	// sqrt(2 * yi * log(yi / (offset + A * exp(-(((x0 - x) / r)^2 + ((y0 - y) / r)^2)))))
	
	double sqrt_2 = 1.414213562373095;
	double measured_intensity;
	double denominator;
	
	for (unsigned long i = 0; i < number_of_intensities; i++) {
		measured_intensity = measured_intensities_array[i];
		return_xy_from_array(i, x_size, y_size, x_offset, y_offset, x, y);
		
		exp_factor = exp(- ((x0 - x)/ r) * ((x0 - x)/ r) - ((y0 - y) / r) * ((y0 - y) / r));
		denominator = (exp_factor * amplitude  + offset) * sqrt(measured_intensity * log(measured_intensity / (exp_factor * amplitude + offset)));
		
		dfdA = - (sqrt_2 * exp_factor * measured_intensity) / (2.0 * denominator);
		dfdr = - sqrt_2 * (2.0 * (y - y0) * (y - y0) / r / r / r + 2 * (x - x0) * (x - x0) / r / r / r) * exp_factor * amplitude * measured_intensity / (2 * denominator);
		dfdx0 = (sqrt_2 * (x0 - x) * exp_factor * amplitude * measured_intensity) / (r * r * denominator);
		dfdy0 = (sqrt_2 * (y0 - y) * exp_factor * amplitude * measured_intensity) / (r * r * denominator);
		dfdoffset = (- sqrt_2 * measured_intensity) / (2.0 * denominator);
		
		gsl_matrix_set(jacobian, i, 0, dfdA);
		gsl_matrix_set(jacobian, i, 1, dfdr);
		gsl_matrix_set(jacobian, i, 2, dfdx0);
		gsl_matrix_set(jacobian, i, 3, dfdy0);
		gsl_matrix_set(jacobian, i, 4, dfdoffset);
	}
	
	return GSL_SUCCESS;
	
}

int Gauss_2D_Poissonian_fit_function_and_Jacobian(const gsl_vector *params, void *measured_intensities_struct, gsl_vector *model_values, gsl_matrix *jacobian) {
	
	Gauss_2D_Poissonian_fit_function(params, measured_intensities_struct, model_values);
	Gauss_2D_Poissonian_fit_function_Jacobian(params, measured_intensities_struct, jacobian);
	
	return GSL_SUCCESS;
}

gsl_vector * convert_gsl_matrix_to_vector(gsl_matrix *matrix) {
	unsigned long number_of_rows, number_of_columns;
	unsigned long number_of_elements;
	unsigned long offset = 0;
	
	number_of_rows = matrix->size1;
	number_of_columns = matrix->size2;
	number_of_elements = number_of_rows * number_of_columns;
	gsl_vector *vector = gsl_vector_alloc(number_of_elements);
	if (vector == NULL) {
		return NULL;
	}
	
	for (unsigned long j = 0; j < number_of_columns; j++) {
		for (unsigned long i = 0; i < number_of_rows; i++) {
			gsl_vector_set(vector, offset, gsl_matrix_get(matrix, i, j));
			offset++;
		}
	}
	return vector;
}	


double * convert_gsl_matrix_to_array(gsl_matrix *matrix, unsigned long &number_of_elements) {
	unsigned long number_of_rows, number_of_columns;
	unsigned long offset = 0;
	
	number_of_rows = matrix->size1;
	number_of_columns = matrix->size2;
	number_of_elements = number_of_rows * number_of_columns;
	double *array = new double[number_of_elements];
	if (array == NULL)
		return NULL;	// we will handle this condition in the calling function
	
	for (unsigned long j = 0; j < number_of_rows; j++) {
		for (unsigned long i = 0; i < number_of_columns; i++) {
			array[offset] = gsl_matrix_get(matrix, i, j);
			offset++;
		}
	}
	return array;
}

inline int return_xy_from_array(unsigned long pos, unsigned long x_size, unsigned long y_size, double x_offset, double y_offset, double &x, double &y) {
	unsigned long number_of_pixels = x_size * y_size;
	
	if (pos > number_of_pixels)
		return -1;
	
	y = (double)(pos / x_size); // integer division
	x = (double)(pos % x_size); // modulus operator
	
	x += x_offset;
	y += y_offset;
	
	return 0;
}


boost::shared_ptr<encap_gsl_matrix> FitPositionsGaussian::fit_positions(const boost::shared_ptr<encap_gsl_matrix> image, boost::shared_ptr<encap_gsl_matrix> positions) {
	
	unsigned long startPosition, endPosition;
	boost::shared_ptr<encap_gsl_matrix> fittedPositions;
	
	startPosition = 0;
	endPosition = positions->get_x_size() - 1;
	
	fittedPositions = fit_positions(image, positions, startPosition, endPosition);
	
	return fittedPositions;
	
}

boost::shared_ptr<encap_gsl_matrix> FitPositionsGaussian::fit_positions(const boost::shared_ptr<encap_gsl_matrix> image, boost::shared_ptr<encap_gsl_matrix> positions, 
																		unsigned long startPos, unsigned long endPos) {
	
	// some safety checks
	if ((endPos >= positions->get_x_size()) || (startPos >= positions->get_x_size())) {
		string error;
		error = "Requested start and end positions are outside the range of positions supplied in FitPositionsGaussian::fit_positions";
		throw std::range_error(error);
	}
	
	if (startPos > endPos) {
		string error;
		error = "Start is beyond end in FitPositionsGaussian::fit_positions";
		throw std::range_error(error);
	}
	
	unsigned long number_of_positions = endPos - startPos + 1;
	unsigned long size_of_subset = 2 * cutoff_radius + 1;
	unsigned long x_offset, y_offset, x_max, y_max;
	unsigned long number_of_intensities = size_of_subset * size_of_subset;
	unsigned long xSize = image->get_x_size();
	unsigned long ySize = image->get_y_size();
	
	double x0_initial, y0_initial, initial_intensity, amplitude;
	double chi, degreesOfFreedom, c;
	long iterations = 0;
	int status;
	
	boost::shared_ptr<encap_gsl_matrix> image_subset;
	boost::shared_ptr<encap_gsl_matrix> fitted_positions;
	
	image_subset = boost::shared_ptr<encap_gsl_matrix> (new encap_gsl_matrix(size_of_subset, size_of_subset));
	fitted_positions = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(number_of_positions, N_OUTPUT_PARAMS_PER_FITTED_POSITION));
	
	fitted_positions->set_all(0);
	
	// initialize the solver
	const gsl_multifit_fdfsolver_type *solver;
	gsl_multifit_fdfsolver *fit_iterator;
	
	solver = gsl_multifit_fdfsolver_lmsder;
	measured_data_Gauss_fits measured_data;
	gsl_multifit_function_fdf f;
	
	gsl_vector *fit_parameters = gsl_vector_alloc(5);
	if (fit_parameters == NULL) {
		string error;
		error = "unable to allocate fit_parameters in FitPositionsGaussian::fit_positions()\r";
		throw OUT_OF_MEMORY(error);
	}
	
	fit_iterator = gsl_multifit_fdfsolver_alloc(solver, number_of_intensities, 5);
	if (fit_iterator == NULL) {
		gsl_vector_free(fit_parameters);
		string error;
		error = "unable to allocate fit_iterator in FitPositionsGaussian::fit_positions()\r";
		throw OUT_OF_MEMORY(error);
	}
	
	gsl_matrix *covarianceMatrix = gsl_matrix_alloc(5, 5);
	if (covarianceMatrix == NULL) {
		gsl_vector_free(fit_parameters);
		gsl_multifit_fdfsolver_free(fit_iterator);
		string error;
		error = "unable to allocate covarianceMatrix in FitPositionsGaussian::fit_positions()\r";
		throw OUT_OF_MEMORY(error);
	}
	
	f.f = &Gauss_2D_fit_function;
	f.df = &Gauss_2D_fit_function_Jacobian;
	f.fdf = &Gauss_2D_fit_function_and_Jacobian;
	f.n = number_of_intensities;
	f.p = 5;
	f.params = (void *)&measured_data;
	
	
	// iterate over all the determined positions
	for (unsigned long i = startPos; i <= endPos; i++) {
		iterations = 0;
		
		initial_intensity = positions->get(i, 0);
		x0_initial = positions->get(i, 1);
		y0_initial = positions->get(i, 2);
		
		amplitude = initial_intensity - background;
		
		x_offset = x0_initial - cutoff_radius;
		y_offset = y0_initial - cutoff_radius;
		x_max = x0_initial + cutoff_radius;
		y_max = y0_initial + cutoff_radius;
		
		if ((x_offset < 0) || (x_max > (xSize - 1)) || (y_offset < 0) || (y_max > (ySize - 1))) {	// this positions is too close to the edge of the image
																									// we cannot include it
			continue;
		}
		
		for (unsigned long k = y_offset; k <= y_max; k++) {
			for (unsigned long j = x_offset; j <= x_max; j++) {
				image_subset->set(j - x_offset, k - y_offset, image->get(j, k));
			}
		}
		
		// ironically, for the fitting we need to provide the data as an array, so we'll convert again
		double *intensity_array = convert_gsl_matrix_to_array(image_subset->get_ptr(), number_of_intensities);
		if (intensity_array == NULL) {
			gsl_vector_free(fit_parameters);
			gsl_multifit_fdfsolver_free(fit_iterator);
			gsl_matrix_free(covarianceMatrix);
			string error;
			error = "unable to allocate intensity_array in FitPositionsGaussian::fit_positions()\r";
			throw OUT_OF_MEMORY(error);
		}
		
		double *sigma_array = new double[number_of_intensities];
		if (sigma_array == NULL) {
			gsl_vector_free(fit_parameters);
			gsl_multifit_fdfsolver_free(fit_iterator);
			gsl_matrix_free(covarianceMatrix);
			delete[] intensity_array;
			string error;
			error = "unable to allocate sigma_array in FitPositionsGaussian::fit_positions()\r";
			throw OUT_OF_MEMORY(error);
		}
		
		for (unsigned long j = 0; j < number_of_intensities; j++) {
			sigma_array[j] = sigma;
		}
		
		measured_data.set_number_of_intensities(number_of_intensities);
		measured_data.set_intensities(intensity_array);
		measured_data.set_sigma(sigma_array);
		measured_data.set_x_size(size_of_subset);
		measured_data.set_y_size(size_of_subset);
		measured_data.set_x_offset(x_offset);
		measured_data.set_y_offset(y_offset);
		
		// provide the initial parameters
		gsl_vector_set(fit_parameters, 0, amplitude);
		gsl_vector_set(fit_parameters, 1, r_initial * 1.414213562373095);	// because the fitting function is of the form 1/r^2, but standard deviation is 1/(2 r^2), we have to correct by sqrt(2)
		gsl_vector_set(fit_parameters, 2, x0_initial);
		gsl_vector_set(fit_parameters, 3, y0_initial);
		gsl_vector_set(fit_parameters, 4, background);
		
		// set the solver
		gsl_multifit_fdfsolver_set(fit_iterator, &f, fit_parameters);
		
		// run the iterations
		do {
			iterations++;
			status = gsl_multifit_fdfsolver_iterate(fit_iterator);
			if (status != 0)
				break;
			//			print_fit_state(iterations, fit_iterator);
			status = gsl_multifit_test_delta(fit_iterator->dx, fit_iterator->x, 10, 10);
		} while ((status = GSL_CONTINUE) && (iterations < 200));
		
		if ((status != GSL_SUCCESS) && (iterations == 200)) {
			// max number of iterations reached
		}
		if ((status != GSL_SUCCESS) && (iterations < 200)) {
			// some error occurred
			//			cout << gsl_strerror(status) << "\n";
		}
		
		// calculate the covariance matrix
		gsl_multifit_covar(fit_iterator->J, 0.0, covarianceMatrix);
		chi = gsl_blas_dnrm2(fit_iterator->f);
		degreesOfFreedom = (2 * cutoff_radius - 1) * (2 * cutoff_radius - 1) - 5;
		c = GSL_MAX_DBL(1, chi / sqrt(degreesOfFreedom));
		
		// store the data
		for (unsigned long j = 0; j < 5; ++j) {
			fitted_positions->set(i - startPos, j, gsl_vector_get(fit_iterator->x, j));
		}
		
		// store the errors
		for (unsigned long j = 5; j < 10; ++j) {
			fitted_positions->set(i - startPos, j, c * sqrt(gsl_matrix_get(covarianceMatrix, j - 5, j - 5)));
		}
		
		// store the number of iterations
		if ((status == GSL_SUCCESS) || (status == GSL_ETOLF) || (status == GSL_ETOLX)) {
			fitted_positions->set(i - startPos, 10, (double)iterations);
		} else {
			fitted_positions->set(i - startPos, 10, (double)(-1 * status));
		}
		
		// the width returned by the fit function is not equal to the standard deviation (a factor of sqrt 2 is missing)
		// so we correct for that
		
		fitted_positions->set(i - startPos, 1, fitted_positions->get(i - startPos, 1) / 1.414213562373095);
		fitted_positions->set(i - startPos, 6, fitted_positions->get(i - startPos, 6) / 1.414213562373095);	// the same for the error
		
		delete[] intensity_array;
		delete[] sigma_array;
	}
	
	gsl_multifit_fdfsolver_free(fit_iterator);
	gsl_vector_free(fit_parameters);
	
	return fitted_positions;
	
}


boost::shared_ptr<encap_gsl_matrix> FitPositionsMultiplication::fit_positions(const boost::shared_ptr<encap_gsl_matrix> image, boost::shared_ptr<encap_gsl_matrix> positions) {
	unsigned long startPosition, endPosition;
	boost::shared_ptr<encap_gsl_matrix> fittedPositions;
	
	startPosition = 0;
	endPosition = positions->get_x_size() - 1;
	
	fittedPositions = fit_positions(image, positions, startPosition, endPosition);
	
	return fittedPositions;
}
		

boost::shared_ptr<encap_gsl_matrix> FitPositionsMultiplication::fit_positions(const boost::shared_ptr<encap_gsl_matrix> image, boost::shared_ptr<encap_gsl_matrix> positions, 
																			  unsigned long startPos, unsigned long endPos) {
	
	// some safety checks
	if ((endPos >= positions->get_x_size()) || (startPos >= positions->get_x_size())) {
		string error;
		error = "Requested start and end positions are outside the range of positions supplied in FitPositionsGaussian::fit_positions";
		throw std::range_error(error);
	}
	
	if (startPos > endPos) {
		string error;
		error = "Start is beyond end in FitPositionsGaussian::fit_positions";
		throw std::range_error(error);
	}
	
	unsigned long number_of_positions = endPos - startPos + 1;
	unsigned long size_of_subset = 2 * cutoff_radius + 1;
	unsigned long xSize = image->get_x_size();
	unsigned long ySize = image->get_y_size();
	unsigned long x_offset, y_offset, x_max, y_max;
	
	double x0_initial, y0_initial, initial_intensity, amplitude;
	unsigned long iterations = 0;
	
	double convergence_treshold_squared = convergence_threshold * convergence_threshold;
	double delta_squared = 10 * convergence_treshold_squared;	// this test the convergence of the position determined by the iteration
	// it is the distance between (xn-1, yn-1) and (xn, yn)
	// we initialize it to a value well over the treshold so that we will run at least two iterations
	double previous_position_x, previous_position_y;
	double current_x, current_y;
	
	boost::shared_ptr<encap_gsl_matrix> image_subset;
	boost::shared_ptr<encap_gsl_matrix> image_subset_mask;
	boost::shared_ptr<encap_gsl_matrix> fitted_positions;
	
	image_subset = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(size_of_subset, size_of_subset));
	image_subset_mask = boost::shared_ptr<encap_gsl_matrix> (new encap_gsl_matrix(size_of_subset, size_of_subset));
	fitted_positions = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(number_of_positions, N_OUTPUT_PARAMS_PER_FITTED_POSITION));
	
	fitted_positions->set_all(0);
	
	for (unsigned long i = startPos; i <= endPos; ++i) {
		initial_intensity = positions->get(i, 0);
		x0_initial = positions->get(i, 1);
		y0_initial = positions->get(i, 2);
		
		amplitude = initial_intensity - background;
		
		x_offset = (unsigned long)x0_initial - cutoff_radius;
		y_offset = (unsigned long)y0_initial - cutoff_radius;
		x_max = (unsigned long)x0_initial + cutoff_radius;
		y_max = (unsigned long)y0_initial + cutoff_radius;
		
		if ((x_offset < 0) || (x_max > (xSize - 1)) || (y_offset < 0) || (y_max > (ySize - 1))) {	// this position is too close to the edge of the image, we cannot include it
			continue;
		}
		
		for (unsigned long k = y_offset; k <= y_max; ++k) {
			for (unsigned long j = x_offset; j <= x_max; ++j) {
				image_subset->set(j - x_offset, k - y_offset, image->get(j, k));
			}
		}
		
		iterations = 0;
		
		current_x = x0_initial - (double)x_offset;	// correct the x- and y-values for the fact that we analyze in a subset of the image rather than the complete frame
		current_y = y0_initial - (double)y_offset;
		
		while (delta_squared > convergence_treshold_squared) {
			previous_position_x = current_x;
			previous_position_y = current_y;
			
			++iterations;
			
			if (iterations > 100) {	// the multiplication is not converging, we should stop
				SetNaN64(&current_x);
				SetNaN64(&current_y);
				break;
			}
			
			multiply_with_gaussian(image_subset, image_subset_mask, current_x, current_y, r_initial, background, amplitude);
			determine_x_y_position(image_subset_mask, current_x, current_y);
			
			if (iterations == 1)	// this is the first iteration, we should not check for termination
				continue;
			
			delta_squared = (current_x - previous_position_x) * (current_x - previous_position_x) + (current_y - previous_position_y) * (current_y - previous_position_y);
		}
		
		delta_squared = 10 * convergence_treshold_squared;
		
		fitted_positions->set(i - startPos, 2, current_x + (double)x_offset);
		fitted_positions->set(i - startPos, 3, current_y + (double)y_offset);
		fitted_positions->set(i - startPos, 10, (double)iterations);
	}
	
	return fitted_positions;
}
		

int FitPositionsMultiplication::multiply_with_gaussian(boost::shared_ptr<encap_gsl_matrix> original_image, boost::shared_ptr<encap_gsl_matrix> masked_image, double x, double y, 
													   double std_dev, double background, double amplitude) {
	// we will replace the contents of masked_image with the multiplication of original_image and a gaussian centered at position (x,y)
	
	unsigned long x_size = masked_image->get_x_size();
	unsigned long y_size = masked_image->get_y_size();
	
	double gaussian_value, distance_squared;
	
	if ((original_image->get_x_size() != x_size) || (original_image->get_y_size() != y_size)) {
		throw DIMENSIONS_SHOULD_BE_EQUAL();
	}
	
	for (unsigned long j = 0; j < y_size; j++) {
		for (unsigned long i = 0; i < x_size; i++) {
			distance_squared = (x - (double)i) * (x - (double)i) + (y - (double)j) * (y - (double)j);
			
			gaussian_value = amplitude * exp(- distance_squared / (2 * std_dev * std_dev)) + background;
			
			masked_image->set(i, j, gaussian_value * original_image->get(i, j));
		}
	}
	
	return 0;
}


int FitPositionsMultiplication::determine_x_y_position(boost::shared_ptr<encap_gsl_matrix> masked_image, double &x, double &y) {
	// based on eq (3) in Thompson Biophys J 2002
	
	unsigned long x_size = (unsigned long)masked_image->get_x_size();
	unsigned long y_size = (unsigned long)masked_image->get_y_size();
	
	double numerator_x = 0, denominator = 0;
	double numerator_y = 0;
	
	// start with determining the x-position
	for (unsigned long j = 0; j < y_size; j++) {
		for (unsigned long i = 0; i < x_size; i++) {
			numerator_x += (double)i * masked_image->get(i, j);
			numerator_y += (double)j * masked_image->get(i, j);
			denominator += masked_image->get(i, j);
		}
	}
	
	x = numerator_x / denominator;
	
	y = numerator_y / denominator;
	
	return 0;
}


boost::shared_ptr<encap_gsl_matrix> FitPositionsCentroid::fit_positions(const boost::shared_ptr<encap_gsl_matrix> image, boost::shared_ptr<encap_gsl_matrix> positions) {
	unsigned long startPosition, endPosition;
	boost::shared_ptr<encap_gsl_matrix> fittedPositions;
	
	startPosition = 0;
	endPosition = positions->get_x_size() - 1;
	
	fittedPositions = fit_positions(image, positions, startPosition, endPosition);
	
	return fittedPositions;
}


boost::shared_ptr<encap_gsl_matrix> FitPositionsCentroid::fit_positions(const boost::shared_ptr<encap_gsl_matrix> image, boost::shared_ptr<encap_gsl_matrix> positions, 
																			  unsigned long startPos, unsigned long endPos) {
	
	// some safety checks
	if ((endPos >= positions->get_x_size()) || (startPos >= positions->get_x_size())) {
		string error;
		error = "Requested start and end positions are outside the range of positions supplied in FitPositionsCentroid::fit_positions";
		throw std::range_error(error);
	}
	
	if (startPos > endPos) {
		string error;
		error = "Start is beyond end in FitPositionsCentroid::fit_positions";
		throw std::range_error(error);
	}
	
	unsigned long number_of_positions = endPos - startPos + 1;
	unsigned long size_of_subset = 2 * cutoff_radius + 1;
	unsigned long xSize = image->get_x_size();
	unsigned long ySize = image->get_y_size();
	unsigned long x_offset, y_offset, x_max, y_max;
	
	unsigned long x0_initial, y0_initial;
	double current_x, current_y;
	double denominator;
	
	boost::shared_ptr<encap_gsl_matrix> fitted_positions;
	
	fitted_positions = boost::shared_ptr<encap_gsl_matrix>(new encap_gsl_matrix(number_of_positions, N_OUTPUT_PARAMS_PER_FITTED_POSITION));
	
	fitted_positions->set_all(0);
	
	for (unsigned long i = startPos; i <= endPos; ++i) {
		x0_initial = positions->get(i, 1);
		y0_initial = positions->get(i, 2);
		current_x = 0;
		current_y = 0;
		denominator = 0;
		
		x_offset = x0_initial - cutoff_radius;
		y_offset = y0_initial - cutoff_radius;
		x_max = x0_initial + cutoff_radius;
		y_max = y0_initial + cutoff_radius;
		
		if ((x_offset > xSize) || (x_max > (xSize - 1)) || (y_offset > ySize) || (y_max > (ySize - 1))) {	// the point is too close to the edge
																											// because all the variables are unsigned, a negative value will
																											// actually end up being larger than xSize or ySize
			continue;
		}
		
		for (unsigned long j = x_offset; j <= x_max; ++j) {
			for (unsigned long k = y_offset; k <= y_max; ++k) {
				current_x += (double)j * image->get(j, k);
				current_y += (double)k * image->get(j, k);
				denominator += image->get(j, k);
			}
		}
		
		current_x /= denominator;
		current_y /= denominator;		
		
		
		fitted_positions->set(i - startPos, 2, current_x);
		fitted_positions->set(i - startPos, 3, current_y);
	}
	
	return fitted_positions;
}


OutputWriter::OutputWriter() {
	file_path.assign("");
	n_images_written = 0;
}

OutputWriter::OutputWriter(const string &rhs, int overwrite) {
	file_path.assign("");
	n_images_written = 0;
}



SimpleOutputWriter::SimpleOutputWriter(const string &rhs,int overwrite) {
	// if overwrite is non-zero then we overwrite any file that exists at the output path
	// if it is set to zero then we throw an error and abort instead of overwriting
	file_path = rhs;
	 int header_length = 3 * sizeof(unsigned long);
	
	char *header = new char[header_length];
	if (header == NULL) {
		string error;
		error = "unable to allocate header in SimpleOutputWriter::SimpleOutputWriter()\r";
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


SimpleOutputWriter::~SimpleOutputWriter() {
	
	if (file.is_open() != 0) {
		flush_and_close();
	}
}



void SimpleOutputWriter::write_image(boost::shared_ptr<encap_gsl_matrix> new_image) {
	
	// check whether we should write the cache to disk
	if (image_buffer.size() == N_SIMULTANEOUS_IMAGE_WRITES) {	// we need to flush the cache before accepting a new image
		flush_cache();
	}
	
	// add the new image to the queue
	image_buffer.push(new_image);
	n_images_written++;
}


void SimpleOutputWriter::flush_cache() {
	boost::shared_ptr<encap_gsl_matrix> current_image;
	unsigned long n_pixels, offset;
	
	if (image_buffer.size() == 0) {
		return;
	}
	
	// determine the size of the frames
	current_image = image_buffer.front();
	x_size = current_image->get_x_size();
	y_size = current_image->get_y_size();
	n_pixels = x_size * y_size;
	
	float *single_image_buffer = new float[n_pixels];	// a temporary buffer for writing a single image
	if (single_image_buffer == NULL) {
		string error;
		error = "unable to allocate buffer in SimpleOutputWriter::flush_cache()\r";
		throw OUT_OF_MEMORY(error);
	}
	
	while (image_buffer.size() != 0) {
		offset = 0;
		
		current_image = image_buffer.front();
		
		
		for (unsigned long j = 0; j < y_size; j++) {
			for (unsigned long i = 0; i < x_size; i++) {
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

int SimpleOutputWriter::flush_and_close() {
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
	file.write((char *)&x_size, sizeof(unsigned long));
	file.write((char *)&y_size, sizeof(unsigned long));
	file.write((char *)&n_images_written, sizeof(unsigned long));
	
	file.close();
	
	// we add the size of the image frame and the number of images to the top of the file
/*	fstream file_header_out;
	file_header_out.open(file_path.c_str(), ios::in | ios::out | ios::binary | ios::ate);
	file_header_out.seekp(0);
	file_header_out.seekg(0);	// should be unnecessary
	
	file_header_out.write((char *)&x_size, sizeof(unsigned long));
	file_header_out.write((char *)&y_size, sizeof(unsigned long));
	file_header_out.write((char *)&n_images_written, sizeof(unsigned long));
	
	file_header_out.close();*/
	
	return 0;
}


TIFFOutputWriter::TIFFOutputWriter(const string &rhs,int overwrite) {
	// if overwrite is non-zero then we overwrite any file that exists at the output path
	// if it is set to zero then we throw an error and abort instead of overwriting
	
	file_path = rhs;
	
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


TIFFOutputWriter::~TIFFOutputWriter() {
	
	if (tiff_file != NULL) {
		flush_and_close();
	}
}



void TIFFOutputWriter::write_image(boost::shared_ptr<encap_gsl_matrix> new_image) {
	
	// check whether we should write the cache to disk
	if (image_buffer.size() == N_SIMULTANEOUS_IMAGE_WRITES) {	// we need to flush the cache before accepting a new image
		flush_cache();
	}
	
	// add the new image to the queue
	image_buffer.push(new_image);
}


void TIFFOutputWriter::flush_cache() {
	boost::shared_ptr<encap_gsl_matrix> current_image;
	unsigned long n_pixels, offset;
	
	int result;
	uint16_t current_uint16;
	uint32_t current_uint32;
	
	if (image_buffer.size() == 0) {
		return;
	}
	
	// determine the size of the frames
	current_image = image_buffer.front();
	x_size = current_image->get_x_size();
	y_size = current_image->get_y_size();
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
		
		result = TIFFSetField(tiff_file, TIFFTAG_COMPRESSION, COMPRESSION_ADOBE_DEFLATE);
			if (result != 1) {
				string error;
				error = "Unable to set the compression method for the image at\"";
				error += file_path;
				error += "\"\r";
				throw ERROR_WRITING_FILE_DATA(error);
			}
		
		for (unsigned long j = 0; j < y_size; j++) {
			offset = 0;
			for (unsigned long i = 0; i < x_size; i++) {
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

int TIFFOutputWriter::flush_and_close() {
	if (tiff_file != NULL) {
		flush_cache();
		TIFFClose(tiff_file);
		tiff_file = NULL;
	}
	
	return 0;
}


IgorOutputWriter::IgorOutputWriter(const string rhs) {
	wave_name = rhs;
	total_number_of_positions = 0;
}

IgorOutputWriter::~IgorOutputWriter() {
	/* unsigned long number_of_matrices = positions_array.size();
	boost::shared_ptr<encap_gsl_matrix> current_array;
	
	for (unsigned long i = 0; i < number_of_matrices; i++) {
		current_array = positions_array[i];
		if (current_array != NULL)
			gsl_matrix_free(current_array);
	}*/
}

int IgorOutputWriter::append_new_positions(boost::shared_ptr<encap_gsl_matrix> positions) {
	positions_array.push_back(positions);
	if (positions != NULL) {	// if it NULL then this no positions were found
		total_number_of_positions += positions->get_x_size();
	}
	return 0;
}

int IgorOutputWriter::write_positions_to_wave() {
	waveHndl output_wave;
	long dim_size[3];
	long indices[2];
	unsigned long vector_size = positions_array.size();
	int status;
	double value;
	boost::shared_ptr<encap_gsl_matrix> positions;
	long number_of_positions_in_matrix, offset = 0;
	
	dim_size[0] = (long)total_number_of_positions;
	dim_size[1] = 12;
	dim_size[2] = 0;
	
	// the output wave will be two-dimensional
	// the first dimension contains the different individual positions
	// the second position is structured as follows:
	// n_image	amplitude	radius	x_position	y_position	offset	amp_error	radius_err	x_err	y_err	offset_err	n_iters
	
	// try to make the wave
	status = MDMakeWave(&output_wave, wave_name.c_str(), NULL, dim_size, NT_FP64, 1);
	if (status != 0)
		return status;
	
	for (long i = 0; i < vector_size; i++) {
		
		positions = positions_array[i];
		
		if (positions.get() == NULL) {
			// no positions were found in this frame
			continue;
		}
		
		number_of_positions_in_matrix = positions->get_x_size();
		
		for (long j = 0; j < number_of_positions_in_matrix; j++) {
			// start by storing the index of the matrix
			// so we can tell what image the positions correspond to
			indices[0] = offset;
			indices[1] = 0;
			value = (double)i;
			MDSetNumericWavePointValue(output_wave, indices, &value);
			
			for (long k = 0; k < 11; ++k) {
				indices[1] = k + 1;
				value = positions->get(j, k);
				MDSetNumericWavePointValue(output_wave, indices, &value);
			}
			offset++;
		}
	}
	return 0;
}

measured_data_Gauss_fits::measured_data_Gauss_fits() {
	intensities = NULL;
	sigma = NULL;
	x_y_positions = NULL;
}

measured_data_Gauss_fits::~measured_data_Gauss_fits() {
	/*	if (intensities != NULL)
	 delete[] intensities;
	 if (sigma != NULL)
	 delete[] sigma;
	 if (x_y_positions != NULL)
	 gsl_matrix_free(x_y_positions);*/
}