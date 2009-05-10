/*
 *  PALM_analysis_storage.cpp
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 30/03/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "PALM_analysis_storage.h"


encap_gsl_matrix::encap_gsl_matrix(size_t x, size_t y) {
	matrix = gsl_matrix_alloc(x, y);
	if (matrix == NULL) {
		string error("unable to allocate matrix in encap_gsl_matrix::encap_gsl_matrix()\r");
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
	assert((x < x_size) && (y < y_size));
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
	assert((x < x_size) && (y < y_size));
	return gsl_matrix_get(matrix, x, y);
}

encap_gsl_matrix_uchar::encap_gsl_matrix_uchar(size_t x, size_t y) {
	matrix = gsl_matrix_uchar_alloc(x, y);
	if (matrix == NULL) {
		string error("unable to allocate matrix in encap_gsl_matrix::encap_gsl_matrix()\r");
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
	assert((x < x_size) && (y < y_size));
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
	assert((x < x_size) && (y < y_size));
	return gsl_matrix_uchar_get(matrix, x, y);
}

encap_gsl_matrix_long::encap_gsl_matrix_long(size_t x, size_t y) {
	matrix = gsl_matrix_long_alloc(x, y);
	if (matrix == NULL) {
		string error("unable to allocate matrix in encap_gsl_matrix::encap_gsl_matrix()\r");
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
	assert((x < x_size) && (y < y_size));
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
	assert((x < x_size) && (y < y_size));
	return gsl_matrix_long_get(matrix, x, y);
}


encap_gsl_volume::encap_gsl_volume(size_t x, size_t y, size_t z) {
	matrices.reserve(z);
	
	for (size_t i = 0; i < z; ++i) {
		matrices.push_back(boost::shared_ptr<encap_gsl_matrix> (new encap_gsl_matrix(x, y)));
	}
	
	x_size = x;
	y_size = y;
	z_size = z;
}

void encap_gsl_volume::set(size_t x, size_t y, size_t z, double value) {
	assert((x < x_size) && (y < y_size) && (z < z_size));
	matrices[z]->set(x, y, value);
}

void encap_gsl_volume::set_all(double value) {
	for (size_t i = 0; i < z_size; ++i) {
		matrices[i]->set_all(value);
	}
}

double encap_gsl_volume::get(size_t x, size_t y, size_t z) {
	assert((x < x_size) && (y < y_size) && (z < z_size));
	return matrices[z]->get(x, y);
}

encap_gsl_volume_ushort::encap_gsl_volume_ushort(size_t x, size_t y, size_t z) {
	matrices.reserve(z);
	
	for (size_t i = 0; i < z; ++i) {
		matrices.push_back(gsl_matrix_ushort_alloc(x, y));
		if (matrices[i] == NULL) {
			throw OUT_OF_MEMORY(string("Unable to alloc a gsl_matrix_ushort in encap_gsl_volume_ushort\r"));
		}
	}
	
	x_size = x;
	y_size = y;
	z_size = z;
}

encap_gsl_volume_ushort::~encap_gsl_volume_ushort() {
	for (size_t i = 0; i < z_size; ++i) {
		if (matrices[i] != NULL) {
			gsl_matrix_ushort_free(matrices[i]);
			matrices[i] = NULL;
		}
	}
}

void encap_gsl_volume_ushort::set_all(unsigned short value) {
	for (size_t i = 0; i < z_size; ++i) {
		gsl_matrix_ushort_set_all(matrices[i], value);
	}
}