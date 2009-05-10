/*
 *  PALM_analysis_storage.h
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 30/03/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef PALM_ANALYSIS_STORAGE
#define PALM_ANALYSIS_STORAGE

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <string>
#include <vector>
#include "boost/smart_ptr.hpp"
#include "PALM_analysis_defines.h"
#include "PALM_analysis_errors.h"

using namespace std;

class position {
public:
	position() {x = 0; y = 0; intensity = 0;}
	position(double xLoc, double yLoc) {x = xLoc; y = yLoc;}
	position(double xLoc, double yLoc, double intensity_rhs) {x = xLoc; y = yLoc; intensity = intensity_rhs;}
	~position() {;}
	
	void set_x(double xLoc) {x = xLoc;}
	void set_y(double yLoc) {y = yLoc;}
	void set_intensity(double intensity_rhs) {intensity = intensity_rhs;}
	
	double get_x() {return x;}
	double get_y() {return y;}
	double get_intensity() {return intensity;}
	
protected:
	double x;
	double y;
	double intensity;
};


template <typename T> class PALMMatrix {
public:
	PALMMatrix(size_t xSize_rhs, size_t ySize_rhs);
	PALMMatrix(const PALMMatrix &rhs);
	~PALMMatrix();
	
	T & operator() (size_t x, size_t y);
	T & get(size_t x, size_t y);
	void set(size_t x, size_t y, T &value);
	void set_all(T &value);
	
	PALMMatrix & operator=(const PALMMatrix &rhs);
	
	const PALMMatrix & operator+(PALMMatrix &rhs);
	const PALMMatrix & operator-(PALMMatrix &rhs);
	const PALMMatrix & operator/(PALMMatrix &rhs);
	const PALMMatrix & operator*(PALMMatrix &rhs);
	
	size_t getXSize() {return xSize;}
	size_t getYSize() {return ySize;}
	
protected:
	size_t xSize;
	size_t ySize;
	T *data;
};

template <typename T> PALMMatrix<T>::PALMMatrix(size_t xSize_rhs, size_t ySize_rhs) {
	data = NULL;
	data = new T[xSize * ySize];
	xSize = xSize_rhs;
	ySize = ySize_rhs;
}

template <typename T> PALMMatrix<T>::PALMMatrix(const PALMMatrix &rhs) {
	data = NULL;
	xSize = rhs.getXSize();
	ySize = rhs.getYSize();
	
	data = new T[xSize * ySize];
	
	for (size_t i = 0; i < xSize; ++i) {
		for (size_t j = 0; j < ySize; ++j) {
			(*this)(i, j) = rhs(i, j);
		}
	}
}

template <typename T> PALMMatrix<T>::~PALMMatrix() {
	if (data != NULL) {
		delete[] data;
	}
}

template <typename T> inline T & PALMMatrix<T>::operator() (size_t x, size_t y) {
	assert ((x < xSize) && (y < ySize));
	return data[x * ySize + y];
}

template <typename T> inline T & PALMMatrix<T>::get(size_t x, size_t y) {
	assert ((x < xSize) && (y < ySize));
	return data[x * ySize + y];
}

template <typename T> inline void PALMMatrix<T>::set(size_t x, size_t y, T &value) {
	assert ((x < xSize) && (y < ySize));
	data[x * ySize + y] = value;
}

template <typename T> inline void PALMMatrix<T>::set_all(T &value) {
	size_t nItems = xSize * ySize;
	for (size_t i = 0; i < nItems; ++i) {
		data[i] = value;
	}
}

template <typename T> PALMMatrix<T> & PALMMatrix<T>::operator=(const PALMMatrix &rhs) {
	if (this == &rhs) {
		return *this;
	}
	
	delete[] data;
	xSize = rhs.getXSize();
	ySize = rhs.getYSize();
	
	for (size_t i = 0; i < xSize; ++i) {
		for (size_t j = 0; j < ySize; ++j) {
			(*this)(i, j) = rhs(i, j);
		}
	}
}

template <typename T> const PALMMatrix<T> & PALMMatrix<T>::operator+(PALMMatrix &rhs) {
	assert ((xSize == rhs.getXSize()) && (ySize = rhs.getYSize()));
	
	PALMMatrix <T> result(*this);	// make a copy of the current matrix
	for (size_t i = 0; i < xSize; ++i) {
		for (size_t j = 0; j < ySize; ++j) {
			result(i, j) = (*this)(i, j) + rhs(i, j);
		}
	}
	
	return result;
}

template <typename T> const PALMMatrix<T> & PALMMatrix<T>::operator-(PALMMatrix &rhs) {
	assert ((xSize == rhs.getXSize()) && (ySize = rhs.getYSize()));
	
	PALMMatrix <T> result(*this);	// make a copy of the current matrix
	for (size_t i = 0; i < xSize; ++i) {
		for (size_t j = 0; j < ySize; ++j) {
			result(i, j) = (*this)(i, j) - rhs(i, j);
		}
	}
	
	return result;
}

template <typename T> const PALMMatrix<T> & PALMMatrix<T>::operator*(PALMMatrix &rhs) {
	assert ((xSize == rhs.getXSize()) && (ySize = rhs.getYSize()));
	
	PALMMatrix <T> result(*this);	// make a copy of the current matrix
	for (size_t i = 0; i < xSize; ++i) {
		for (size_t j = 0; j < ySize; ++j) {
			result(i, j) = (*this)(i, j) * rhs(i, j);
		}
	}
	
	return result;
}

template <typename T> const PALMMatrix<T> & PALMMatrix<T>::operator/(PALMMatrix &rhs) {
	assert ((xSize == rhs.getXSize()) && (ySize = rhs.getYSize()));
	
	PALMMatrix <T> result(*this);	// make a copy of the current matrix
	for (size_t i = 0; i < xSize; ++i) {
		for (size_t j = 0; j < ySize; ++j) {
			result(i, j) = (*this)(i, j) / rhs(i, j);
		}
	}
	
	return result;
}


class encap_gsl_matrix {
public:
	encap_gsl_matrix() {matrix = NULL;}
	encap_gsl_matrix(size_t x, size_t y);
	~encap_gsl_matrix();
	
	gsl_matrix* get_ptr() const {return matrix;}
	
	void set(size_t x, size_t y, double value);
	void set_all(double value) {gsl_matrix_set_all(matrix, value);}
	void operator=(gsl_matrix* rhs);
	
	double get(size_t x, size_t y);
	size_t get_x_size() const {return x_size;}
	size_t get_y_size() const {return y_size;}
	
protected:
	gsl_matrix *matrix;
	size_t x_size, y_size;
};

class encap_gsl_matrix_uchar {
public:
	encap_gsl_matrix_uchar() {matrix = NULL;}
	encap_gsl_matrix_uchar(size_t x, size_t y);
	~encap_gsl_matrix_uchar();
	
	gsl_matrix_uchar* get_ptr() const {return matrix;}
	
	void set(size_t x, size_t y, unsigned char value);
	void set_all(unsigned char value) {gsl_matrix_uchar_set_all(matrix, value);}
	void operator=(gsl_matrix_uchar* rhs);
	
	unsigned char get(size_t x, size_t y);
	size_t get_x_size() const {return x_size;}
	size_t get_y_size() const {return y_size;}
	
protected:
	gsl_matrix_uchar *matrix;
	size_t x_size, y_size;
};

class encap_gsl_matrix_long {
public:
	encap_gsl_matrix_long() {matrix = NULL;}
	encap_gsl_matrix_long(size_t x, size_t y);
	~encap_gsl_matrix_long();
	
	gsl_matrix_long* get_ptr() const {return matrix;}
	
	void set(size_t x, size_t y, long value);
	void set_all(long value) {gsl_matrix_long_set_all(matrix, value);}
	void operator=(gsl_matrix_long* rhs);
	
	long get(size_t x, size_t y);
	size_t get_x_size() const {return x_size;}
	size_t get_y_size() const {return y_size;}
	
protected:
	gsl_matrix_long *matrix;
	size_t x_size, y_size;
};

class encap_gsl_volume {	// extension of an encap_gsl_matrix to three dimensions
public:
	encap_gsl_volume(size_t x, size_t y, size_t z);
	~encap_gsl_volume() {;}
	
	void set(size_t x, size_t y, size_t z, double value);
	void set_all(double value);
	
	double get(size_t x, size_t y, size_t z);
	size_t get_x_size() const {return x_size;}
	size_t get_y_size() const {return y_size;}
	size_t get_z_size() const {return z_size;}
	
protected:
	vector<boost::shared_ptr<encap_gsl_matrix> > matrices;
	size_t x_size, y_size, z_size;
};

class encap_gsl_volume_ushort {	// extension of an encap_gsl_matrix to three dimensions
public:
	encap_gsl_volume_ushort(size_t x, size_t y, size_t z);
	~encap_gsl_volume_ushort();
	
	void set(size_t x, size_t y, size_t z, unsigned short value) {assert((x < x_size) && (y < y_size) && (z < z_size)); gsl_matrix_ushort_set(matrices[z], x, y, value);}
	void set_all(unsigned short value);
	
	unsigned short get(size_t x, size_t y, size_t z) {assert((x < x_size) && (y < y_size) && (z < z_size)); return gsl_matrix_ushort_get(matrices[z], x, y);}
	size_t get_x_size() const {return x_size;}
	size_t get_y_size() const {return y_size;}
	size_t get_z_size() const {return z_size;}
	
protected:
	vector <gsl_matrix_ushort*> matrices;
	size_t x_size, y_size, z_size;
};

#endif
