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
	T operator() (size_t x, size_t y) const;
	T & get(size_t x, size_t y);
	void set(size_t x, size_t y, T &value);
	void set_all(T &value);
	void set_all(T value);
	
	PALMMatrix & operator=(const PALMMatrix &rhs);
	
	PALMMatrix & operator+=(const PALMMatrix &rhs);
	PALMMatrix & operator-=(const PALMMatrix &rhs);
	PALMMatrix & operator*=(const PALMMatrix &rhs);
	PALMMatrix & operator/=(const PALMMatrix &rhs);
	
	const PALMMatrix & operator+(const PALMMatrix &rhs) const;
	const PALMMatrix & operator-(const PALMMatrix &rhs) const;
	const PALMMatrix & operator/(const PALMMatrix &rhs) const;
	const PALMMatrix & operator*(const PALMMatrix &rhs) const;
	
	size_t getXSize() const {return xSize;}
	size_t getYSize() const {return ySize;}
	
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

template <typename T> inline T PALMMatrix<T>::operator() (size_t x, size_t y) const {
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

template <typename T> inline void PALMMatrix<T>::set_all(T value) {
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
	
	return (*this);
}

template <typename T> PALMMatrix<T> & PALMMatrix<T>::operator+=(const PALMMatrix &rhs) {
	assert ((xSize == rhs.getXSize()) && (ySize == rhs.getYSize()));
	
	for (size_t i = 0; i < xSize; ++i) {
		for (size_t j = 0; j < ySize; ++j) {
			(*this)(i, j) += rhs(i, j);
		}
	}
	
	return *this;
}

template <typename T> PALMMatrix<T> & PALMMatrix<T>::operator-=(const PALMMatrix &rhs) {
	assert ((xSize == rhs.getXSize()) && (ySize == rhs.getYSize()));
	
	for (size_t i = 0; i < xSize; ++i) {
		for (size_t j = 0; j < ySize; ++j) {
			(*this)(i, j) -= rhs(i, j);
		}
	}
	
	return *this;
}

template <typename T> PALMMatrix<T> & PALMMatrix<T>::operator*=(const PALMMatrix &rhs) {
	assert ((xSize == rhs.getXSize()) && (ySize == rhs.getYSize()));
	for (size_t i = 0; i < xSize; ++i) {
		for (size_t j = 0; j < ySize; ++j) {
			(*this)(i, j) *= rhs(i, j);
		}
	}
	
	return *this;
}

template <typename T> PALMMatrix<T> & PALMMatrix<T>::operator/=(const PALMMatrix &rhs) {
	assert ((xSize == rhs.getXSize()) && (ySize == rhs.getYSize()));
	for (size_t i = 0; i < xSize; ++i) {
		for (size_t j = 0; j < ySize; ++j) {
			(*this)(i, j) /= rhs(i, j);
		}
	}
	
	return *this;
}

template <typename T> const PALMMatrix<T> & PALMMatrix<T>::operator+(const PALMMatrix &rhs) const {
	assert ((xSize == rhs.getXSize()) && (ySize == rhs.getYSize()));
	
	PALMMatrix <T> result(*this);	// make a copy of the current matrix
	result += rhs;
	
	return result;
}

template <typename T> const PALMMatrix<T> & PALMMatrix<T>::operator-(const PALMMatrix &rhs) const {
	assert ((xSize == rhs.getXSize()) && (ySize == rhs.getYSize()));
	
	PALMMatrix <T> result(*this);	// make a copy of the current matrix
	result -= rhs;
	
	return result;
}

template <typename T> const PALMMatrix<T> & PALMMatrix<T>::operator*(const PALMMatrix &rhs) const {
	assert ((xSize == rhs.getXSize()) && (ySize == rhs.getYSize()));
	
	PALMMatrix <T> result(*this);	// make a copy of the current matrix
	result *= rhs;
	
	return result;
}

template <typename T> const PALMMatrix<T> & PALMMatrix<T>::operator/(const PALMMatrix &rhs) const {
	assert ((xSize == rhs.getXSize()) && (ySize == rhs.getYSize()));
	
	PALMMatrix <T> result(*this);	// make a copy of the current matrix
	result /= rhs;
	
	return result;
}

template <typename T> class PALMVolume {
public:
	PALMVolume(size_t xSize_rhs, size_t ySize_rhs, size_t zSize_rhs);
	PALMVolume(const PALMVolume &rhs);
	~PALMVolume() {;}
	
	T & operator() (size_t x, size_t y, size_t z);
	T operator() (size_t x, size_t y, size_t z) const;
	T & get(size_t x, size_t y, size_t z);
	void set(size_t x, size_t y, size_t z, T &value);
	void set_all(T &value);
	void set_all(T value);
	
	PALMVolume & operator=(const PALMVolume &rhs);
	
	PALMVolume & operator+=(const PALMVolume &rhs);
	PALMVolume & operator-=(const PALMVolume &rhs);
	PALMVolume & operator*=(const PALMVolume &rhs);
	PALMVolume & operator/=(const PALMVolume &rhs);
	
	const PALMVolume & operator+(const PALMVolume &rhs) const;
	const PALMVolume & operator-(const PALMVolume &rhs) const;
	const PALMVolume & operator/(const PALMVolume &rhs) const;
	const PALMVolume & operator*(const PALMVolume &rhs) const;
	
	size_t getXSize() const {return xSize;}
	size_t getYSize() const {return ySize;}
	size_t getZSize() const {return zSize;}
	
protected:
	size_t xSize;
	size_t ySize;
	size_t zSize;
	vector <PALMMatrix<T> > data;
};

template <typename T> PALMVolume<T>::PALMVolume(size_t xSize_rhs, size_t ySize_rhs, size_t zSize_rhs) {
	data.resize(zSize_rhs, PALMMatrix<T>(xSize, ySize));
	xSize = xSize_rhs;
	ySize = ySize_rhs;
	zSize = zSize_rhs;
}

template <typename T> PALMVolume<T>::PALMVolume(const PALMVolume &rhs) {
	xSize = rhs.getXSize();
	ySize = rhs.getYSize();
	zSize = rhs.getZSize();
	data.resize(zSize, PALMMatrix<T>(xSize, ySize));
	
	for (size_t k = 0; k < zSize; ++k) {
		for (size_t i = 0; i < xSize; ++i) {
			for (size_t j = 0; j < ySize; ++j) {
				(*this)(i, j, k) = rhs(i, j, k);
			}
		}
	}
}

template <typename T> inline T & PALMVolume<T>::operator() (size_t x, size_t y, size_t z) {
	assert ((x < xSize) && (y < ySize) && (z < zSize));
	return data[z](x, y);
}

template <typename T> inline T PALMVolume<T>::operator() (size_t x, size_t y, size_t z) const {
	assert ((x < xSize) && (y < ySize) && (z < zSize));
	return data[z](x, y);
}

template <typename T> inline T & PALMVolume<T>::get(size_t x, size_t y, size_t z) {
	assert ((x < xSize) && (y < ySize) && (z < zSize));
	return data[z](x, y);
}

template <typename T> inline void PALMVolume<T>::set(size_t x, size_t y, size_t z, T &value) {
	assert ((x < xSize) && (y < ySize) && (z < zSize));
	data[z](x, y) = value;
}

template <typename T> inline void PALMVolume<T>::set_all(T &value) {
	for (size_t k = 0; k < zSize; ++k) {
		for (size_t i = 0; i < xSize; ++i) {
			for (size_t j = 0; j < ySize; ++j) {
				(*this)(i, j, k) = value;
			}
		}
	}
}

template <typename T> inline void PALMVolume<T>::set_all(T value) {
	for (size_t k = 0; k < zSize; ++k) {
		for (size_t i = 0; i < xSize; ++i) {
			for (size_t j = 0; j < ySize; ++j) {
				(*this)(i, j, k) = value;
			}
		}
	}
}

template <typename T> PALMVolume<T> & PALMVolume<T>::operator=(const PALMVolume &rhs) {
	if (this == &rhs) {
		return (*this);
	}
	
	data.clear();
	xSize = rhs.getXSize();
	ySize = rhs.getYSize();
	zSize = rhs.getZSize();
	data.resize(zSize, PALMMatrix<T>(xSize, ySize));
	
	for (size_t k = 0; k < zSize; ++k) {
		for (size_t i = 0; i < xSize; ++i) {
			for (size_t j = 0; j < ySize; ++j) {
				(*this)(i, j, k) = rhs(i, j, k);
			}
		}
	}
	
	return (*this);
}

template <typename T> PALMVolume<T> & PALMVolume<T>::operator+=(const PALMVolume &rhs) {
	assert ((xSize == rhs.getXSize()) && (ySize == rhs.getYSize()) && (zSize == rhs.getZSize()));
	
	for (size_t k = 0; k < zSize; ++k) {
		for (size_t i = 0; i < xSize; ++i) {
			for (size_t j = 0; j < ySize; ++j) {
				(*this)(i, j, k) += rhs(i, j, k);
			}
		}
	}
	
	return *this;
}

template <typename T> PALMVolume<T> & PALMVolume<T>::operator-=(const PALMVolume &rhs) {
	assert ((xSize == rhs.getXSize()) && (ySize == rhs.getYSize()) && (zSize == rhs.getZSize()));
	
	for (size_t k = 0; k < zSize; ++k) {
		for (size_t i = 0; i < xSize; ++i) {
			for (size_t j = 0; j < ySize; ++j) {
				(*this)(i, j, k) -= rhs(i, j, k);
			}
		}
	}
	
	return *this;
}

template <typename T> PALMVolume<T> & PALMVolume<T>::operator*=(const PALMVolume &rhs) {
	assert ((xSize == rhs.getXSize()) && (ySize == rhs.getYSize()) && (zSize == rhs.getZSize()));
	
	for (size_t k = 0; k < zSize; ++k) {
		for (size_t i = 0; i < xSize; ++i) {
			for (size_t j = 0; j < ySize; ++j) {
				(*this)(i, j, k) *= rhs(i, j, k);
			}
		}
	}
	
	return *this;
}

template <typename T> PALMVolume<T> & PALMVolume<T>::operator/=(const PALMVolume &rhs) {
	assert ((xSize == rhs.getXSize()) && (ySize == rhs.getYSize()) && (zSize == rhs.getZSize()));
	
	for (size_t k = 0; k < zSize; ++k) {
		for (size_t i = 0; i < xSize; ++i) {
			for (size_t j = 0; j < ySize; ++j) {
				(*this)(i, j, k) /= rhs(i, j, k);
			}
		}
	}
	
	return *this;
}

template <typename T> const PALMVolume<T> & PALMVolume<T>::operator+(const PALMVolume &rhs) const {
	assert ((xSize == rhs.getXSize()) && (ySize == rhs.getYSize()) && (zSize == rhs.getZSize()));
	
	PALMVolume <T> result(*this);	// make a copy of the current matrix
	result += rhs;
	
	return result;
}

template <typename T> const PALMVolume<T> & PALMVolume<T>::operator-(const PALMVolume &rhs) const {
	assert ((xSize == rhs.getXSize()) && (ySize == rhs.getYSize()) && (zSize == rhs.getZSize()));
	
	PALMVolume <T> result(*this);	// make a copy of the current matrix
	result -= rhs;
	
	return result;
}

template <typename T> const PALMVolume<T> & PALMVolume<T>::operator*(const PALMVolume &rhs) const {
	assert ((xSize == rhs.getXSize()) && (ySize == rhs.getYSize()) && (zSize == rhs.getZSize()));
	
	PALMVolume <T> result(*this);	// make a copy of the current matrix
	result *= rhs;
	
	return result;
}

template <typename T> const PALMVolume<T> & PALMVolume<T>::operator/(const PALMVolume &rhs) const {
	assert ((xSize == rhs.getXSize()) && (ySize == rhs.getYSize()) && (zSize == rhs.getZSize()));
	
	PALMVolume <T> result(*this);	// make a copy of the current matrix
	result /= rhs;
	
	return result;
}

#endif
