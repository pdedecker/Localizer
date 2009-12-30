/*
 *  PALM_analysis_storage.h
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 30/03/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef PALM_ANALYSIS_STORAGE_H
#define PALM_ANALYSIS_STORAGE_H

#include <string>
#include <vector>
#include <cmath>
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
	
	void set_x(const double xLoc) {x = xLoc;}
	void set_y(const double yLoc) {y = yLoc;}
	void set_intensity(const double intensity_rhs) {intensity = intensity_rhs;}
	void set_background(const double background_rhs) {background = background_rhs;}
	
	double get_x() const {return x;}
	double get_y() const {return y;}
	double get_intensity() const {return intensity;}
	double get_background() const {return background;}
	
protected:
	double x;
	double y;
	double intensity;
	double background;
};


template <typename T> class PALMMatrix {
public:
	PALMMatrix(size_t xSize_rhs, size_t ySize_rhs);
	PALMMatrix(const PALMMatrix &rhs);
	~PALMMatrix();
	
	T & operator() (const size_t x, const size_t y);
	T operator() (const size_t x, const size_t y) const;
	T get(size_t x, size_t y) const;
	void set(size_t x, size_t y, const T &value);
	void set_all(const T &value);
	
	T Sum() const;
	double Average() const;
	double StandardDeviation() const;
	
	PALMMatrix operator=(const PALMMatrix &rhs);
	
	PALMMatrix & operator+=(const PALMMatrix &rhs);
	PALMMatrix & operator-=(const PALMMatrix &rhs);
	PALMMatrix & operator*=(const PALMMatrix &rhs);
	PALMMatrix & operator/=(const PALMMatrix &rhs);
	
	PALMMatrix operator+(const PALMMatrix &rhs) const;
	PALMMatrix operator-(const PALMMatrix &rhs) const;
	PALMMatrix operator/(const PALMMatrix &rhs) const;
	PALMMatrix operator*(const PALMMatrix &rhs) const;
	
	const PALMMatrix AddScalar(const double scalar);
	const PALMMatrix SubtractScalar(const double scalar);
	const PALMMatrix MultiplyWithScalar(const double scalar);
	const PALMMatrix DivideByScalar(const double scalar);
	
	const PALMMatrix RaiseToPower(const double power);
	
	size_t getXSize() const {return xSize;}
	size_t getYSize() const {return ySize;}
	
protected:
	size_t xSize;
	size_t ySize;
	T *data;
};

template <typename T> PALMMatrix<T>::PALMMatrix(size_t xSize_rhs, size_t ySize_rhs) {
	data = NULL;
	xSize = xSize_rhs;
	ySize = ySize_rhs;
	data = new T[xSize * ySize];
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

template <typename T> inline T & PALMMatrix<T>::operator() (const size_t x, const size_t y) {
	assert ((x < xSize) && (y < ySize));
	return data[x * ySize + y];
}

template <typename T> inline T PALMMatrix<T>::operator() (const size_t x, const size_t y) const {
	assert ((x < xSize) && (y < ySize));
	return data[x * ySize + y];
}

template <typename T> inline T PALMMatrix<T>::get(size_t x, size_t y) const {
	assert ((x < xSize) && (y < ySize));
	return data[x * ySize + y];
}

template <typename T> inline void PALMMatrix<T>::set(size_t x, size_t y, const T &value) {
	assert ((x < xSize) && (y < ySize));
	data[x * ySize + y] = value;
}

template <typename T> inline void PALMMatrix<T>::set_all(const T &value) {
	size_t nItems = xSize * ySize;
	
	int i;
	#pragma omp parallel for private(i) num_threads(2)
	for (i = 0; i < nItems; ++i) {
		data[i] = value;
	}
}

template <typename T> T PALMMatrix<T>::Sum() const {
	size_t nItems = xSize * ySize;
	T sum = 0;
	
	int i;
	#pragma omp parallel for private(i) num_threads(2)
	for (i = 0; i < nItems; ++i) {
		sum += data[i];
	}
}

template <typename T> double PALMMatrix<T>::Average() const {
	size_t nItems = xSize * ySize;
	double sum = 0;
	
	int i;
	#pragma omp parallel for private(i) num_threads(2)
	for (i = 0; i < nItems; ++i) {
		sum += (double)data[i];
	}
	
	sum /= (double)nItems;
	return sum;
}

template <typename T> double PALMMatrix<T>::StandardDeviation() const {
	size_t nItems = xSize * ySize;
	double average = 0;
	double standardDeviation = 0;
	
	average = this->Average();
	
	int i;
	#pragma omp parallel for private(i) num_threads(2)
	for (i = 0; i < nItems; ++i) {
		standardDeviation += ((double)data[i] - average) * ((double)data[i] - average);
	}
	
	standardDeviation /= (double)nItems;
	return sqrt(standardDeviation);
}

template <typename T> PALMMatrix<T> PALMMatrix<T>::operator=(const PALMMatrix &rhs) {
	if (this == &rhs) {
		return *this;
	}
	
	delete[] data;
	xSize = rhs.getXSize();
	ySize = rhs.getYSize();
	data = new T[xSize * ySize];
	
	int i, j;
	#pragma omp parallel for private(i, j) num_threads(2)
	for (i = 0; i < xSize; ++i) {
		for (j = 0; j < ySize; ++j) {
			(*this)(i, j) = rhs(i, j);
		}
	}
	
	return (*this);
}

template <typename T> PALMMatrix<T> & PALMMatrix<T>::operator+=(const PALMMatrix &rhs) {
	assert ((xSize == rhs.getXSize()) && (ySize == rhs.getYSize()));
	
	int i, j;
	#pragma omp parallel for private(i, j) num_threads(2)
	for (i = 0; i < xSize; ++i) {
		for (j = 0; j < ySize; ++j) {
			(*this)(i, j) += rhs(i, j);
		}
	}
	
	return *this;
}

template <typename T> PALMMatrix<T> & PALMMatrix<T>::operator-=(const PALMMatrix &rhs) {
	assert ((xSize == rhs.getXSize()) && (ySize == rhs.getYSize()));
	
	int i, j;
	#pragma omp parallel for private(i, j) num_threads(2)
	for (i = 0; i < xSize; ++i) {
		for (j = 0; j < ySize; ++j) {
			(*this)(i, j) -= rhs(i, j);
		}
	}
	
	return *this;
}

template <typename T> PALMMatrix<T> & PALMMatrix<T>::operator*=(const PALMMatrix &rhs) {
	assert ((xSize == rhs.getXSize()) && (ySize == rhs.getYSize()));
	
	int i, j;
	#pragma omp parallel for private(i, j) num_threads(2)
	for (i = 0; i < xSize; ++i) {
		for (j = 0; j < ySize; ++j) {
			(*this)(i, j) *= rhs(i, j);
		}
	}
	
	return *this;
}

template <typename T> PALMMatrix<T> & PALMMatrix<T>::operator/=(const PALMMatrix &rhs) {
	assert ((xSize == rhs.getXSize()) && (ySize == rhs.getYSize()));
	
	int i, j;
	#pragma omp parallel for private(i, j) num_threads(2)
	for (i = 0; i < xSize; ++i) {
		for (j = 0; j < ySize; ++j) {
			(*this)(i, j) /= rhs(i, j);
		}
	}
	
	return *this;
}

template <typename T> PALMMatrix<T> PALMMatrix<T>::operator+(const PALMMatrix &rhs) const {
	assert ((xSize == rhs.getXSize()) && (ySize == rhs.getYSize()));
	
	PALMMatrix <T> result(*this);	// make a copy of the current matrix
	result += rhs;
	
	return result;
}

template <typename T> PALMMatrix<T> PALMMatrix<T>::operator-(const PALMMatrix &rhs) const {
	assert ((xSize == rhs.getXSize()) && (ySize == rhs.getYSize()));
	
	PALMMatrix <T> result(*this);	// make a copy of the current matrix
	result -= rhs;
	
	return result;
}

template <typename T> PALMMatrix<T> PALMMatrix<T>::operator*(const PALMMatrix &rhs) const {
	assert ((xSize == rhs.getXSize()) && (ySize == rhs.getYSize()));
	
	PALMMatrix <T> result(*this);	// make a copy of the current matrix
	result *= rhs;
	
	return result;
}

template <typename T> PALMMatrix<T> PALMMatrix<T>::operator/(const PALMMatrix &rhs) const {
	assert ((xSize == rhs.getXSize()) && (ySize == rhs.getYSize()));
	
	PALMMatrix <T> result(*this);	// make a copy of the current matrix
	result /= rhs;
	
	return result;
}

template <typename T> const PALMMatrix<T> PALMMatrix<T>::AddScalar(const double scalar) {
	PALMMatrix <T> result(*this);	// make a copy of the current matrix
	
	int i, j;
	#pragma omp parallel for private(i, j) shared(result) num_threads(2)
	for (i = 0; i < xSize; ++i) {
		for (j = 0; j < ySize; ++j) {
			result(i, j) += scalar;
		}
	}
	return result;
}

template <typename T> const PALMMatrix<T> PALMMatrix<T>::SubtractScalar(const double scalar) {
	PALMMatrix <T> result(*this);	// make a copy of the current matrix
	
	int i, j;
	#pragma omp parallel for private(i, j) shared(result) num_threads(2)
	for (i = 0; i < xSize; ++i) {
		for (j = 0; j < ySize; ++j) {
			result(i, j) -= scalar;
		}
	}
	return result;
}

template <typename T> const PALMMatrix<T> PALMMatrix<T>::MultiplyWithScalar(const double scalar) {
	PALMMatrix <T> result(*this);	// make a copy of the current matrix
	
	int i, j;
	#pragma omp parallel for private(i, j) shared(result) num_threads(2)
	for (i = 0; i < xSize; ++i) {
		for (j = 0; j < ySize; ++j) {
			result(i, j) *= scalar;
		}
	}
	return result;
}

template <typename T> const PALMMatrix<T> PALMMatrix<T>::DivideByScalar(const double scalar) {
	PALMMatrix <T> result(*this);	// make a copy of the current matrix
	
	int i, j;
	#pragma omp parallel for private(i, j) shared(result) num_threads(2)
	for (i = 0; i < xSize; ++i) {
		for (j = 0; j < ySize; ++j) {
			result(i, j) /= scalar;
		}
	}
	return result;
}

template <typename T> const PALMMatrix<T> PALMMatrix<T>::RaiseToPower(const double power) {
	PALMMatrix <T> result(*this);	// make a copy of the current matrix
	int i, j;
	#pragma omp parallel for private(i, j) shared(result) num_threads(2)
	for (i = 0; i < xSize; ++i) {
		for (j = 0; j < ySize; ++j) {
			result(i, j) = pow(result(i, j), power);
		}
	}
	return result;
}


template <typename T> class PALMVolume {
public:
	PALMVolume(size_t xSize_rhs, size_t ySize_rhs, size_t zSize_rhs);
	PALMVolume(const PALMVolume &rhs);
	~PALMVolume() {;}
	
	T & operator() (size_t x, size_t y, size_t z);
	T operator() (size_t x, size_t y, size_t z) const;
	T get(size_t x, size_t y, size_t z) const;
	void set(size_t x, size_t y, size_t z, const T &value);
	void set_all(const T &value);
	
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
	xSize = xSize_rhs;
	ySize = ySize_rhs;
	zSize = zSize_rhs;
	data.resize(zSize, PALMMatrix<T>(xSize, ySize));
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

template <typename T> inline T PALMVolume<T>::get(size_t x, size_t y, size_t z) const {
	assert ((x < xSize) && (y < ySize) && (z < zSize));
	return (data[z])(x, y);
}

template <typename T> inline void PALMVolume<T>::set(size_t x, size_t y, size_t z, const T &value) {
	assert ((x < xSize) && (y < ySize) && (z < zSize));
	data[z](x, y) = value;
}

template <typename T> inline void PALMVolume<T>::set_all(const T &value) {
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
