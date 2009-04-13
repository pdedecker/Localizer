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

#endif