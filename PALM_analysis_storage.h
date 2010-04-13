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

/**
 * @brief A simple container holding the coordinates of a single point in 3D space
 */
class Point3D {
public:
	double xPosition;
	double yPosition;
	double zPosition;
};

/**
 * @brief A more complex container containing the coordinates of a single point in 2D space
 * with its estimated intensity and background
 */
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

#endif
