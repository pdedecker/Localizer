/*
 *  PALM_analysis_PALMImages.h
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 01/10/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef PALM_ANALYSIS_PALMIMAGES_H
#define PALM_ANALYSIS_PALMIMAGES_H

#include "boost/smart_ptr.hpp"
#include "boost/thread.hpp"
#include "PALM_analysis_storage.h"


class PALMBitmapImageDeviationCalculator {
public:
	PALMBitmapImageDeviationCalculator() {;}
	virtual ~PALMBitmapImageDeviationCalculator() {;}
	
	virtual double getDeviation(boost::shared_ptr<PALMMatrix<double> > positions, size_t index) = 0;
};

class PALMBitmapImageDeviationCalculator_FitUncertainty : public PALMBitmapImageDeviationCalculator {
public:
	PALMBitmapImageDeviationCalculator_FitUncertainty(double scaleFactor_rhs, double upperLimit_rhs) {scaleFactor = scaleFactor_rhs; upperLimit = upperLimit_rhs;}
	~PALMBitmapImageDeviationCalculator_FitUncertainty() {;}
	
	double getDeviation(boost::shared_ptr<PALMMatrix<double> > positions, size_t index) {return ((positions->get(index, 7) < upperLimit) ? (positions->get(index, 7) * scaleFactor) : -1);}
	
private:
	double scaleFactor;
	double upperLimit;
};

class PALMBitmapImageDeviationCalculator_Constant : public PALMBitmapImageDeviationCalculator {
public:
	PALMBitmapImageDeviationCalculator_Constant(double deviation_rhs) {deviation = deviation_rhs;}
	PALMBitmapImageDeviationCalculator_Constant() {;}
	
	double getDeviation(boost::shared_ptr<PALMMatrix<double> > positions, size_t index) {return deviation;}
	
private:
	double deviation;
};

class PALMBitmapImageDeviationCalculator_AmplitudeSquareRoot : public PALMBitmapImageDeviationCalculator {
public:
	PALMBitmapImageDeviationCalculator_AmplitudeSquareRoot(double PSFWidth_rhs, double scaleFactor_rhs) {PSFWidth = PSFWidth_rhs; scaleFactor = scaleFactor_rhs;}
	PALMBitmapImageDeviationCalculator_AmplitudeSquareRoot() {;}
	
	double getDeviation(boost::shared_ptr<PALMMatrix<double> > positions, size_t index) {return PSFWidth / (scaleFactor * sqrt(positions->get(index, 1)));}
	
private:
	double PSFWidth;
	double scaleFactor;
};

class calculate_PALM_bitmap_image_ThreadStartParameters {
public:
	calculate_PALM_bitmap_image_ThreadStartParameters() {;}
	~calculate_PALM_bitmap_image_ThreadStartParameters() {;}
	
	boost::shared_ptr<PALMMatrix<double> > positions;
	boost::shared_ptr<PALMVolume <unsigned short> > image;
	boost::shared_ptr<PALMMatrix<double> > totalIntensities;
	
	boost::shared_ptr<PALMMatrix<double> > colors;
	boost::shared_ptr<PALMBitmapImageDeviationCalculator> deviationCalculator;
	
	int normalizeColors;
	size_t nFrames;
	size_t startIndex;
	size_t endIndex;
	size_t imageWidth;
	size_t imageHeight;
	size_t xSize;
	size_t ySize;
	double maxAmplitude;	// maximum amplitude of a fitted Gaussian over the entire fitted positions
	double scaleFactor;
};

boost::shared_ptr<PALMVolume <unsigned short> > calculate_PALM_bitmap_image(boost::shared_ptr<PALMMatrix<double> > positions, boost::shared_ptr<PALMMatrix<double> > colors, boost::shared_ptr<PALMBitmapImageDeviationCalculator> deviationCalculator,
																			size_t xSize, size_t ySize, size_t imageWidth, size_t imageHeight, int normalizeColors);

boost::shared_ptr<PALMVolume <unsigned short> > calculate_PALM_bitmap_image_parallel(boost::shared_ptr<PALMMatrix<double> > positions, boost::shared_ptr<PALMMatrix<double> > colors, boost::shared_ptr<PALMBitmapImageDeviationCalculator> deviationCalculator,
																					 size_t xSize, size_t ySize, size_t imageWidth, size_t imageHeight, int normalizeColors);

void calculate_PALM_bitmap_image_ThreadStart(boost::shared_ptr<calculate_PALM_bitmap_image_ThreadStartParameters> startParameters);

#endif
