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

#include "PALM_analysis_defines.h"
#include "boost/smart_ptr.hpp"
#include "boost/thread.hpp"
#include "PALM_analysis_storage.h"

class PALMBitmapImageDeviationCalculator;

/**
 * Given a set of positions, calculate a bitmap image adding an individual Gaussian for every position to the image.
 * Assumes that the positions given are a copy of the positions wave as it is produced in Igor by the fitting routine.
 */
class PALMBitmapImageCalculator {
public:
	PALMBitmapImageCalculator(boost::shared_ptr<PALMBitmapImageDeviationCalculator> devationCalculator_rhs) {
		devationCalculator = devationCalculator_rhs;
	}
	~PALMBitmapImageCalculator() {;}
	
	boost::shared_ptr<PALMMatrix<float> > CalculateImage(boost::shared_ptr<PALMMatrix<double> > positions, size_t xSize, 
														 size_t ySize, size_t imageWidth, size_t imageHeight);
	
	
protected:
	boost::shared_ptr<PALMBitmapImageDeviationCalculator> devationCalculator;
};


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

#endif
