/*
 Copyright 2008-2011 Peter Dedecker.
 
 This file is part of Localizer.
 
 Localizer is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 Localizer is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with Localizer.  If not, see <http://www.gnu.org/licenses/>.
 
 
 Additional permission under GNU GPL version 3 section 7
 
 If you modify this Program, or any covered work, by 
 linking or combining it with libraries required for interaction 
 with analysis programs such as Igor Pro or Matlab, or to acquire
 data from or control hardware related to an experimental measurement,
 the licensors of this Program grant you additional permission 
 to convey the resulting work.
 */

#ifndef PALM_ANALYSIS_SOFI_H
#define PALM_ANALYSIS_SOFI_H

#include <queue>
#include <vector>

#include <eigen3/Eigen/Eigen>
#include "boost/smart_ptr.hpp"
#include <gsl/gsl_min.h>

#include "PALM_analysis_FileIO.h"
#include "PALM_analysis_ProgressReporting.h"

class SOFIFrameVerifier;

// exception classes
class SOFINoImageInCalculation;

void DoSOFIAnalysis(boost::shared_ptr<ImageLoader> imageLoader, boost::shared_ptr<ImageOutputWriter> outputWriter,
					std::vector<boost::shared_ptr<SOFIFrameVerifier> > frameVerifiers, boost::shared_ptr<ProgressReporter> progressReporter,
					size_t nFramesToSkip, size_t nFramesToInclude, int lagTime, int order, int crossCorrelate, int nFramesToGroup);

/* The precise SOFI calculation depends on the type of calculation (order, crosscorrelation or not, etc)
 * In addition it's possible that the calculation will only operate on subranges
 * determined by a max number of allowable frames or a bleaching limit.
 * to avoid lots of boilerplate code implement each calculation in a separate class,
 * where each class accepts a new image passed in, and returns the accumulated result
 * via getResult(). Adding a new image after a getResult() call leads to the start of a new
 * calculation.
 */
class SOFICalculator {
public:
	SOFICalculator() {;}
	virtual ~SOFICalculator() {;}
	
	virtual void addNewImage(ImagePtr image) = 0;
	virtual ImagePtr getResult() = 0;
	
protected:
};

class SOFICalculator_Order2_auto : public SOFICalculator {
public:
	SOFICalculator_Order2_auto(int lagTime);
	~SOFICalculator_Order2_auto() {;}
	
	void addNewImage(ImagePtr image);
	ImagePtr getResult();
	
protected:
	size_t lagTime;
	std::queue<ImagePtr> imageQueue;
	ImagePtr outputImage;
	ImagePtr averageImage;
	size_t nEvaluations;
};

class SOFICalculator_Order2_cross : public SOFICalculator {
public:
	SOFICalculator_Order2_cross(int lagTime);
	~SOFICalculator_Order2_cross() {;}
	
	void addNewImage(ImagePtr image);
	ImagePtr getResult();
	
protected:
	size_t lagTime;
	std::queue<ImagePtr> imageQueue;
	
	ImagePtr outputImageCrossCorrelation;
	ImagePtr outputImageHorizontalAutoCorrelation;
	ImagePtr outputImageVerticalAutoCorrelation;
	ImagePtr averageImage;
	
	size_t nEvaluations;
};

class SOFICorrector_Order2 {
public:
	SOFICorrector_Order2() {;}
	~SOFICorrector_Order2() {;}
	
	ImagePtr doImageCorrection(ImagePtr imageToCorrect);
	
protected:
	//double determinePSFStdDev(ImagePtr imageToCorrect);
	//static ImagePtr performPSFCorrection(Image *image, double psfStdDev);
	//static double functionToMinimize(double psfStdDev, void *params);
};

class SOFIFrameVerifier {
public:
	SOFIFrameVerifier() {;}
	virtual ~SOFIFrameVerifier() {;}
	
	virtual int isValidFrame(ImagePtr frame) = 0;
};

class SOFIFrameVerifier_NoSaturation : public SOFIFrameVerifier {
public:
	SOFIFrameVerifier_NoSaturation(int storageType_rhs);
	~SOFIFrameVerifier_NoSaturation() {;}
	
	int isValidFrame(ImagePtr frame);
	
protected:
	double saturationValue;
};

class SOFIFrameVerifier_MaxPixelValue : public SOFIFrameVerifier {
public:
	SOFIFrameVerifier_MaxPixelValue(int maxPixelValue_rhs) {maxPixelValue = maxPixelValue_rhs;}
	~SOFIFrameVerifier_MaxPixelValue() {;}
	
	int isValidFrame(ImagePtr frame);
	
protected:
	double maxPixelValue;
};

class SOFINoImageInCalculation : public std::runtime_error {
public:
	SOFINoImageInCalculation(const std::string& error_message) :
	std::runtime_error(error_message) {}
};

#endif
