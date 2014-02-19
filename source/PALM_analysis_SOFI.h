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
#include <list>
#include <vector>

#include <eigen3/Eigen/Eigen>
#include "boost/smart_ptr.hpp"
#include <gsl/gsl_min.h>

#include "PALM_analysis_FileIO.h"
#include "PALM_analysis_ProgressReporting.h"
#include "SOFIFrameVerifiers.h"

class SingleSOFICalculation;
class XCSOFIKernelProvider;
class SOFIFrameVerifier;

// exception classes
class SOFINoImageInCalculation;

void DoSOFIAnalysis(std::shared_ptr<ImageLoader> imageLoader, std::vector<std::shared_ptr<SOFIFrameVerifier> > frameVerifiers, 
					std::shared_ptr<ProgressReporter> progressReporter,
					int nFramesToSkip, int nFramesToInclude, std::vector<int> lagTimes, int order, int crossCorrelate, int nFramesToGroup, 
					std::vector<ImagePtr>& sofiOutputImages, std::vector<ImagePtr>& averageOutputImages);

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
	SOFICalculator(int order, std::vector<int> lagTimes_rhs, size_t batchSize) : _order(order), _lagTimes(lagTimes_rhs), _batchSize(batchSize), sumOfWeightsIncluded(0.0) {}
	virtual ~SOFICalculator() {;}
	
    void addNewImage(ImagePtr image);
    void getResult(ImagePtr &calculatedSOFIImage, ImagePtr &calculatedAverageImage);
	
protected:
    void addNewBlockFromStoredImages();
    virtual void performCalculation(ImagePtr &calculatedSOFIImage, ImagePtr &calculatedAverageImage) = 0;
    virtual void performCorrection(ImagePtr imageToCorrect) = 0;
    
    int _order;
    std::vector<int> _lagTimes;
    std::vector<ImagePtr> imageVector;
    ImagePtr sofiImage;
    ImagePtr averageImage;
    
    size_t _batchSize;
    double sumOfWeightsIncluded;
};

class SOFICalculator_AutoCorrelation : public SOFICalculator {
public:
	SOFICalculator_AutoCorrelation(int order, std::vector<int> lagTimes, size_t batchSize);
	~SOFICalculator_AutoCorrelation() {;}
	
protected:
    void performCalculation(ImagePtr &calculatedSOFIImage, ImagePtr &calculatedAverageImage);
    void performCorrection(ImagePtr imageToCorrect) {;}
};

class SOFICalculator_CrossCorrelation : public SOFICalculator {
public:
	SOFICalculator_CrossCorrelation(int order, std::vector<int> lagTimes, size_t batchSize);
	~SOFICalculator_CrossCorrelation() {;}
	
protected:
    void performCalculation(ImagePtr &calculatedSOFIImage, ImagePtr &calculatedAverageImage);
    void performCorrection(ImagePtr imageToCorrect);
};

class SOFIPixelCalculation {
public:
    std::vector<int> inputRowDeltas;
    std::vector<int> inputColDeltas; // relative (delta) coordinates of the INPUT pixels that are to be included
    void getOutputPixelCoordinates(int order, int inputRow, int inputCol, int &outputRow, int &outputCol) const;
    int outputRowDelta; // used by getOutputPixelCoordinates
    int outputColDelta;
};

class XCSOFIKernelProvider {
public:
    XCSOFIKernelProvider() {;}
    ~XCSOFIKernelProvider() {;}
    
    void getCoordinatesOfFirstUsableInputPixel(int order, int nRowsInput, int nColsInput, int &firstRow, int &firstCol) const;
    void getCoordinatesOfLastUsableInputPixel(int order, int nRowsInput, int nColsInput, int &lastRow, int &lastCol) const;
    void getSizeOfOutputImage(int order, int nRowsInput, int nColsInput, int &nRows, int &nCols) const;
    boost::shared_array<std::vector<SOFIPixelCalculation> > getKernel(int order, int lagTime, int &nRows, int &nCols) const;
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

class SOFINoImageInCalculation : public std::runtime_error {
public:
	SOFINoImageInCalculation(const std::string& error_message) :
	std::runtime_error(error_message) {}
};

#endif
