/*
 *  PALM_analysis_SOFI.h
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 01/06/2011.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <queue>
#include <vector>

#include <Eigen/Eigen>
#include "boost/smart_ptr.hpp"

#include "PALM_analysis_FileIO.h"

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

void DoSOFIAnalysis(boost::shared_ptr<ImageLoader> imageLoader, boost::shared_ptr<ImageOutputWriter> outputWriter,
					int lagTime, int order, int crossCorrelate);
