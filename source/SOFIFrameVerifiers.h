#ifndef SOFIFRAMEVERIFIERS_H
#define SOFIFRAMEVERIFIERS_H

#include "PALM_analysis_defines.h"

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
	bool validLimitsSet;
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

#endif
