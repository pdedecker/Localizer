#ifndef SOFIFRAMEVERIFIERS_H
#define SOFIFRAMEVERIFIERS_H

#include "PALM_analysis_defines.h"

class SOFIFrameVerifier {
public:
	SOFIFrameVerifier() {;}
	virtual ~SOFIFrameVerifier() {;}
	
	virtual bool isValidFrame(ImagePtr frame) const = 0;
};

class SOFIFrameVerifier_NoSaturation : public SOFIFrameVerifier {
public:
	SOFIFrameVerifier_NoSaturation(int storageType_rhs);
	~SOFIFrameVerifier_NoSaturation() {;}
	
	bool isValidFrame(ImagePtr frame) const;
	
private:
	bool validLimitsSet;
	double saturationValue;
};

class SOFIFrameVerifier_MaxPixelValue : public SOFIFrameVerifier {
public:
	SOFIFrameVerifier_MaxPixelValue(int maxPixelValue_rhs) {maxPixelValue = maxPixelValue_rhs;}
	~SOFIFrameVerifier_MaxPixelValue() {;}
	
	bool isValidFrame(ImagePtr frame) const;
	
private:
	double maxPixelValue;
};

#endif
