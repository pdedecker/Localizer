/*
 Copyright 2008-2014 Peter Dedecker.
 
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
