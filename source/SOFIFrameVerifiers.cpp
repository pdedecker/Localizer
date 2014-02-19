#include "SOFIFrameVerifiers.h"

SOFIFrameVerifier_NoSaturation::SOFIFrameVerifier_NoSaturation(int storageType_rhs) {
	
	validLimitsSet = true;
	switch (storageType_rhs) {
            // some systems (e.g. the Nikon TIRF) never report absolutely saturated pixel
            // values, but instead report a different value close to it
            // so add in some arbitrary margin for the saturation detection
		case STORAGE_TYPE_INT8:
			this->saturationValue = 127.0 - 1.0;
			break;
		case STORAGE_TYPE_UINT8:
			this->saturationValue = 255.0 - 2.0;
			break;
		case STORAGE_TYPE_INT16:
			this->saturationValue = 32767.0 - 200.0;
			break;
		case STORAGE_TYPE_UINT16:
			this->saturationValue = 65535.0 - 200.0;
			break;
		case STORAGE_TYPE_INT32:
			this->saturationValue = 2147483647.0 - 500.0;
			break;
		case STORAGE_TYPE_UINT32:
			this->saturationValue = 4294967295.0 - 500.0;
			break;
		case STORAGE_TYPE_INT64:
			this->saturationValue = 9223372036854775807.0 - 1000.0;
			break;
		case STORAGE_TYPE_UINT64:
			this->saturationValue = 18446744073709551615.0 - 1000.0;
			break;
		default:
			validLimitsSet = false;
			break;
	}
}

int SOFIFrameVerifier_NoSaturation::isValidFrame(ImagePtr frame) {
	if (validLimitsSet == false)
		return 1;
	
	size_t nRows = frame->rows();
	size_t nCols = frame->cols();
	
	size_t nPixels = nRows * nCols;
	double localSaturationValue = this->saturationValue;
	
	double *dataPtr = frame->data();
	for (size_t i = 0; i < nPixels; i+=1) {
		if (dataPtr[i] >= localSaturationValue) {
			return 0;
		}
	}
	
	return 1;
}

int SOFIFrameVerifier_MaxPixelValue::isValidFrame(ImagePtr image) {
	size_t nRows = image->rows();
	size_t nCols = image->cols();
	
	size_t nPixels = nRows * nCols;
	double localMaxPixelValue = this->maxPixelValue;
	
	double *dataPtr = image->data();
	for (size_t i = 0; i < nPixels; i+=1) {
		if (dataPtr[i] > localMaxPixelValue) {
			return 0;
		}
	}
	
	return 1;
}
