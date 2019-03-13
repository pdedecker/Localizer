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

#ifndef PALM_ANALYSIS_DEFINES_H
#define PALM_ANALYSIS_DEFINES_H

#ifdef WITH_IGOR
#define PALM_ANALYSIS_XOP_ERROR 1 + FIRST_XOP_ERR
#define ROLLING_AVERAGE_NEEDS_ODD_NUMBER_OF_FRAMES 2 + FIRST_XOP_ERR
#endif // WITH_IGOR

#include <algorithm>

#include <boost/smart_ptr.hpp>
#include <eigen3/Eigen/Eigen>
#include <fftw3.h>

// the typedef for an image
typedef Eigen::MatrixXd Image;
typedef std::shared_ptr<Image> ImagePtr;

typedef double fftw_complex__;	// work around issues with fftw_complex seemingly not working with std::shared_ptr.

const int kMinSofiOrder = 1;
const int kMaxSofiOrder = 6;

template<typename T>
T square(T val) {
	return val * val;
}

template<typename T>
T Clip (T a, T b, T c) {
    return std::max(b, std::min(a, c));
}

template <typename T>
bool Within(T a, T b, T c) {
    return (Clip(a, b, c) == a);
}

template<typename T>
void MinAndMax(const T* vals, size_t nVals, T& minVal, T& maxVal) {
    minVal = vals[0];
    maxVal = vals[0];
    for (size_t i = 1; i < nVals; i++) {
        minVal = std::min(minVal, vals[i]);
        maxVal = std::max(maxVal, vals[i]);
    }
}

// some 'reasonable' upper limits for various quantities
// values over these limits will be considered as errors
const int kMaxImageDimension = 50000;
const int kMaxNFrames = 10000000;

enum LocalizerStorageType {
    kInt4 = 0,
    kUInt4,
    kInt8,
    kUInt8,
    kInt16,
    kUInt16,
    kInt32,
    kUInt32,
    kInt64,
    kUInt64,
    kFP32,
    kFP64
};

const double PI = 3.1415926535897932384626433;
const double SQRT2 = 1.4142135623730950488;
const double SQRT2PI = 2.506628274631;

const size_t LOCALIZED_POSITIONS_TYPE_SIMULATION = -1;
const size_t LOCALIZED_POSITIONS_TYPE_2DGAUSS = 0;
const size_t LOCALIZED_POSITIONS_TYPE_2DGAUSS_FIXED_WIDTH = 1;
const size_t LOCALIZED_POSITIONS_TYPE_MULTIPLICATION = 2;
const size_t LOCALIZED_POSITIONS_TYPE_CENTROID = 3;
const size_t LOCALIZED_POSITIONS_TYPE_ZEISSPALM = 4;
const size_t LOCALIZED_POSITIONS_TYPE_ELLIPSOIDAL2DGAUSS = 5;
const size_t LOCALIZED_POSITIONS_TYPE_MLEWG = 6;

const int PALMBITMAP_EMITTERWEIGHING_SAME = 0;
const int PALMBITMAP_EMITTERWEIGHING_INTEGRAL = 1;

const int PALMBITMAP_DEVIATION_SAME = 0;
const int PALMBITMAP_DEVIATION_FITUNCERTAINTY = 1;
const int PALMBITMAP_DEVIATION_GAUSSIANMASK = 2;

const int CAMERA_TYPE_WINSPEC = 0;
const int CAMERA_TYPE_ANDOR = 1;
const int CAMERA_TYPE_HAMAMATSU = 2;
const int CAMERA_TYPE_TIFF = 3;
const int CAMERA_TYPE_PDE = 4;	// a custom, very simple image format. Not currently used anywhere
const int CAMERA_TYPE_ZEISS = 5;	// Zeiss .lsm files. Not working currently
const int CAMERA_TYPE_IGOR_WAVE = 6;
const int CAMERA_TYPE_MATLAB_MATRIX = 7;
const int CAMERA_TYPE_MULTIFILE_TIFF = 8;
const int CAMERA_TYPE_RAW_POINTER = 8;

const int PREPROCESSOR_NONE = 0;
const int PREPROCESSOR_3X3MEDIAN = 1;
const int PREPROCESSOR_5X5MEDIAN = 2;
const int PREPROCESSOR_1X1GAUSSIAN = 3;
const int PREPROCESSOR_2X2GAUSSIAN = 4;
const int PREPROCESSOR_3X3MEAN = 5;
const int PREPROCESSOR_5X5MEAN = 6;

const int THRESHOLD_METHOD_GLRT = 0;
const int THRESHOLD_METHOD_ISODATA = 1;
const int THRESHOLD_METHOD_TRIANGLE = 2;
const int THRESHOLD_METHOD_DIRECT = 3;
const int THRESHOLD_METHOD_SMOOTHSIGMA = 4;

const int POSTPROCESSOR_NONE = 0;
const int POSTPROCESSOR_REMOVE_ISOLATED_PIXELS = 1;

const int PARTICLEFINDER_ADJACENT4 = 0;
const int PARTICLEFINDER_ADJACENT8 = 1;
const int PARTICLEFINDER_RADIUS = 2;

const int PARTICLEVERIFIER_NONE = 0;
const int PARTICLEVERIFIER_SYMMETRICGAUSS = 1;
const int PARTICLEVERIFIER_ELLIPSOIDALGAUSS_SYMM = 2;
const int PARTICLEVERIFIER_REMOVEOVERLAPPINGPARTICLES = 3;
const int PARTICLEVERIFIER_ELLIPSOIDALGAUSS_ASTIG = 4;
const int PARTICLEVERIFIER_FIXEDWIDTHGAUSS = 5;


const int LOCALIZATION_METHOD_2DGAUSS = 0;
const int LOCALIZATION_METHOD_2DGAUSS_FIXEDWIDTH = 1;
const int LOCALIZATION_METHOD_MULTIPLICATION = 2;
const int LOCALIZATION_METHOD_CENTROID = 3;
const int LOCALIZATION_METHOD_2DGAUSS_ELLIPSOIDAL = 5; // 5 since Igor reserves 4 for positions fitted with Zeiss software
const int LOCALIZATION_METHOD_MLEwG = 6;
const int LOCALIZATION_METHOD_2DGAUSS_ELLIPSOIDAL_ASTIGMATISM = 7;

const int IMAGE_OUTPUT_TYPE_TIFF = 0;
const int IMAGE_OUTPUT_TYPE_COMPRESSED_TIFF = 1;
const int IMAGE_OUTPUT_TYPE_IGOR = 2;
const int IMAGE_OUTPUT_TYPE_PDE = 3;
const int IMAGE_OUTPUT_TYPE_MULTIFILE_TIFF = 4;

const int PROCESSING_AVERAGESUBTRACTION = 0;
const int PROCESSING_DIFFERENCEIMAGE = 1;
const int PROCESSING_CHANGEFORMAT = 2;
const int PROCESSING_CROP = 3;
const int PROCESSING_CONVERTTOPHOTONS = 4;

const int ANALYZING_SUMMEDTRACE = 0;
const int ANALYZING_AVERAGETRACE = 1;
const int ANALYZING_AVERAGEIMAGE = 2;
const int ANALYZING_VARIANCEIMAGE = 3;

#endif
