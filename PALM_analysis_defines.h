/*
 *  PALM_analysis_defines.h
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 07/03/2008.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef PALM_ANALYSIS_DEFINES_H
#define PALM_ANALYSIS_DEFINES_H

#define N_OUTPUT_PARAMS_PER_FITTED_POSITION 11

#define ENCAP_GSL_RANGE_CHECK_OFF
#define HAVE_INLINE	// gsl uses inline functions

#define N_SIMULTANEOUS_IMAGE_LOADS 40	// determines the extent of the caching
#define N_SIMULTANEOUS_IMAGE_WRITES 40	// determines the extent of the caching

#define GET_NTH_IMAGE_NULL_DEF 1 + FIRST_XOP_ERR
#define END_SHOULD_BE_LARGER_THAN_START_DEF 2 + FIRST_XOP_ERR
#define IMAGE_INDEX_BEYOND_N_IMAGES_DEF 3 + FIRST_XOP_ERR
#define GET_NTH_IMAGE_FILE_NOT_OPEN_DEF 4 + FIRST_XOP_ERR
#define CANNOT_DETERMINE_SPE_STORAGE_TYPE_DEF 5 + FIRST_XOP_ERR
#define CANNOT_OPEN_FILE_DEF 6 + FIRST_XOP_ERR
#define CANNOT_OPEN_OUTPUT_FILE_DEF 7 + FIRST_XOP_ERR
#define SIZE_OF_CHAR_IS_NOT_ONE_BYTE_DEF 8 + FIRST_XOP_ERR
#define SIZE_OF_FLOAT_IS_NOT_FOUR_BYTES_DEF 9 + FIRST_XOP_ERR
#define OUTPUT_FILE_ALREADY_EXISTS_DEF 10 + FIRST_XOP_ERR
#define NUMBER_OF_AVERAGING_FRAMES_SHOULD_BE_ODD_DEF 11 + FIRST_XOP_ERR
#define UNKNOWN_CCD_IMAGES_PROCESSING_METHOD 12 + FIRST_XOP_ERR
#define UNSUPPORTED_CCD_FILE_TYPE 13 + FIRST_XOP_ERR
#define UNKNOWN_CCD_IMAGES_ANALYSIS_METHOD 14 + FIRST_XOP_ERR
#define UNKNOWN_THRESHOLD_PREPROCESSING_METHOD 15 + FIRST_XOP_ERR
#define UNKNOWN_THRESHOLD_POSTPROCESSING_METHOD 16 + FIRST_XOP_ERR
#define ERROR_READING_FILE_DATA_DEF 17 + FIRST_XOP_ERR
#define ERROR_WRITING_FILE_DATA_DEF 18 + FIRST_XOP_ERR

const int STORAGE_TYPE_INT4 = 0;
const int STORAGE_TYPE_UINT4 = 1;
const int STORAGE_TYPE_INT8 = 2;
const int STORAGE_TYPE_UINT8 = 3;
const int STORAGE_TYPE_INT16 = 4;
const int STORAGE_TYPE_UINT16 = 5;
const int STORAGE_TYPE_INT32 = 6;
const int STORAGE_TYPE_UINT32 = 7;
const int STORAGE_TYPE_INT64 = 8;
const int STORAGE_TYPE_UINT64 = 9;
const int STORAGE_TYPE_FP32 = 10;
const int STORAGE_TYPE_FP64 = 11;

const double PI = 3.1415926535897932384626433;
const double SQRT2 = 1.4142135623730950488;
const double SQRT2PI = 2.506628274631;

const size_t LOCALIZED_POSITIONS_TYPE_2DGAUSS = 0;
const size_t LOCALIZED_POSITIONS_TYPE_2DGAUSS_FIXED_WIDTH = 1;
const size_t LOCALIZED_POSITIONS_TYPE_CENTROID = 2;
const size_t LOCALIZED_POSITIONS_TYPE_MULTIPLICATION = 3;

const int PALMBITMAP_EMITTERWEIGHING_SAME = 0;
const int PALMBITMAP_EMITTERWEIGHING_INTEGRAL = 1;

const int CAMERA_TYPE_WINSPEC = 0;
const int CAMERA_TYPE_ANDOR = 1;
const int CAMERA_TYPE_HAMAMATSU = 2;
const int CAMERA_TYPE_TIFF = 3;
const int CAMERA_TYPE_SIMPLE = 4;	// a custom, very simple image format. Not currently used anywhere
const int CAMERA_TYPE_ZEISS = 5;	// Zeiss .lsm files. Not working currently
const int CAMERA_TYPE_IGOR_WAVE = 6;

const int IMAGE_OUTPUT_TYPE_TIFF = 0;
const int IMAGE_OUTPUT_TYPE_COMPRESSED_TIFF = 1;
const int IMAGE_OUTPUT_TYPE_IGOR = 2;


#endif
