/*
 *  PALM_analysis_errors.h
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 30/03/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef PALM_ANALYSIS_ERRORS_H
#define PALM_ANALYSIS_ERRORS_H

#include <string>

using namespace std;

class OUT_OF_MEMORY {
public:
	OUT_OF_MEMORY(string rhs) {error = rhs;}
	~OUT_OF_MEMORY() {;}
	
	string get_error() {return error;}
	
protected:
	string error;
};
class CANNOT_OPEN_FILE {};
class DIMENSIONS_SHOULD_BE_EQUAL {};
class KERNEL_SIZE_SHOULD_BE_ODD {};
class IMAGE_INDEX_BEYOND_N_IMAGES {};
class GET_NTH_IMAGE_FILE_NOT_OPEN {};
class CANNOT_DETERMINE_SPE_STORAGE_TYPE {};
class CANNOT_OPEN_OUTPUT_FILE {};
class SIZE_OF_CHAR_IS_NOT_ONE_BYTE {};
class SIZE_OF_FLOAT_IS_NOT_FOUR_BYTES {};
class OUTPUT_FILE_ALREADY_EXISTS {};
class NUMBER_OF_AVERAGING_FRAMES_SHOULD_BE_ODD {};
class NUMBER_OF_AVERAGING_LARGER_THAN_N_FRAMES_IN_FILE {};

class ERROR_READING_FILE_DATA {
public:
	ERROR_READING_FILE_DATA(string rhs) {error = rhs;}
	~ERROR_READING_FILE_DATA() {;}
	
	string get_error() {return error;}
	
protected:
	string error;
};

class ERROR_WRITING_FILE_DATA {
public:
	ERROR_WRITING_FILE_DATA(string rhs) {error = rhs;}
	~ERROR_WRITING_FILE_DATA() {;}
	
	string get_error() {return error;}
	
protected:
	string error;
};

#endif
