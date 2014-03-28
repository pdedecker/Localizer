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

#ifndef PALM_ANALYSIS_ERRORS_H
#define PALM_ANALYSIS_ERRORS_H

#include <string>
#include <stdexcept>

class CANNOT_OPEN_FILE : public std::runtime_error {
public:
	CANNOT_OPEN_FILE(const std::string& error_message) :
		std::runtime_error(error_message) {}
};

class DIMENSIONS_SHOULD_BE_EQUAL : public std::runtime_error {
public:
	DIMENSIONS_SHOULD_BE_EQUAL(const std::string& error_message) :
	std::runtime_error(error_message) {}
};

class KERNEL_SIZE_SHOULD_BE_ODD : public std::runtime_error {
public:
	KERNEL_SIZE_SHOULD_BE_ODD(const std::string& error_message) :
	std::runtime_error(error_message) {}
};

class IMAGE_INDEX_BEYOND_N_IMAGES : public std::runtime_error {
public:
	IMAGE_INDEX_BEYOND_N_IMAGES(const std::string& error_message) :
	std::runtime_error(error_message) {}
};

class GET_NTH_IMAGE_FILE_NOT_OPEN : public std::runtime_error {
public:
	GET_NTH_IMAGE_FILE_NOT_OPEN(const std::string& error_message) :
	std::runtime_error(error_message) {}
};

class CANNOT_DETERMINE_SPE_STORAGE_TYPE : public std::runtime_error {
public:
	CANNOT_DETERMINE_SPE_STORAGE_TYPE(const std::string& error_message) :
	std::runtime_error(error_message) {}
};

class CANNOT_OPEN_OUTPUT_FILE : public std::runtime_error {
public:
	CANNOT_OPEN_OUTPUT_FILE(const std::string& error_message) :
	std::runtime_error(error_message) {}
};

class SIZE_OF_CHAR_IS_NOT_ONE_BYTE {};
class SIZE_OF_FLOAT_IS_NOT_FOUR_BYTES {};

class OUTPUT_FILE_ALREADY_EXISTS : public std::runtime_error {
public:
	OUTPUT_FILE_ALREADY_EXISTS(const std::string& error_message) :
	std::runtime_error(error_message) {}
};

class NUMBER_OF_AVERAGING_FRAMES_SHOULD_BE_ODD : public std::runtime_error {
public:
	NUMBER_OF_AVERAGING_FRAMES_SHOULD_BE_ODD(const std::string& error_message) :
	std::runtime_error(error_message) {}
};

class NUMBER_OF_AVERAGING_LARGER_THAN_N_FRAMES_IN_FILE : public std::runtime_error {
public:
	NUMBER_OF_AVERAGING_LARGER_THAN_N_FRAMES_IN_FILE(const std::string& error_message) :
	std::runtime_error(error_message) {}
};

class END_SHOULD_BE_LARGER_THAN_START : public std::runtime_error {
public:
	END_SHOULD_BE_LARGER_THAN_START(const std::string& error_message) :
	std::runtime_error(error_message) {}
};

class ERROR_READING_FILE_DATA : public std::runtime_error {
public:
	ERROR_READING_FILE_DATA(const std::string& error_message) :
	std::runtime_error(error_message) {}
};

class ERROR_WRITING_FILE_DATA : public std::runtime_error {
public:
	ERROR_WRITING_FILE_DATA(const std::string& error_message) :
	std::runtime_error(error_message) {}
};

class ERROR_RUNNING_THREADED_ANALYSIS : public std::runtime_error {
public:
	ERROR_RUNNING_THREADED_ANALYSIS(const std::string& error_message) :
	std::runtime_error(error_message) {}
};

class USER_ABORTED : public std::runtime_error {
public:
	USER_ABORTED(const std::string& error_message) :
	std::runtime_error(error_message) {}
};

#endif
