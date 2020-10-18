#ifndef TIFFTAG_H
#define TIFFTAG_H

#include <fstream>
#include <map>
#include <vector>

#include "TIFFDefs.h"

class TIFFTag {
public:
	TIFFTag(std::ifstream& f, bool isBigTiff, bool shouldEndianSwap);	// leads next tag from current position in f. On return f.tellg() will be set to where the next tag would logically come.

	~TIFFTag() { ; }

	int tagCode() const { return _tagCode; }
	std::uint64_t nValues() const { return _nValues; };

	bool isNumeric() const;

	std::vector<std::uint64_t> getNumericValues() const;
	std::vector<std::string> getStringValues() const;
	const std::vector<std::uint8_t>& getTagData() const { return _tagData; }

	const char* tagSymbolicName() const { return _getTagDescriptor(_tagCode); }
	std::string descriptor() const;

private:
	const char* _getTagDescriptor(int tagCode) const;

	bool _isBigTiff;
	bool _shouldEndianSwap;

	int _tagCode;
	TagType _tagType;
	std::uint64_t _nValues;
	std::vector<std::uint8_t> _tagData;
	static std::map<int, const char*> _gTIFFTagsMap;
};

#endif