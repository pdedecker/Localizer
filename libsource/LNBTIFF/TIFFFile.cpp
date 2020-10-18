#include "TIFFFile.h"

#include <fstream>
#include <iostream>

#include "TIFFUtils.h"

TIFFFile::TIFFFile(const std::string & filePath, std::function<bool(std::uint64_t)> ifdSearchProgressFunc) :
	_filePath(filePath),
	_currentIFDIndex(-1),
	_isBigTIFF(false),
	_shouldEndianSwap(false) {
	_f.exceptions(std::ifstream::failbit | std::ifstream::badbit | std::ifstream::eofbit);
	_f.open(_filePath.c_str(), std::ios::binary | std::ios::in);

	_readTIFFHeader(_f);
	_ifdOffsets = _findIFDOffsets(_f, filePath, ifdSearchProgressFunc);
	_readIFDAtIndex(0);
}

TIFFFile::~TIFFFile() {
	_f.close();
}

void TIFFFile::getImageDimensions(const std::uint64_t imageIndex, std::uint64_t & imageLength, std::uint64_t& imageWidth,
	TagType & pixelType, std::uint64_t & nBytesInImage) {
	TIFFIFD& ifd = _readIFDAtIndex(imageIndex);
	ifd.getImageDimensions(imageLength, imageWidth, pixelType, nBytesInImage);
}

void TIFFFile::loadImageData(const std::uint64_t imageIndex, std::uint8_t* data, std::uint64_t bufSizeInBytes) {
	TIFFIFD& ifd = _readIFDAtIndex(imageIndex);
	ifd.loadImageData(_f, data, bufSizeInBytes);
}

bool TIFFFile::haveTag(const std::uint64_t imageIndex, const std::uint64_t tagCode) {
	TIFFIFD& ifd = _readIFDAtIndex(imageIndex);
	return ifd.haveTag(tagCode);
}

std::vector<std::uint64_t> TIFFFile::getNumericTagValues(const std::uint64_t imageIndex, const int tagCode) {
	TIFFIFD& ifd = _readIFDAtIndex(imageIndex);
	return ifd.getNumericTagValues(tagCode);
}

std::uint64_t TIFFFile::getTagNumericValueOrError(const std::uint64_t imageIndex, const int tagCode) {
	TIFFIFD& ifd = _readIFDAtIndex(imageIndex);
	return ifd.tagNumericValueOrError(tagCode);
}

std::vector<std::string> TIFFFile::getStringValues(const std::uint64_t imageIndex, const int tagCode) {
	TIFFIFD& ifd = _readIFDAtIndex(imageIndex);
	return ifd.getStringValues(tagCode);
}

std::string TIFFFile::getTagStringValueOrError(const std::uint64_t imageIndex, const int tagCode) {
	TIFFIFD& ifd = _readIFDAtIndex(imageIndex);
	return ifd.tagStringValueOrError(tagCode);
}

void TIFFFile::printDescriptor() {
	for (size_t i = 0; i < _ifdOffsets.size(); i += 1) {
		TIFFIFD& ifd = _readIFDAtIndex(i);
		std::cout << ifd.descriptor();
	}
}

void TIFFFile::_readTIFFHeader(std::ifstream& f) {
	std::uint16_t byteOrder = 0;
	f.read(reinterpret_cast<char*>(&byteOrder), sizeof(byteOrder));
	if ((byteOrder != 0x4949) && (byteOrder != 0x4D4D)) {
		throw std::runtime_error("Not a TIFF header");
	}

	bool dataIsBigEndian = (byteOrder == 0x4D4D);
	if ((dataIsBigEndian && !_machineIsBigEndian()) || (!dataIsBigEndian && _machineIsBigEndian())) {
		_shouldEndianSwap = true;
	} else {
		_shouldEndianSwap = false;
	}

	std::uint16_t magic = 0;
	ReadAndSwapIfNeeded(f, magic, _shouldEndianSwap);
	if ((magic != 42) && (magic != 43)) {
		throw std::runtime_error("Not small or big tiff header");
	}

	_isBigTIFF = (magic == 43);
	if (_isBigTIFF) {
		std::uint16_t offsetSize, constant0;
		ReadAndSwapIfNeeded(f, offsetSize, _shouldEndianSwap);
		ReadAndSwapIfNeeded(f, constant0, _shouldEndianSwap);
		if ((offsetSize != 8) || (constant0 != 0)) {
			throw std::runtime_error("not valid bigtiff header");
		}
	}

	std::uint64_t ifdOffset = 0;
	std::uint32_t ifdOffset32 = 0;
	if (_isBigTIFF) {
		ReadAndSwapIfNeeded(f, ifdOffset, _shouldEndianSwap);
	} else {
		ReadAndSwapIfNeeded(f, ifdOffset32, _shouldEndianSwap);
		ifdOffset = ifdOffset32;
	}

	f.seekg(ifdOffset);
}

std::map<std::string, TIFFFile::IFDOffsets> TIFFFile::_gIFDOffsetsMap;// = std::map<std::string, TIFFFile::IFDOffsets>();

std::vector<std::uint64_t> TIFFFile::_findIFDOffsets(std::ifstream& f, const std::string& filePath, std::function<bool(std::uint64_t)> ifdSearchProgressFunc) const {
#ifdef LNBTIFF_CACHEIFDOFFSETS
	std::uint64_t lastModificationTime = GetLastModificationTime(filePath);
	
	// do we already have offsets for this file?
	if (_gIFDOffsetsMap.count(filePath) == 1) {
		// we have information, now check if the timestamp is still okay.
		if (lastModificationTime != _gIFDOffsetsMap[filePath].modificationTime) {
			// looks like the file was modified
			// the offsets information will be overwritten below.
		} else {
			// we still have valid offsets
			return _gIFDOffsetsMap[_filePath].offsets;
		}
	}
#endif

	// if we are still here then we didn't have valid offsets stored.

	// assumes f.tellg() is at the first IFD offset, and also adds this offset to the list.
	std::vector<std::uint64_t> offsets;
	std::unordered_set<std::uint64_t> ifdOffsetsSeenSet;
	offsets.push_back(f.tellg());
	ifdOffsetsSeenSet.insert(f.tellg());
	for (; ; ) {
		if (ifdSearchProgressFunc) {
			bool shouldAbort = ifdSearchProgressFunc(offsets.size());
			if (shouldAbort) {
				throw std::runtime_error("user abort while loading TIFF file");
			}
		}
		std::uint64_t nextOffset = TIFFIFD::FindNextIFDOffset(f, _isBigTIFF, _shouldEndianSwap);
		if (nextOffset == 0) {
			break;
		}
		if (ifdOffsetsSeenSet.count(nextOffset) != 0) {
			throw std::runtime_error("circular IFDs");
		}

		offsets.push_back(nextOffset);
		ifdOffsetsSeenSet.insert(nextOffset);
		f.seekg(nextOffset);
	}

#ifdef LNBTIFF_CACHEIFDOFFSETS
	_gIFDOffsetsMap.insert({ filePath, IFDOffsets(offsets, lastModificationTime) });
#endif

	return offsets;
}

TIFFIFD& TIFFFile::_readIFDAtIndex(const std::uint64_t index) {
	if (index != _currentIFDIndex) {
		std::uint64_t offset = _ifdOffsets.at(index);
		_f.seekg(offset);
		_currentIFDPtr = std::shared_ptr<TIFFIFD>(new TIFFIFD(_f, _isBigTIFF, _shouldEndianSwap));
		_currentIFDIndex = index;
	}
	return (*_currentIFDPtr);
}

/*std::vector<TIFFIFD> TIFFFile::_readIFDs(std::ifstream & f) {
for ( ; ; ) {
_IFDs.emplace_back(f, _isBigTIFF, _shouldEndianSwap);
std::uint64_t nextIFDOffset = _IFDs.at(_IFDs.size() - 1).nextIFDOffset();
if (nextIFDOffset == 0) {
break;
} else {
f.seekg(nextIFDOffset);
continue;
}
}
}*/

/*void TIFFFile::_removeInvalidIFDs() {
for (size_t idx = _IFDs.size() - 1; idx >= 0; idx -= 1) {
bool keepIt = true;
const TIFFIFD& ifd = _IFDs.at(idx);

std::uint64_t photometricInterpretation = ifd.tagNumericValueOrError(TIFFTAG_PHOTOMETRIC);
if ((photometricInterpretation != PHOTOMETRIC_MINISBLACK) && (photometricInterpretation != PHOTOMETRIC_MINISWHITE)) {
keepIt = false;
}
std::vector<std::uint64_t> subFileType = ifd.getNumericTagValues(TIFFTAG_SUBFILETYPE);
if (!subFileType.empty() && ((subFileType.at(0) == FILETYPE_REDUCEDIMAGE) || (subFileType.at(0) == FILETYPE_MASK))) {
keepIt = false;
}
std::uint64_t compression = ifd.tagNumericValueOrError(TIFFTAG_COMPRESSION);
if (compression != 1) {
keepIt = false;
}

if (!keepIt) {
_IFDs.erase(_IFDs.begin() + idx);
}
}
}*/

bool TIFFFile::_machineIsBigEndian() const {
	// https://stackoverflow.com/questions/1001307/detecting-endianness-programmatically-in-a-c-program
	union {
		uint32_t i;
		char c[4];
	} bint = { 0x01020304 };

	return bint.c[0] == 1;
}
