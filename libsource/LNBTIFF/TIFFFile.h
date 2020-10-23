#ifndef TIFFFILE_H
#define TIFFFILE_H

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

#include "TIFFIFD.h"

class TIFFFile {
public:
	class IFDOffsets {
	public:
		IFDOffsets() :
			modificationTime(0) 
		{}
		IFDOffsets(std::vector<std::uint64_t>& offsetsA, std::int64_t modificationTimeA) :
			offsets(offsetsA),
			modificationTime(modificationTimeA) 
		{}
		std::vector<std::uint64_t> offsets;
		std::int64_t modificationTime;
	};

	TIFFFile(const std::string& filePath, std::function<bool(std::uint64_t)> ifdSearchProgressFunc = std::function<bool(std::uint64_t)>());
	~TIFFFile();

	std::uint64_t nImages() const { return _ifdOffsets.size(); }

	void getImageDimensions(const std::uint64_t imageIndex, std::uint64_t& imageLength, std::uint64_t& imageWidth,
							TagType& pixelType, std::uint64_t& nBytesInImage);
	void loadImageData(const std::uint64_t imageIndex, std::uint8_t* data, std::uint64_t bufSizeInBytes);

	bool haveTag(const std::uint64_t imageIndex, const std::uint64_t tagCode);
	std::vector<std::uint64_t> getNumericTagValues(const std::uint64_t imageIndex, const int tagCode);
	std::uint64_t getTagNumericValueOrError(const std::uint64_t imageIndex, const int tagCode);
	std::vector<std::string> getStringValues(const std::uint64_t imageIndex, const int tagCode);
	std::string getTagStringValueOrError(const std::uint64_t imageIndex, const int tagCode);

	void printDescriptor();

private:
	// verifies that this is a tiff file, sets _isBigTIFF and _needsEndianSwapping,
	// and calls f.seekg() with offset to first IFD.
	void _readTIFFHeader(std::ifstream& f);
	std::vector<std::uint64_t> _findIFDOffsets(std::ifstream& f, const std::string& filePath, std::function<bool(std::uint64_t)> ifdSearchProgressFunc) const;
	std::vector<std::uint64_t> _readIFDOffsetsFromLNBTag(std::ifstream& f) const;
	TIFFIFD& _readIFDAtIndex(const std::uint64_t index);

	//std::vector<TIFFIFD> _readIFDs(std::ifstream& f);
	//void _removeInvalidIFDs();

	bool _machineIsBigEndian() const;

	std::string _filePath;
	std::ifstream _f;
	std::vector<std::uint64_t> _ifdOffsets;
	std::shared_ptr<TIFFIFD> _currentIFDPtr;
	std::uint64_t _currentIFDIndex;
	bool _isBigTIFF;
	bool _shouldEndianSwap;
	
	static std::map<std::string, TIFFFile::IFDOffsets> _gIFDOffsetsMap;
};

#endif
