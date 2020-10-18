#ifndef TIFFIFD_H
#define TIFFIFD_H

#include <fstream>
#include <vector>

#include "TIFFTag.h"

class TIFFIFD {
public:
	TIFFIFD(std::ifstream& f, bool isBigTiff, bool shouldEndianSwap);	// read the IFD at f.tellg()

	bool haveTag(int tagCode) const;
	std::vector<std::uint64_t> getNumericTagValues(int tagCode) const;
	std::uint64_t tagNumericValueOrError(int tagCode) const;
	std::vector<std::string> getStringValues(int tagCode) const;
	std::string tagStringValueOrError(int tagCode) const;
	std::string descriptor() const;

	void getImageDimensions(std::uint64_t& imageLength, std::uint64_t& imageWidth,
							TagType& pixelType, std::uint64_t& nBytesInImage) const;
	void loadImageData(std::ifstream& f, std::uint8_t* data, std::uint64_t bufSizeInBytes) const;

	std::uint64_t nextIFDOffset() const { return _nextIFDOffset; }

	static std::uint64_t FindNextIFDOffset(std::ifstream& f, const bool isBigTIFF, const bool shouldEndianSwap);

private:
	TagType _decodeTagType(int bitsPerPixel, int sampleFormat)const ;

	std::vector<TIFFTag> _tags;
	std::uint64_t _ifdOffset;
	std::uint64_t _nextIFDOffset;
	bool _isBigTiff;
	bool _shouldEndianSwap;
};

#endif
