#include "TIFFIFD.h"

#include <algorithm>
#include <sstream>

#include "TIFFCompression.h"
#include "TIFFUtils.h"

template <typename T>
void EndianSwapTag(T& tag) {
	tag.tagCode = EndianSwap(tag.tagCode);
	tag.dataType = EndianSwap(tag.dataType);
	tag.nValues = EndianSwap(tag.nValues);
	tag.tagData = EndianSwap(tag.tagData);
}

TIFFIFD::TIFFIFD(std::ifstream& f, bool isBigTiff, bool shouldEndianSwap) :
	_isBigTiff(isBigTiff),
	_shouldEndianSwap(shouldEndianSwap),
	_ifdOffset(0),
	_nextIFDOffset(0) {
	std::uint64_t nTags = 0;
	_ifdOffset = f.tellg();
	if (isBigTiff) {
		ReadAndSwapIfNeeded(f, nTags, _shouldEndianSwap);
	} else {
		std::uint16_t nTags16 = 0;
		ReadAndSwapIfNeeded(f, nTags16, _shouldEndianSwap);
		nTags = nTags16;
	}

	for (std::uint32_t tagIdx = 0; tagIdx < nTags; tagIdx += 1) {
		_tags.emplace_back(f, _isBigTiff, _shouldEndianSwap);
	}

	if (isBigTiff) {
		ReadAndSwapIfNeeded(f, _nextIFDOffset, _shouldEndianSwap);
	} else {
		std::uint32_t nextIFDOffset32 = 0;
		ReadAndSwapIfNeeded(f, nextIFDOffset32, _shouldEndianSwap);
		_nextIFDOffset = nextIFDOffset32;
	}
}

bool TIFFIFD::haveTag(int tagCode) const {
	auto it = std::find_if(_tags.cbegin(), _tags.cend(), [=](const TIFFTag&t) -> bool {
		return t.tagCode() == tagCode;
	});

	return (it != _tags.cend());
}

std::vector<std::uint64_t> TIFFIFD::getNumericTagValues(int tagCode) const {
	auto it = std::find_if(_tags.cbegin(), _tags.cend(), [=](const TIFFTag&t) -> bool {
		return t.tagCode() == tagCode;
	});
	if (it == _tags.cend()) {
		return std::vector<std::uint64_t>();
	}

	return it->getNumericValues();
}

std::uint64_t TIFFIFD::tagNumericValueOrError(int tagCode) const {
	std::vector<std::uint64_t> values = getNumericTagValues(tagCode);
	if (values.size() != 1) {
		throw std::runtime_error("tagNumericValueOrError() has error");
	}
	return values.at(0);
}

std::vector<std::string> TIFFIFD::getStringValues(int tagCode) const {
	auto it = std::find_if(_tags.cbegin(), _tags.cend(), [=](const TIFFTag&t) -> bool {
		return t.tagCode() == tagCode;
	});
	if (it == _tags.cend()) {
		throw std::runtime_error("asked for missing tag");
	}

	return it->getStringValues();
}

std::string TIFFIFD::tagStringValueOrError(int tagCode) const {
	std::vector<std::string> strings = getStringValues(tagCode);
	if (strings.size() != 1) {
		throw std::runtime_error("tagStringValueOrError() has error");
	}
	return strings.at(0);
}

std::string TIFFIFD::descriptor() const {
	std::stringstream ss;
	ss << "IFD at offset " << _ifdOffset << " with " << _tags.size() << " tags\n";
	for (const auto& t : _tags) {
		ss << t.descriptor();
	}
	return ss.str();
}

void TIFFIFD::getImageDimensions(std::uint64_t& imageLength, std::uint64_t& imageWidth,
	TagType& pixelType, std::uint64_t& nBytesInImage) const {
	imageLength = tagNumericValueOrError(TIFFTAG_IMAGELENGTH);
	imageWidth = tagNumericValueOrError(TIFFTAG_IMAGEWIDTH);
	int bitsPerSample = tagNumericValueOrError(TIFFTAG_BITSPERSAMPLE);
	if (bitsPerSample < 8) {
		throw std::runtime_error("bitsPerSample < 8");
	}
	auto photometricInterpretation = tagNumericValueOrError(TIFFTAG_PHOTOMETRIC);
	if ((photometricInterpretation != PHOTOMETRIC_MINISBLACK) && (photometricInterpretation != PHOTOMETRIC_MINISWHITE)) {
		throw std::runtime_error("not a grayscale image");
	}
	std::vector<std::uint64_t> values = getNumericTagValues(TIFFTAG_SAMPLEFORMAT);
	int sampleFormat = (values.size() > 0) ? values.at(0) : SAMPLEFORMAT_UINT;
	pixelType = _decodeTagType(bitsPerSample, sampleFormat);
	nBytesInImage = (std::uint64_t)imageLength * imageWidth * (bitsPerSample / 8);
}

void TIFFIFD::loadImageData(std::ifstream & f, std::uint8_t * data, std::uint64_t bufSizeInBytes) const {
	std::uint64_t imageLength, imageWidth, nBytesInImage;
	TagType pixelType;
	getImageDimensions(imageLength, imageWidth, pixelType, nBytesInImage);
	size_t bytesPerPixel = nBytesInImage / (imageLength * imageWidth);
	if (nBytesInImage != bufSizeInBytes) {
		throw std::runtime_error("size in bytes doesn't match");
	}

	std::uint64_t compression = COMPRESSION_NONE;
	if (haveTag(TIFFTAG_COMPRESSION)) {
		compression = tagNumericValueOrError(TIFFTAG_COMPRESSION);
	}
	if ((compression != COMPRESSION_NONE) && (compression != COMPRESSION_ADOBE_DEFLATE) && (compression != COMPRESSION_DEFLATE)) {
		throw std::runtime_error("unsupported compression method");
	}
	bool isDeflateCompressed = (compression != COMPRESSION_NONE);
	std::vector<std::uint64_t> stripOffsets = getNumericTagValues(TIFFTAG_STRIPOFFSETS);
	std::vector<std::uint64_t> stripByteCounts = getNumericTagValues(TIFFTAG_STRIPBYTECOUNTS);
	if (stripOffsets.size() != stripByteCounts.size()) {
		throw std::runtime_error("stripOffsets.size() != stripByteCounts.size()");
	}
	std::uint64_t nStrips = stripOffsets.size();
	std::uint64_t rowsPerStrip = 0;
	if (haveTag(TIFFTAG_ROWSPERSTRIP)) {
		rowsPerStrip = tagNumericValueOrError(TIFFTAG_ROWSPERSTRIP);
	} else {
		if (nStrips == 1) {
			rowsPerStrip = imageLength;
		} else {
			rowsPerStrip = (std::uint64_t)std::ceil((double)imageLength / (stripOffsets.size() - 1)) - 1;	// max possible rowsPerStrip that can result in this many strips
		}
	}

	std::vector<std::uint8_t> compressedStrip;
	if (isDeflateCompressed) {
		compressedStrip.resize(*std::max_element(stripByteCounts.cbegin(), stripByteCounts.cend()));
	}

	std::uint64_t nBytesStoredThusFar = 0;
	std::uint8_t* dataOffsetPtr = data;
	for (std::uint64_t stripIdx = 0; stripIdx < stripOffsets.size(); stripIdx += 1) {
		f.seekg(stripOffsets.at(stripIdx));
		size_t nBytesStoredFromStrip = 0;
		if (!isDeflateCompressed) {
			size_t nBytesToRead = std::min(stripByteCounts.at(stripIdx), bufSizeInBytes - nBytesStoredThusFar);
			f.read(reinterpret_cast<char*>(dataOffsetPtr), nBytesToRead);
			nBytesStoredFromStrip = f.gcount();
		} else {
			f.read(reinterpret_cast<char*>(compressedStrip.data()), stripByteCounts.at(stripIdx));
			nBytesStoredFromStrip = Inflate_zlib(compressedStrip, dataOffsetPtr, bufSizeInBytes - nBytesStoredThusFar);
		}
		dataOffsetPtr += nBytesStoredFromStrip;
		nBytesStoredThusFar += nBytesStoredFromStrip;
	}
	if (nBytesStoredThusFar != nBytesInImage) {
		throw std::runtime_error("not enough data to fill complete image");
	}

	if (_shouldEndianSwap && (bytesPerPixel > 1)) {
		EndianSwapTIFFArray(data, pixelType, nBytesInImage);
	}
}

std::uint64_t TIFFIFD::FindNextIFDOffset(std::ifstream& f, const bool isBigTIFF, const bool shouldEndianSwap) {
	if (isBigTIFF) {
		std::uint64_t nTags = 0;
		ReadAndSwapIfNeeded(f, nTags, shouldEndianSwap);
		f.seekg(nTags * 20, std::ios_base::cur);
		std::uint64_t nextOffset = 0;
		ReadAndSwapIfNeeded(f, nextOffset, shouldEndianSwap);
		return nextOffset;
	} else {
		std::uint16_t nTags16 = 0;
		ReadAndSwapIfNeeded(f, nTags16, shouldEndianSwap);
		f.seekg((int)nTags16 * 12, std::ios_base::cur);
		std::uint32_t nextOffset32 = 0;
		ReadAndSwapIfNeeded(f, nextOffset32, shouldEndianSwap);
		return nextOffset32;
	}
}

TagType TIFFIFD::_decodeTagType(int bitsPerPixel, int sampleFormat) const {
	switch (sampleFormat) {
	case SAMPLEFORMAT_UINT:
		switch (bitsPerPixel) {
		case 8:
			return TIFF_BYTE;
			break;
		case 16:
			return TIFF_SHORT;
			break;
		case 32:
			return TIFF_LONG;
			break;
		case 64:
			return TIFF_LONG8;
			break;
		default:
			throw std::runtime_error("unknown pixel sample format");
		}
		break;
	case SAMPLEFORMAT_INT:
		switch (bitsPerPixel) {
		case 8:
			return TIFF_SBYTE;
			break;
		case 16:
			return TIFF_SSHORT;
			break;
		case 32:
			return TIFF_SLONG;
			break;
		case 64:
			return TIFF_SLONG8;
			break;
		default:
			throw std::runtime_error("unknown pixel sample format");
		}
		break;
		break;

	case SAMPLEFORMAT_IEEEFP:
		switch (bitsPerPixel) {
		case 32:
			return TIFF_FLOAT;
			break;
		case 64:
			return TIFF_DOUBLE;
			break;
		default:
			throw std::runtime_error("unknown pixel sample format");
		}
		break;
		break;
	default:
		throw std::runtime_error("unknown pixel sample format");
	}
}

