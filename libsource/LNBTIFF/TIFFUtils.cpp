#include "TIFFUtils.h"

size_t NBytesPerValue(TagType tagType) {
	size_t bytesPerValue = 0;
	switch (tagType) {
	case TIFF_BYTE:
	case TIFF_SBYTE:
	case TIFF_UNDEFINED:
	case TIFF_ASCII:
		bytesPerValue = 1;
		break;
	case TIFF_SHORT:
	case TIFF_SSHORT:
		bytesPerValue = 2;
		break;
	case TIFF_LONG:
	case TIFF_SLONG:
	case TIFF_FLOAT:
		bytesPerValue = 4;
		break;
	case TIFF_RATIONAL:
	case TIFF_SRATIONAL:
	case TIFF_DOUBLE:
	case TIFF_LONG8:
	case TIFF_SLONG8:
		bytesPerValue = 8;
		break;
	default:
		throw std::runtime_error("unknown numeric tag type in NumericTagDataSize()");
		break;
	}
	return bytesPerValue;
}

/*std::int64_t GetLastModificationTime(const std::string& path) {
	std::int64_t modTime = 1;

#ifdef _WIN32
	struct _stat64 buffer;
	int err = _stati64(path.c_str(), &buffer);
	if ((err != 0) && (errno != 0)) {
		throw std::runtime_error("error from _stat()");
	}

	modTime = static_cast<std::int64_t>(buffer.st_mtime);

#else	// not _WIN32
	struct stat buffer;
	int err = stat(path.c_str(), &buffer);
	if (err != 0)
		throw std::runtime_error("non-zero return from _stat()");

	modTime = static_cast<int64_t>(buffer.st_mtime);

#endif	// _WIN32

	return modTime;
}
*/