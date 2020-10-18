#include "TIFFCompression.h"

#include "zlib.h"

size_t Inflate_zlib(std::vector<std::uint8_t>& compressed, std::uint8_t* inflateBuffer, size_t maxNBytes) {
	z_stream strm;
	strm.zalloc = Z_NULL;
	strm.zfree = Z_NULL;
	strm.opaque = Z_NULL;
	strm.avail_in = 0;
	strm.next_in = Z_NULL;
	int ret = inflateInit(&strm);
	if (ret != Z_OK) {
		throw std::runtime_error("error initing zlib");
	}

	strm.avail_in = compressed.size();
	strm.next_in = compressed.data();
	strm.avail_out = maxNBytes;
	strm.next_out = inflateBuffer;
	ret = inflate(&strm, Z_FINISH);
	if (ret != Z_STREAM_END) {
		throw std::runtime_error("zlib did not return Z_STREAM_END");
	}
	return strm.total_out;
}
