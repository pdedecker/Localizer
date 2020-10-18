#ifndef TIFFCOMPRESSION_H
#define TIFFCOMPRESSION_H

#include <vector>

// return number of bytes in decompressed output.
// throws std::runtime_error if inflated is too small.
size_t Inflate_zlib(std::vector<std::uint8_t>& compressed, std::uint8_t* inflateBuffer, size_t maxNBytes);

#endif
