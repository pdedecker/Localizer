#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <fstream>

#include "TIFFDefs.h"

template <typename T>
T ReadAndSwapIfNeeded(std::ifstream& f, T& val, bool needsSwap) {
	f.read(reinterpret_cast<char*>(&val), sizeof(val));
	val = MaybeEndianSwap(val, needsSwap);
	return val;
}

template <typename T>
T SwapBinTagIfNeeded(T& tag, bool needsSwap) {
	tag.tagCode = MaybeEndianSwap(tag.tagCode, needsSwap);
	tag.dataType = MaybeEndianSwap(tag.dataType, needsSwap);
	tag.nValues = MaybeEndianSwap(tag.nValues, needsSwap);
	tag.tagData = MaybeEndianSwap(tag.tagData, needsSwap);
	return tag;
}

template <typename T>
std::vector<std::uint64_t> ExtractNumericValues(const std::vector<std::uint8_t>& rawValues, bool shouldEndianSwap) {
	const T* vPtr = reinterpret_cast<const T*>(rawValues.data());
	size_t nValues = rawValues.size() / sizeof(T);
	std::vector<std::uint64_t> values(nValues);
	for (size_t i = 0; i < nValues; i += 1) {
		values[i] = vPtr[i];
		if (shouldEndianSwap) {
			values[i] = EndianSwap(values[i]);
		}
	}
	return values;
}

template <typename T>
std::vector<std::uint64_t> ExtractRationalValues(const std::vector<std::uint8_t>& rawValues, bool shouldEndianSwap) {
	const T* vPtr = reinterpret_cast<const T*>(rawValues.data());
	size_t nValues = rawValues.size() / (2 * sizeof(T));
	std::vector<std::uint64_t> values(nValues);
	for (size_t i = 0; i < nValues; i += 1) {
		values[i] = (double)MaybeEndianSwap(vPtr[2 * i], shouldEndianSwap) / (double)MaybeEndianSwap(vPtr[2 * i + 1], shouldEndianSwap);
	}
	return values;
}

size_t NBytesPerValue(TagType tagType);

std::int64_t GetLastModificationTime(const std::string& path);

// https://stackoverflow.com/questions/105252/how-do-i-convert-between-big-endian-and-little-endian-values-in-c
template <typename T>
T EndianSwap(T u) {
	static_assert (CHAR_BIT == 8, "CHAR_BIT != 8");

	union {
		T u;
		unsigned char u8[sizeof(T)];
	} source, dest;

	source.u = u;

	for (size_t k = 0; k < sizeof(T); k++)
		dest.u8[k] = source.u8[sizeof(T) - k - 1];

	return dest.u;
}

template <typename T>
T MaybeEndianSwap(T u, bool doIt) {
	if (doIt) {
		return EndianSwap(u);
	} else {
		return u;
	}
}

#endif
