#include "TIFFTag.h"

#include <algorithm>
#include <sstream>
#include <vector>

#include "TIFFDefs.h"
#include "TIFFUtils.h"

#pragma pack(push, 1)
struct BinTag {
	std::uint16_t tagCode;
	std::uint16_t dataType;
	std::uint64_t nValues;
	std::uint64_t tagData;
};

struct LittleBinTag {
	std::uint16_t tagCode;
	std::uint16_t dataType;
	std::uint32_t nValues;
	std::uint32_t tagData;
};
#pragma pack(pop)

BinTag LittleBinTagToBinTag(const LittleBinTag& littleTag) {
	BinTag tag;
	tag.tagCode = littleTag.tagCode;
	tag.dataType = littleTag.dataType;
	tag.nValues = littleTag.nValues;
	tag.tagData = littleTag.tagData;
	return tag;
}

template <typename T>
T SwapBinTagMetaDataIfNeeded(T& tag, bool needsSwap) {
	// don't swap the tagData field.
	tag.tagCode = MaybeEndianSwap(tag.tagCode, needsSwap);
	tag.dataType = MaybeEndianSwap(tag.dataType, needsSwap);
	tag.nValues = MaybeEndianSwap(tag.nValues, needsSwap);
	return tag;
}

TIFFTag::TIFFTag(std::ifstream & f, bool isBigTiff, bool shouldEndianSwap) :
	_isBigTiff(isBigTiff),
	_shouldEndianSwap(shouldEndianSwap) {
	BinTag binTag;
	LittleBinTag littleBinTag;
	if (isBigTiff) {
		f.read(reinterpret_cast<char*>(&binTag), sizeof(binTag));
		SwapBinTagMetaDataIfNeeded(binTag, _shouldEndianSwap);
	} else {
		f.read(reinterpret_cast<char*>(&littleBinTag), sizeof(littleBinTag));
		SwapBinTagMetaDataIfNeeded(littleBinTag, _shouldEndianSwap);
		binTag = LittleBinTagToBinTag(littleBinTag);
	}
	std::streampos nextTagPos = f.tellg();

	size_t nBytesPerValue = NBytesPerValue((TagType)binTag.dataType);
	size_t nBytesRequired = nBytesPerValue * binTag.nValues;
	_tagData.resize(nBytesRequired);
	if ((_isBigTiff && (nBytesRequired <= 8)) || (!_isBigTiff && (nBytesRequired <= 4))) {
		if (_isBigTiff) {
			memcpy(_tagData.data(), &(binTag.tagData), nBytesRequired);
		} else {
			memcpy(_tagData.data(), &(littleBinTag.tagData), nBytesRequired);
		}
	} else {
		std::uint64_t dataOffset = 0;
		if (_isBigTiff) {
			dataOffset = MaybeEndianSwap(binTag.tagData, shouldEndianSwap);
		} else {
			dataOffset = MaybeEndianSwap(littleBinTag.tagData, shouldEndianSwap);
		}
		f.seekg(dataOffset);
		f.read(reinterpret_cast<char*>(_tagData.data()), nBytesRequired);
	}

	if (_shouldEndianSwap && (nBytesPerValue > 1)) {
		EndianSwapTIFFArray(_tagData.data(), (TagType)binTag.dataType, _tagData.size());
	}

	f.seekg(nextTagPos);
	_tagCode = binTag.tagCode;
	_tagType = (TagType)binTag.dataType;
	_nValues = binTag.nValues;
}

std::vector<std::uint64_t> TIFFTag::getNumericValues() const {
	switch (_tagType) {
	case TIFF_BYTE:
		return ExtractNumericValues<std::uint8_t>(_tagData);
		break;
	case TIFF_ASCII:
		throw std::runtime_error("Numeric access on ascii tag");
		break;
	case TIFF_SHORT:
		return ExtractNumericValues<std::uint16_t>(_tagData);
		break;
	case TIFF_LONG:
		return ExtractNumericValues<std::uint32_t>(_tagData);
		break;
	case TIFF_RATIONAL:
		return ExtractRationalValues<std::uint32_t>(_tagData);
		break;
	case TIFF_SBYTE:
		return ExtractNumericValues<std::int8_t>(_tagData);
		break;
	case TIFF_UNDEFINED:
		return ExtractNumericValues<std::uint8_t>(_tagData);
		break;
	case TIFF_SSHORT:
		return ExtractNumericValues<std::int16_t>(_tagData);
		break;
	case TIFF_SLONG:
		return ExtractNumericValues<std::int32_t>(_tagData);
		break;
	case TIFF_SRATIONAL:
		return ExtractRationalValues<std::int32_t>(_tagData);
		break;
	case TIFF_FLOAT:
		return ExtractNumericValues<float>(_tagData);
		break;
	case TIFF_DOUBLE:
		return ExtractNumericValues<double>(_tagData);
		break;
	case TIFF_LONG8:
		return ExtractNumericValues<std::uint64_t>(_tagData);
		break;
	case TIFF_SLONG8:
		return ExtractNumericValues<std::int64_t>(_tagData);
		break;
	case TIFF_IFD8:
		return ExtractNumericValues<std::uint64_t>(_tagData);
		break;
	default:
		throw std::runtime_error("unknown TIFF data type while reading numeric tag");
		break;
	}
	return std::vector<std::uint64_t>();
}

std::vector<std::string> TIFFTag::getStringValues() const {
	std::vector<std::string> strings;

	int start = 0;
	for (int i = 0; i < _tagData.size(); i += 1) {
		if (_tagData[i] == 0) {
			if ((i - start) > 0) {
				strings.emplace_back(reinterpret_cast<const char*>(_tagData.data() + start), i - start);
			}
			start = i + 1;
		}
	}

	return strings;
}

bool TIFFTag::isNumeric() const {
	return _tagType != TIFF_ASCII;
}

std::map<int, const char*> TIFFTag::_gTIFFTagsMap = std::map<int, const char*>();

std::string TIFFTag::descriptor() const {
	std::stringstream ss;
	ss << "tagcode: " << _tagCode << "(" << tagSymbolicName() << ") ";
	if (!isNumeric()) {
		ss << "string:\n";
		for (const auto& s : getStringValues()) {
			ss << "\t" << s << "\n";
		}
	} else {
		ss << "numeric:\n";
		std::vector<std::uint64_t> values = getNumericValues();
		for (size_t i = 0; i < std::min(values.size(), (size_t)10); i += 1) {
			ss << values.at(i) << ",";
		}
		if (values.size() > 10) {
			ss << " and " << values.size() - 10 << " other values";
		}
		ss << "\n";
	}
	return ss.str();
}

const char * TIFFTag::_getTagDescriptor(int tagCode) const {
	if (_gTIFFTagsMap.empty()) {
		std::map<int, const char*>&m = _gTIFFTagsMap;
		m.insert({ 254,"TIFFTAG_SUBFILETYPE" });
		m.insert({ 255,"TIFFTAG_OSUBFILETYPE" });
		m.insert({ 256,"TIFFTAG_IMAGEWIDTH" });
		m.insert({ 257,"TIFFTAG_IMAGELENGTH" });
		m.insert({ 258,"TIFFTAG_BITSPERSAMPLE" });
		m.insert({ 259,"TIFFTAG_COMPRESSION" });
		m.insert({ 262,"TIFFTAG_PHOTOMETRIC" });
		m.insert({ 263,"TIFFTAG_THRESHHOLDING" });
		m.insert({ 264,"TIFFTAG_CELLWIDTH" });
		m.insert({ 265,"TIFFTAG_CELLLENGTH" });
		m.insert({ 266,"TIFFTAG_FILLORDER" });
		m.insert({ 269,"TIFFTAG_DOCUMENTNAME" });
		m.insert({ 270,"TIFFTAG_IMAGEDESCRIPTION" });
		m.insert({ 271,"TIFFTAG_MAKE" });
		m.insert({ 272,"TIFFTAG_MODEL" });
		m.insert({ 273,"TIFFTAG_STRIPOFFSETS" });
		m.insert({ 274,"TIFFTAG_ORIENTATION" });
		m.insert({ 277,"TIFFTAG_SAMPLESPERPIXEL" });
		m.insert({ 278,"TIFFTAG_ROWSPERSTRIP" });
		m.insert({ 279,"TIFFTAG_STRIPBYTECOUNTS" });
		m.insert({ 280,"TIFFTAG_MINSAMPLEVALUE" });
		m.insert({ 281,"TIFFTAG_MAXSAMPLEVALUE" });
		m.insert({ 282,"TIFFTAG_XRESOLUTION" });
		m.insert({ 283,"TIFFTAG_YRESOLUTION" });
		m.insert({ 284,"TIFFTAG_PLANARCONFIG" });
		m.insert({ 285,"TIFFTAG_PAGENAME" });
		m.insert({ 286,"TIFFTAG_XPOSITION" });
		m.insert({ 287,"TIFFTAG_YPOSITION" });
		m.insert({ 288,"TIFFTAG_FREEOFFSETS" });
		m.insert({ 289,"TIFFTAG_FREEBYTECOUNTS" });
		m.insert({ 290,"TIFFTAG_GRAYRESPONSEUNIT" });
		m.insert({ 291,"TIFFTAG_GRAYRESPONSECURVE" });
		m.insert({ 292,"TIFFTAG_GROUP3OPTIONS" });
		m.insert({ 292,"TIFFTAG_T4OPTIONS" });
		m.insert({ 293,"TIFFTAG_GROUP4OPTIONS" });
		m.insert({ 293,"TIFFTAG_T6OPTIONS" });
		m.insert({ 296,"TIFFTAG_RESOLUTIONUNIT" });
		m.insert({ 297,"TIFFTAG_PAGENUMBER" });
		m.insert({ 300,"TIFFTAG_COLORRESPONSEUNIT" });
		m.insert({ 301,"TIFFTAG_TRANSFERFUNCTION" });
		m.insert({ 305,"TIFFTAG_SOFTWARE" });
		m.insert({ 306,"TIFFTAG_DATETIME" });
		m.insert({ 315,"TIFFTAG_ARTIST" });
		m.insert({ 316,"TIFFTAG_HOSTCOMPUTER" });
		m.insert({ 317,"TIFFTAG_PREDICTOR" });
		m.insert({ 318,"TIFFTAG_WHITEPOINT" });
		m.insert({ 319,"TIFFTAG_PRIMARYCHROMATICITIES" });
		m.insert({ 320,"TIFFTAG_COLORMAP" });
		m.insert({ 321,"TIFFTAG_HALFTONEHINTS" });
		m.insert({ 322,"TIFFTAG_TILEWIDTH" });
		m.insert({ 323,"TIFFTAG_TILELENGTH" });
		m.insert({ 324,"TIFFTAG_TILEOFFSETS" });
		m.insert({ 325,"TIFFTAG_TILEBYTECOUNTS" });
		m.insert({ 326,"TIFFTAG_BADFAXLINES" });
		m.insert({ 327,"TIFFTAG_CLEANFAXDATA" });
		m.insert({ 328,"TIFFTAG_CONSECUTIVEBADFAXLINES" });
		m.insert({ 330,"TIFFTAG_SUBIFD" });
		m.insert({ 332,"TIFFTAG_INKSET" });
		m.insert({ 333,"TIFFTAG_INKNAMES" });
		m.insert({ 334,"TIFFTAG_NUMBEROFINKS" });
		m.insert({ 336,"TIFFTAG_DOTRANGE" });
		m.insert({ 337,"TIFFTAG_TARGETPRINTER" });
		m.insert({ 338,"TIFFTAG_EXTRASAMPLES" });
		m.insert({ 339,"TIFFTAG_SAMPLEFORMAT" });
		m.insert({ 340,"TIFFTAG_SMINSAMPLEVALUE" });
		m.insert({ 341,"TIFFTAG_SMAXSAMPLEVALUE" });
		m.insert({ 343,"TIFFTAG_CLIPPATH" });
		m.insert({ 344,"TIFFTAG_XCLIPPATHUNITS" });
		m.insert({ 345,"TIFFTAG_YCLIPPATHUNITS" });
		m.insert({ 346,"TIFFTAG_INDEXED" });
		m.insert({ 347,"TIFFTAG_JPEGTABLES" });
		m.insert({ 351,"TIFFTAG_OPIPROXY" });
		m.insert({ 400,"TIFFTAG_GLOBALPARAMETERSIFD" });
		m.insert({ 401,"TIFFTAG_PROFILETYPE" });
		m.insert({ 402,"TIFFTAG_FAXPROFILE" });
		m.insert({ 403,"TIFFTAG_CODINGMETHODS" });
		m.insert({ 404,"TIFFTAG_VERSIONYEAR" });
		m.insert({ 405,"TIFFTAG_MODENUMBER" });
		m.insert({ 433,"TIFFTAG_DECODE" });
		m.insert({ 434,"TIFFTAG_IMAGEBASECOLOR" });
		m.insert({ 435,"TIFFTAG_T82OPTIONS" });
		m.insert({ 512,"TIFFTAG_JPEGPROC" });
		m.insert({ 513,"TIFFTAG_JPEGIFOFFSET" });
		m.insert({ 514,"TIFFTAG_JPEGIFBYTECOUNT" });
		m.insert({ 515,"TIFFTAG_JPEGRESTARTINTERVAL" });
		m.insert({ 517,"TIFFTAG_JPEGLOSSLESSPREDICTORS" });
		m.insert({ 518,"TIFFTAG_JPEGPOINTTRANSFORM" });
		m.insert({ 519,"TIFFTAG_JPEGQTABLES" });
		m.insert({ 520,"TIFFTAG_JPEGDCTABLES" });
		m.insert({ 521,"TIFFTAG_JPEGACTABLES" });
		m.insert({ 529,"TIFFTAG_YCBCRCOEFFICIENTS" });
		m.insert({ 530,"TIFFTAG_YCBCRSUBSAMPLING" });
		m.insert({ 531,"TIFFTAG_YCBCRPOSITIONING" });
		m.insert({ 532,"TIFFTAG_REFERENCEBLACKWHITE" });
		m.insert({ 559,"TIFFTAG_STRIPROWCOUNTS" });
		m.insert({ 700,"TIFFTAG_XMLPACKET" });
		m.insert({ 32781,"TIFFTAG_OPIIMAGEID" });
		m.insert({ 32953,"TIFFTAG_REFPTS" });
		m.insert({ 32954,"TIFFTAG_REGIONTACKPOINT" });
		m.insert({ 32955,"TIFFTAG_REGIONWARPCORNERS" });
		m.insert({ 32956,"TIFFTAG_REGIONAFFINE" });
		m.insert({ 32995,"TIFFTAG_MATTEING" });
		m.insert({ 32996,"TIFFTAG_DATATYPE" });
		m.insert({ 32997,"TIFFTAG_IMAGEDEPTH" });
		m.insert({ 32998,"TIFFTAG_TILEDEPTH" });
		m.insert({ 33300,"TIFFTAG_PIXAR_IMAGEFULLWIDTH" });
		m.insert({ 33301,"TIFFTAG_PIXAR_IMAGEFULLLENGTH" });
		m.insert({ 33302,"TIFFTAG_PIXAR_TEXTUREFORMAT" });
		m.insert({ 33303,"TIFFTAG_PIXAR_WRAPMODES" });
		m.insert({ 33304,"TIFFTAG_PIXAR_FOVCOT" });
		m.insert({ 33305,"TIFFTAG_PIXAR_MATRIX_WORLDTOSCREEN" });
		m.insert({ 33306,"TIFFTAG_PIXAR_MATRIX_WORLDTOCAMERA" });
		m.insert({ 33405,"TIFFTAG_WRITERSERIALNUMBER" });
		m.insert({ 33421,"TIFFTAG_CFAREPEATPATTERNDIM" });
		m.insert({ 33422,"TIFFTAG_CFAPATTERN" });
		m.insert({ 33432,"TIFFTAG_COPYRIGHT" });
		m.insert({ 33723,"TIFFTAG_RICHTIFFIPTC" });
		m.insert({ 34016,"TIFFTAG_IT8SITE" });
		m.insert({ 34017,"TIFFTAG_IT8COLORSEQUENCE" });
		m.insert({ 34018,"TIFFTAG_IT8HEADER" });
		m.insert({ 34019,"TIFFTAG_IT8RASTERPADDING" });
		m.insert({ 34020,"TIFFTAG_IT8BITSPERRUNLENGTH" });
		m.insert({ 34021,"TIFFTAG_IT8BITSPEREXTENDEDRUNLENGTH" });
		m.insert({ 34022,"TIFFTAG_IT8COLORTABLE" });
		m.insert({ 34023,"TIFFTAG_IT8IMAGECOLORINDICATOR" });
		m.insert({ 34024,"TIFFTAG_IT8BKGCOLORINDICATOR" });
		m.insert({ 34025,"TIFFTAG_IT8IMAGECOLORVALUE" });
		m.insert({ 34026,"TIFFTAG_IT8BKGCOLORVALUE" });
		m.insert({ 34027,"TIFFTAG_IT8PIXELINTENSITYRANGE" });
		m.insert({ 34028,"TIFFTAG_IT8TRANSPARENCYINDICATOR" });
		m.insert({ 34029,"TIFFTAG_IT8COLORCHARACTERIZATION" });
		m.insert({ 34030,"TIFFTAG_IT8HCUSAGE" });
		m.insert({ 34031,"TIFFTAG_IT8TRAPINDICATOR" });
		m.insert({ 34032,"TIFFTAG_IT8CMYKEQUIVALENT" });
		m.insert({ 34232,"TIFFTAG_FRAMECOUNT" });
		m.insert({ 34377,"TIFFTAG_PHOTOSHOP" });
		m.insert({ 34665,"TIFFTAG_EXIFIFD" });
		m.insert({ 34675,"TIFFTAG_ICCPROFILE" });
		m.insert({ 34732,"TIFFTAG_IMAGELAYER" });
		m.insert({ 34750,"TIFFTAG_JBIGOPTIONS" });
		m.insert({ 34853,"TIFFTAG_GPSIFD" });
		m.insert({ 34908,"TIFFTAG_FAXRECVPARAMS" });
		m.insert({ 34909,"TIFFTAG_FAXSUBADDRESS" });
		m.insert({ 34910,"TIFFTAG_FAXRECVTIME" });
		m.insert({ 34911,"TIFFTAG_FAXDCS" });
		m.insert({ 37439,"TIFFTAG_STONITS" });
		m.insert({ 34929,"TIFFTAG_FEDEX_EDR" });
		m.insert({ 40965,"TIFFTAG_INTEROPERABILITYIFD" });
		m.insert({ 50706,"TIFFTAG_DNGVERSION" });
		m.insert({ 50707,"TIFFTAG_DNGBACKWARDVERSION" });
		m.insert({ 50708,"TIFFTAG_UNIQUECAMERAMODEL" });
		m.insert({ 50709,"TIFFTAG_LOCALIZEDCAMERAMODEL" });
		m.insert({ 50710,"TIFFTAG_CFAPLANECOLOR" });
		m.insert({ 50711,"TIFFTAG_CFALAYOUT" });
		m.insert({ 50712,"TIFFTAG_LINEARIZATIONTABLE" });
		m.insert({ 50713,"TIFFTAG_BLACKLEVELREPEATDIM" });
		m.insert({ 50714,"TIFFTAG_BLACKLEVEL" });
		m.insert({ 50715,"TIFFTAG_BLACKLEVELDELTAH" });
		m.insert({ 50716,"TIFFTAG_BLACKLEVELDELTAV" });
		m.insert({ 50717,"TIFFTAG_WHITELEVEL" });
		m.insert({ 50718,"TIFFTAG_DEFAULTSCALE" });
		m.insert({ 50719,"TIFFTAG_DEFAULTCROPORIGIN" });
		m.insert({ 50720,"TIFFTAG_DEFAULTCROPSIZE" });
		m.insert({ 50721,"TIFFTAG_COLORMATRIX1" });
		m.insert({ 50722,"TIFFTAG_COLORMATRIX2" });
		m.insert({ 50723,"TIFFTAG_CAMERACALIBRATION1" });
		m.insert({ 50724,"TIFFTAG_CAMERACALIBRATION2" });
		m.insert({ 50725,"TIFFTAG_REDUCTIONMATRIX1" });
		m.insert({ 50726,"TIFFTAG_REDUCTIONMATRIX2" });
		m.insert({ 50727,"TIFFTAG_ANALOGBALANCE" });
		m.insert({ 50728,"TIFFTAG_ASSHOTNEUTRAL" });
		m.insert({ 50729,"TIFFTAG_ASSHOTWHITEXY" });
		m.insert({ 50730,"TIFFTAG_BASELINEEXPOSURE" });
		m.insert({ 50731,"TIFFTAG_BASELINENOISE" });
		m.insert({ 50732,"TIFFTAG_BASELINESHARPNESS" });
		m.insert({ 50733,"TIFFTAG_BAYERGREENSPLIT" });
		m.insert({ 50734,"TIFFTAG_LINEARRESPONSELIMIT" });
		m.insert({ 50735,"TIFFTAG_CAMERASERIALNUMBER" });
		m.insert({ 50736,"TIFFTAG_LENSINFO" });
		m.insert({ 50737,"TIFFTAG_CHROMABLURRADIUS" });
		m.insert({ 50738,"TIFFTAG_ANTIALIASSTRENGTH" });
		m.insert({ 50739,"TIFFTAG_SHADOWSCALE" });
		m.insert({ 50740,"TIFFTAG_DNGPRIVATEDATA" });
		m.insert({ 50741,"TIFFTAG_MAKERNOTESAFETY" });
		m.insert({ 50778,"TIFFTAG_CALIBRATIONILLUMINANT1" });
		m.insert({ 50779,"TIFFTAG_CALIBRATIONILLUMINANT2" });
		m.insert({ 50780,"TIFFTAG_BESTQUALITYSCALE" });
		m.insert({ 50781,"TIFFTAG_RAWDATAUNIQUEID" });
		m.insert({ 50827,"TIFFTAG_ORIGINALRAWFILENAME" });
		m.insert({ 50828,"TIFFTAG_ORIGINALRAWFILEDATA" });
		m.insert({ 50829,"TIFFTAG_ACTIVEAREA" });
		m.insert({ 50830,"TIFFTAG_MASKEDAREAS" });
		m.insert({ 50831,"TIFFTAG_ASSHOTICCPROFILE" });
		m.insert({ 50832,"TIFFTAG_ASSHOTPREPROFILEMATRIX" });
		m.insert({ 50833,"TIFFTAG_CURRENTICCPROFILE" });
		m.insert({ 50834,"TIFFTAG_CURRENTPREPROFILEMATRIX" });
		m.insert({ 65535,"TIFFTAG_DCSHUESHIFTVALUES" });
		m.insert({ 65536,"TIFFTAG_FAXMODE" });
		m.insert({ 65537,"TIFFTAG_JPEGQUALITY" });
		m.insert({ 65538,"TIFFTAG_JPEGCOLORMODE" });
		m.insert({ 65539,"TIFFTAG_JPEGTABLESMODE" });
		m.insert({ 65540,"TIFFTAG_FAXFILLFUNC" });
		m.insert({ 65549,"TIFFTAG_PIXARLOGDATAFMT" });
		m.insert({ 65550,"TIFFTAG_DCSIMAGERTYPE" });
		m.insert({ 65551,"TIFFTAG_DCSINTERPMODE" });
		m.insert({ 65552,"TIFFTAG_DCSBALANCEARRAY" });
		m.insert({ 65553,"TIFFTAG_DCSCORRECTMATRIX" });
		m.insert({ 65554,"TIFFTAG_DCSGAMMA" });
		m.insert({ 65555,"TIFFTAG_DCSTOESHOULDERPTS" });
		m.insert({ 65556,"TIFFTAG_DCSCALIBRATIONFD" });
		m.insert({ 65557,"TIFFTAG_ZIPQUALITY" });
		m.insert({ 65558,"TIFFTAG_PIXARLOGQUALITY" });
		m.insert({ 65559,"TIFFTAG_DCSCLIPRECTANGLE" });
		m.insert({ 65560,"TIFFTAG_SGILOGDATAFMT" });
		m.insert({ 65561,"TIFFTAG_SGILOGENCODE" });
		m.insert({ 65562,"TIFFTAG_LZMAPRESET" });
		m.insert({ 65563,"TIFFTAG_PERSAMPLE" });
		m.insert({ 65534,"TIFFTAG_ZSTD_LEVEL" });
	}
	if (_gTIFFTagsMap.count(tagCode) == 1) {
		return _gTIFFTagsMap.at(tagCode);
	} else {
		return "unknown tag";
	}
}
