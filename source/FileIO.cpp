/*
 Copyright 2008-2014 Peter Dedecker.
 
 This file is part of Localizer.
 
 Localizer is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 Localizer is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with Localizer.  If not, see <http://www.gnu.org/licenses/>.
 
 
 Additional permission under GNU GPL version 3 section 7
 
 If you modify this Program, or any covered work, by
 linking or combining it with libraries required for interaction
 with analysis programs such as Igor Pro or Matlab, or to acquire
 data from or control hardware related to an experimental measurement,
 the licensors of this Program grant you additional permission
 to convey the resulting work.
 */

#include "MatrixRecycler.h"
#include "FileIO.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <boost/algorithm/string.hpp>
#include <limits>
#include <zlib.h>

int GetFileStorageType(const std::string &filePath) {
    size_t startOfExtension = filePath.rfind('.');
    if (startOfExtension == size_t(-1)) {
        // the filepath does not appear to contain an extension
        throw std::runtime_error("Unable to deduce the file type");
    }
    
    std::string extension = filePath.substr(startOfExtension + 1);
    if ((extension.length() < 3) || (extension.length() > 4)) {
        throw std::runtime_error("Unable to deduce the file type");
    }
    
    if (boost::algorithm::iequals(extension, "spe"))
        return CAMERA_TYPE_WINSPEC;
    if (boost::algorithm::iequals(extension, "sif"))
        return CAMERA_TYPE_ANDOR;
    if (boost::algorithm::iequals(extension, "his"))
        return CAMERA_TYPE_HAMAMATSU;
    if (boost::algorithm::iequals(extension, "pde"))
        return CAMERA_TYPE_PDE;
    if (boost::algorithm::iequals(extension, "tif"))
        return CAMERA_TYPE_TIFF;
    if (boost::algorithm::iequals(extension, "tiff"))
        return CAMERA_TYPE_TIFF;
	if (boost::algorithm::iequals(extension, "btf"))
        return CAMERA_TYPE_TIFF;
	if (boost::algorithm::iequals(extension, "tf8"))
        return CAMERA_TYPE_TIFF;
    if (boost::algorithm::iequals(extension, "ome"))
        return CAMERA_TYPE_TIFF;
    if (boost::algorithm::iequals(extension, "lsm"))
        return CAMERA_TYPE_TIFF;
    
    // if we're still here then the extension was not recognized
    throw std::runtime_error("Unable to deduce the file type");
}

void GetFilePathAndCameraType(const std::string &inputFilePath, std::string &filePath, size_t &cameraType) {
    std::string possiblyConvertedPath = inputFilePath;
    
#ifdef WITH_IGOR
    int isWavePath = 1;
    
    try {
        FetchWaveUsingFullPath(inputFilePath);
    }
    // if the wave pointed to by the path does not the exist
    // (as would be the case if it was a file path)
    // then FetchWaveUsingFullPath will throw exceptions
    catch (...) {
        isWavePath = 0;
    }
    
    if (isWavePath == 1) {
        filePath = inputFilePath;
        cameraType = CAMERA_TYPE_IGOR_WAVE;
        return;
    }
    
    // if we're still here then it's a path to a file
    // first we need to try to convert the path to the
    // appropriate format, if required
    possiblyConvertedPath = ConvertPathToNativePath(inputFilePath);
#endif
    
    filePath = possiblyConvertedPath;
    cameraType = GetFileStorageType(filePath);
    return;
}

std::shared_ptr<ImageLoader> GetImageLoader(const std::string& data_file_path, int cameraType) {
    std::shared_ptr<ImageLoader> image_loader;
    size_t estimatedCameraType;
    std::string convertedFilePath = data_file_path;
    
    // the camera type might be unknown
    if (cameraType == -1) {
        GetFilePathAndCameraType(data_file_path, convertedFilePath, estimatedCameraType);
    } else {
        estimatedCameraType = cameraType;
#ifdef WITH_IGOR
        if (estimatedCameraType != CAMERA_TYPE_IGOR_WAVE) {
            convertedFilePath = ConvertPathToNativePath(data_file_path);
        }
#endif
    }
    
    switch (estimatedCameraType) {
		case CAMERA_TYPE_WINSPEC:	// spe files
			image_loader = std::shared_ptr<ImageLoader>(new ImageLoaderSPE(convertedFilePath));
			break;
		case CAMERA_TYPE_ANDOR:
			image_loader = std::shared_ptr<ImageLoader>(new ImageLoaderAndor(convertedFilePath));
			break;
		case CAMERA_TYPE_HAMAMATSU:
			image_loader = std::shared_ptr<ImageLoader>(new ImageLoaderHamamatsu(convertedFilePath));
			break;
		case CAMERA_TYPE_TIFF:
			image_loader = std::shared_ptr<ImageLoader>(new ImageLoaderTIFF(convertedFilePath));
			break;
		case CAMERA_TYPE_PDE:
			image_loader = std::shared_ptr<ImageLoader>(new ImageLoaderPDE(convertedFilePath));
			break;
		case CAMERA_TYPE_ZEISS:	// Zeiss lsm files
			image_loader = std::shared_ptr<ImageLoader>(new ImageLoaderTIFF(convertedFilePath));
			break;
#ifdef WITH_IGOR
		case CAMERA_TYPE_IGOR_WAVE: // Matrix wave in Igor
			image_loader = std::shared_ptr<ImageLoader>(new ImageLoaderIgor(convertedFilePath));
			break;
#endif
		case CAMERA_TYPE_MULTIFILE_TIFF:
			image_loader = std::shared_ptr<ImageLoader>(new ImageLoaderMultiFileTIFF(convertedFilePath));
			break;
		default:
			throw std::runtime_error("Unsupported CCD file type (/Y flag)");
			break;
    }
    
    return image_loader;
    
}

int64_t GetLastModificationTime(const std::string& path) {
	int64_t modTime = 1;
	
#ifdef _WIN32
	struct _stat buffer;
	int err = _stat(path.c_str(), &buffer);
	if ((err != 0) && (errno != 0)) {
		// for whatever reason _stat() will return -1 but not set errno,
		// which leads me to believe that this -1 value does not indicate an error.
		// the microsoft docs say that _stat() should return non-zero for error, and should
		// set errno as well. So I'm seeing undocumented behavior.
		char errorMessage[100];
		sprintf(errorMessage, "Return value %d from _stat(), errno is %d", err, errno);
		if (errno != 0) {
			throw std::runtime_error(errorMessage);
		}
	}
	
	modTime = static_cast<int64_t>(buffer.st_mtime);
	
#else	// not _WIN32
	struct stat buffer;
	int err = stat(path.c_str(), &buffer);
	if (err != 0)
		throw std::runtime_error("non-zero return from _stat()");
	
	modTime = static_cast<int64_t>(buffer.st_mtime);
	
#endif	// _WIN32
	
	return modTime;
}

ImagePtr BufferWithFormatToImage(const std::vector<char>& imageBuffer, int nRows, int nCols, int format, int treatAsRowMajor) {
    
	ImagePtr image(GetRecycledMatrix(nRows, nCols), FreeRecycledMatrix);
    size_t nBytesInBuffer = NBytesInImage(nRows, nCols, format);
    
    assert(imageBuffer.size() == nBytesInBuffer);
    
	switch (format) {
        case STORAGE_TYPE_INT8:
            CopyBufferToImage<int8_t>(imageBuffer, image, treatAsRowMajor);
            break;
            
        case STORAGE_TYPE_UINT8:
            CopyBufferToImage<uint8_t>(imageBuffer, image, treatAsRowMajor);
            break;
            
        case STORAGE_TYPE_INT16:
            CopyBufferToImage<int16_t>(imageBuffer, image, treatAsRowMajor);
            break;
            
        case STORAGE_TYPE_UINT16:
            CopyBufferToImage<uint16_t>(imageBuffer, image, treatAsRowMajor);
            break;
            
		case STORAGE_TYPE_UINT32:
            CopyBufferToImage<uint32_t>(imageBuffer, image, treatAsRowMajor);
            break;
            
		case STORAGE_TYPE_FP32:
            CopyBufferToImage<float>(imageBuffer, image, treatAsRowMajor);
            break;
            
		case STORAGE_TYPE_FP64:
            CopyBufferToImage<double>(imageBuffer, image, treatAsRowMajor);
            break;
            
		default:
			throw std::runtime_error("unknown format in BufferWithFormatToImage()");
			break;
	}
	
	return image;
}

void ImageToBufferWithFormat(ImagePtr image, int format, std::vector<char>& imageBuffer, int treatAsRowMajor) {
	switch (format) {
		case STORAGE_TYPE_INT8:
			CopyImageToBuffer<int8_t>(image, imageBuffer, treatAsRowMajor);
			break;
		case STORAGE_TYPE_UINT8:
			CopyImageToBuffer<uint8_t>(image, imageBuffer, treatAsRowMajor);
			break;
		case STORAGE_TYPE_INT16:
			CopyImageToBuffer<int16_t>(image, imageBuffer, treatAsRowMajor);
			break;
		case STORAGE_TYPE_UINT16:
			CopyImageToBuffer<uint16_t>(image, imageBuffer, treatAsRowMajor);
			break;
		case STORAGE_TYPE_INT32:
			CopyImageToBuffer<int32_t>(image, imageBuffer, treatAsRowMajor);
			break;
		case STORAGE_TYPE_UINT32:
			CopyImageToBuffer<uint32_t>(image, imageBuffer, treatAsRowMajor);
			break;
		case STORAGE_TYPE_FP32:
			CopyImageToBuffer<float>(image, imageBuffer, treatAsRowMajor);
			break;
		case STORAGE_TYPE_FP64:
			CopyImageToBuffer<double>(image, imageBuffer, treatAsRowMajor);
			break;
        default:
            throw std::runtime_error("unknown format in ImageToBufferWithFormat()");
			break;
	}
}

size_t NBytesInImage(int nRows, int nCols, int format) {
    size_t nPixels = nRows * nCols;
    size_t nBytesInImage;
    switch (format) {
		case STORAGE_TYPE_INT8:
		case STORAGE_TYPE_UINT8:
			nBytesInImage = nPixels;
			break;
		case STORAGE_TYPE_INT16:
		case STORAGE_TYPE_UINT16:
			nBytesInImage = nPixels * 2;
			break;
		case STORAGE_TYPE_INT32:
		case STORAGE_TYPE_UINT32:
		case STORAGE_TYPE_FP32:
			nBytesInImage = nPixels * 4;
			break;
		case STORAGE_TYPE_FP64:
			nBytesInImage = nPixels * 8;
			break;
		default:
			throw std::runtime_error("Unknown format in NBytesInImage()");
			break;
	}
    
    return nBytesInImage;
}

std::vector<char> Deflate(std::vector<char>& data) {
    int err;
    
    // allocate zlib machinery
    z_stream stream;
    stream.zalloc = Z_NULL;
    stream.zfree = Z_NULL;
    stream.opaque = Z_NULL;
    err = deflateInit(&stream, Z_DEFAULT_COMPRESSION);
    if (err != Z_OK)
        throw std::runtime_error("unable to compress");
    
    // allocate output space. We take a bit of margin to be on the safe side
    size_t nBytesInOutput = std::max(static_cast<size_t>(static_cast<double>(data.size()) * 1.1), static_cast<size_t>(1024)); // min 1024 because I'm paranoid. See http://www.zlib.net/zlib_tech.html for a discussion of worst-case compression.
    std::vector<char> output(nBytesInOutput);
    
    stream.next_in = reinterpret_cast<unsigned char*>(data.data());
    stream.next_out = reinterpret_cast<unsigned char*>(output.data());
    stream.avail_in = data.size();
    stream.avail_out = output.size();
    err = deflate(&stream, Z_FINISH);
    deflateEnd(&stream);
    if (err != Z_STREAM_END) {
        throw std::runtime_error("unable to compress");
    }
    
    size_t nBytesUsed = output.size() - stream.avail_out;
    output.resize(nBytesUsed);
    return output;
}

#ifdef _WIN32
WindowsFileStream::~WindowsFileStream() {
    if (this->fileRef != NULL) {
        fclose(this->fileRef);
        this->fileRef = NULL;
    }
}

void WindowsFileStream::open(const char *path_rhs, std::ios_base::openmode mode) {
	if (this->fileRef != NULL) {
		fclose(this->fileRef);
		this->fileRef = NULL;
	}

	std::string modeStr;
	if (mode & std::ios_base::in) {
		modeStr += "r";
	}
	if (mode & std::ios_base::out) {
		modeStr += "w";
	}
	if (mode & std::ios_base::binary) {
		modeStr += "b";
	}
    int err = fopen_s(&fileRef, path_rhs, modeStr.c_str());
    path = path_rhs;
}

void WindowsFileStream::close() {
    int err;
    if (this->fileRef != NULL) {
        err = fclose(this->fileRef);
        if (err != 0) {
			std::string error;
			std::stringstream ss;
			ss << "Error " << err << " returned using close() on the image file at \"" << path << "\"";
            error = ss.str();
            throw ERROR_READING_FILE_DATA(error);
        }
        this->fileRef = NULL;
    }
}

void WindowsFileStream::get(char& c) {
	assert (this->fileRef != NULL);
	if (this->fileRef == NULL) {
		throw ERROR_READING_FILE_DATA(std::string("\"get\" was called on a NULL FILE*"));
	}
    c = fgetc(this->fileRef);
}

void WindowsFileStream::read(char *buffer, size_t nBytes) {
	assert (this->fileRef != NULL);
    int itemsRead;
    itemsRead = fread(buffer, nBytes, 1, this->fileRef);
    if (itemsRead != 1) {	// signals that not everything was read
		std::string error;
		error = "Could not read a full buffer from the file at \"";
		error += path;
		error += "\"";
        throw ERROR_READING_FILE_DATA(error);
    }
}

void WindowsFileStream::getline(char *buffer, size_t nMax) {
	assert (this->fileRef != NULL);
    fgets(buffer, nMax, this->fileRef);
    if (ferror(fileRef) != 0) {
		std::string error;
		error = "Error returned using getline() on the image file at \"";
		error += path;
		error += "\"";
        throw ERROR_READING_FILE_DATA(error);
    }
}

void WindowsFileStream::write(const char *buffer, size_t nBytes) {
    assert (this->fileRef != NULL);
    
    int err = fwrite(buffer, nBytes, 1, this->fileRef);
    if (ferror(fileRef) != 0) {
		std::string error;
		error = "Error returned using write() on the image file at \"";
		error += path;
		error += "\"";
        throw ERROR_READING_FILE_DATA(error);
    }
}

uint64_t WindowsFileStream::tellg() {
	assert (this->fileRef != NULL);
    uint64_t pos;
    
    pos = _ftelli64(this->fileRef);
    return pos;
}

uint64_t WindowsFileStream::tellp() {
    return tellg();
}

void WindowsFileStream::seekg(uint64_t pos, std::ios_base::seekdir dir) {
	seekp(pos, dir);
}

void WindowsFileStream::seekp(uint64_t pos, std::ios_base::seekdir dir) {
    assert (this->fileRef != NULL);
    
    int origin;
    switch (dir) {
        case std::ios_base::beg:
            origin = SEEK_SET;
            break;
        case std::ios_base::cur:
            origin = SEEK_CUR;
            break;
        case std::ios_base::end:
            origin = SEEK_END;
            break;
        default:
            throw std::runtime_error("unknown seek origin in WindowsFileStream::seekp");
            break;
    }
    int err = _fseeki64(this->fileRef, pos, origin);
    if (err != 0) {
		std::string error;
       	error = "Error returned using seekg() on the image file at \"";
		error += path;
		error += "\"";
        throw ERROR_READING_FILE_DATA(error);
    }
}
#endif // _WIN32


ImageLoader::ImageLoader() {
	this->nImages = 0;
	this->xSize = 0;
	this->ySize = 0;
	this->header_length = 0;
	this->nextImageToRead = 0;
}

ImageLoader::~ImageLoader() {
	if (file.is_open() == 1)
		file.close();
}

void ImageLoader::checkForReasonableValues() {
	if ((this->xSize > kMaxImageDimension) || (this->ySize > kMaxImageDimension)) {
		throw std::runtime_error("the reported frame size is unreasonably large");
	}
	
	if (this->nImages > kMaxNFrames) {
		throw std::runtime_error("the reported number of frames is unreasonably large");
	}
}

ImagePtr ImageLoader::readImage(int index) {
	this->spoolTo(index);
	return this->readNextImage(index);
}

ImagePtr ImageLoader::readNextImage() {
	int dummy;
	return this->readNextImage(dummy);
}

ImagePtr ImageLoader::readNextImageAndLoop(int &index) {
	this->spoolTo(index % this->nImages);
	return this->readNextImage(index);
}

void ImageLoader::spoolTo(int index) {
	if (index >= this->getNImages())
		throw IMAGE_INDEX_BEYOND_N_IMAGES(std::string("Requested more images than there are in the file"));
	
    this->nextImageToRead = index;
}

ImageLoaderWrapper::ImageLoaderWrapper(std::shared_ptr<ImageLoader> baseImageLoader) :
    _baseImageLoader(baseImageLoader),
    _firstImageToInclude(0),
    _minX(0),
    _minY(0),
    _haveCustomROI(false)
{
    if ((_baseImageLoader->getNImages() < 1) || (_baseImageLoader->getXSize() < 1) || (_baseImageLoader->getYSize() < 1))
        throw std::runtime_error("ImageLoaderWrapper with invalid base loader");
    
    _lastImageToInclude = _baseImageLoader->getNImages() - 1;
    _maxX = _baseImageLoader->getXSize() - 1;
    _maxY = _baseImageLoader->getYSize() - 1;
}

int ImageLoaderWrapper::getNImages() const {
    return (_lastImageToInclude - _firstImageToInclude + 1);
}

int ImageLoaderWrapper::getXSize() const {
    return (_maxX - _minX + 1);
}

int ImageLoaderWrapper::getYSize() const {
    return (_maxY - _minY + 1);
}

int ImageLoaderWrapper::getStorageType() const {
    return _baseImageLoader->getStorageType();
}

int ImageLoaderWrapper::getFileType() const {
    return _baseImageLoader->getFileType();
}

ImagePtr ImageLoaderWrapper::readNextImage(int &indexOfImageThatWasRead) {
    ImagePtr image = _baseImageLoader->readNextImage(indexOfImageThatWasRead);
    indexOfImageThatWasRead -= _firstImageToInclude;
    if (_haveCustomROI) {
        int nOutputRows = this->getXSize();
        int nOutputCols = this->getYSize();
        ImagePtr croppedImage(new Image(nOutputRows, nOutputCols));
        for (int col = 0; col < nOutputCols; ++col) {
            for (int row = 0; row < nOutputRows; ++row) {
                (*croppedImage)(row, col) = (*image)(row + _minX, col + _minY);
            }
        }
        return croppedImage;
    } else {
        return image;
    }
}

void ImageLoaderWrapper::spoolTo(int index) {
    _baseImageLoader->spoolTo(index + _firstImageToInclude);
}

void ImageLoaderWrapper::setImageRange(int firstImageToInclude, int lastImageToInclude) {
    if (lastImageToInclude < 0)
        lastImageToInclude = _baseImageLoader->getNImages() - 1;
    if ((firstImageToInclude < 0) || (lastImageToInclude >= _baseImageLoader->getNImages()) || (firstImageToInclude > lastImageToInclude))
        throw std::runtime_error("invalid setImageRange()");
    
    _firstImageToInclude = firstImageToInclude;
    _lastImageToInclude = lastImageToInclude;
    _baseImageLoader->spoolTo(firstImageToInclude);
}

void ImageLoaderWrapper::setROI(int minX, int maxX, int minY, int maxY) {
    if (minX == -1)
        minX = 0;
    if (minY == -1)
        minY = 0;
    if (maxX == -1)
        maxX = _baseImageLoader->getXSize() - 1;
    if (maxY == -1)
        maxY = _baseImageLoader->getYSize() - 1;
    
    _minX = clip(minX, 0, _baseImageLoader->getXSize() - 1);
    _maxX = clip(maxX, 0, _baseImageLoader->getXSize() - 1);
    _minY = clip(minY, 0, _baseImageLoader->getYSize() - 1);
    _maxY = clip(maxY, 0, _baseImageLoader->getYSize() - 1);
    
    if ((_minX != 0) || (_minY != 0) || (_maxX != _baseImageLoader->getXSize() - 1) || (_maxY != _baseImageLoader->getYSize() - 1)) {
        _haveCustomROI = true;
    } else {
        _haveCustomROI = false;
    }
}

ImageLoaderSPE::ImageLoaderSPE(std::string rhs) {
	this->filePath = rhs;
	
	header_length = 4100;
	
	file.open(this->filePath.c_str(), std::ios::binary | std::ios::in);
	if (file.fail() == 1) {
		std::string error ("Unable to open the file at ");
		error += this->filePath;
		throw CANNOT_OPEN_FILE(error);
	}
	
	parse_header_information();
}

ImageLoaderSPE::~ImageLoaderSPE() {
	if (file.is_open() == 1)
		file.close();
}

void ImageLoaderSPE::parse_header_information() {
	// warning: only safe as long as both the writing
	// and reading systems are little-endian
	this->xSize = 0;
	this->ySize = 0;
	this->nImages = 0;
	this->storage_type = 0;
	
	file.seekg(42);
	file.read((char *)&(this->xSize), 2);
	
	file.seekg(656);
	file.read((char *)&(this->ySize), 2);
	
	file.seekg(1446);
	file.read((char *)&(this->nImages), 4);
	
	file.seekg(108);
	file.read((char *)&(this->storage_type), 2);
	
	switch (storage_type) {
		case 0:
			storage_type = STORAGE_TYPE_FP32;
			break;
		case 1:
			storage_type = STORAGE_TYPE_UINT32;
			break;
		case 2:
			storage_type = STORAGE_TYPE_INT16;
			break;
		case 3:
			storage_type = STORAGE_TYPE_UINT16;
			break;
		default:
			std::string error("Unable to determine the storage type used in ");
			error += this->filePath;
			throw CANNOT_DETERMINE_SPE_STORAGE_TYPE(error);
			break;
	}
	
	// was there an error sometime during this procedure that would have caused the reading to fail?
	if (file.fail() != 0) {
		std::string error;
		error = "Error parsing the header information in \"";
		error += this->filePath;
		error += "\" assuming the SPE format";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	file.seekg(0);
	this->checkForReasonableValues();
}

ImagePtr ImageLoaderSPE::readNextImage(int &indexOfImageThatWasRead) {
	uint64_t offset;
	
	uint64_t imageSize, pixelSize;
	
	// determine how big we have to make the single image buffer and the offset
	switch(storage_type) {
		case STORAGE_TYPE_FP32:	// 4 byte float
		case STORAGE_TYPE_UINT32:	// 4-byte long
			pixelSize = 4;
			break;
		case STORAGE_TYPE_INT16:	// 2 byte signed short
		case STORAGE_TYPE_UINT16:	// 2 byte unsigned short
			pixelSize = 2;
			break;
		default:
			std::string error("Unable to determine the storage type used in ");
			error += this->filePath;
			throw CANNOT_DETERMINE_SPE_STORAGE_TYPE(error);
			break;
	}
	
	imageSize = pixelSize * xSize * ySize;
	
    std::vector<char> imageBuffer(imageSize);
	ImagePtr image(GetRecycledMatrix((int)xSize, (int)ySize), FreeRecycledMatrix);
	
	{
		boost::lock_guard<boost::mutex> locker(loadImagesMutex);
		if (this->nextImageToRead >= this->nImages)
			throw std::runtime_error("requested more images than there are in the file");
		
		offset = this->header_length + this->nextImageToRead * imageSize;
		
		file.seekg(offset);
		file.read(imageBuffer.data(), imageSize);
		if (file.fail() != 0) {
			std::string error;
			error = "Error trying to read image data from \"";
			error += this->filePath;
			error += "\" assuming the SPE format";
			throw ERROR_READING_FILE_DATA(error);
		}
		indexOfImageThatWasRead = this->nextImageToRead;
		this->nextImageToRead += 1;
	}
	
	switch(storage_type) {
		case STORAGE_TYPE_FP32:	// 4-byte float
            CopyBufferToImage<float>(imageBuffer, image);
			break;
            
		case STORAGE_TYPE_UINT32:	// 4-byte long
			CopyBufferToImage<uint32_t>(imageBuffer, image);
			break;
            
		case STORAGE_TYPE_INT16:	// 2-byte signed short
			CopyBufferToImage<int16_t>(imageBuffer, image);
			break;
			
		case STORAGE_TYPE_UINT16: // 2-byte unsigned short
			CopyBufferToImage<uint16_t>(imageBuffer, image);
			break;
            
		default:
			std::string error("Unable to determine the storage type used in ");
			error += this->filePath;
			throw CANNOT_DETERMINE_SPE_STORAGE_TYPE(error);
			break;
	}
	
	return image;
}

ImageLoaderAndor::ImageLoaderAndor(std::string rhs) {
	this->filePath = rhs;
	
	file.open(this->filePath.c_str(), std::ios::binary | std::ios::in);
	if (file.fail() == 1) {
		std::string error ("Error opening the file at ");
		error += this->filePath;
		throw CANNOT_OPEN_FILE(error);
	}
	
	parse_header_information();
}

ImageLoaderAndor::~ImageLoaderAndor() {
	if (file.is_open() == 1) {
		file.close();
	}
}

void ImageLoaderAndor::parse_header_information() {
	int temp;
	int xBinning, yBinning;
	int frameXSize, frameYSize, nImagesHeader;
	int xStart, yStart, xEnd, yEnd;
	int result;
    std::unique_ptr<char[]> headerBuffer(new char[500001]);
    std::unique_ptr<char[]> singleLineBuffer(new char[4096]);
	std::string singleLine;
	
	storage_type = STORAGE_TYPE_FP32;
	
	this->file.read(headerBuffer.get(), 500000);
	if (this->file.good() != 1)
		throw std::runtime_error(std::string("Error encountered assuming the Andor format on the file at ") + this->filePath);
	
	headerBuffer[500000] = '\0';
	// look for and replace any intermediate nul characters
	for (int i = 0; i < 500000; ++i) {
		if (headerBuffer[i] == '\0') {
			headerBuffer[i] = '1';
		}
	}
	
	std::stringstream ss(headerBuffer.get(), std::ios::in);
	ss.getline(singleLineBuffer.get(), 4096);
	singleLine = singleLineBuffer.get();
	
	if (singleLine.find("Andor Technology Multi-Channel File") == std::string::npos) {
		throw std::runtime_error(std::string("the file at ") + this->filePath + "does not appear to be an Andor data file");
	}
	
	for (size_t i = 0;; ++i) {
		ss.getline(singleLineBuffer.get(), 4096);
		singleLine = singleLineBuffer.get();
		if ((ss.eof() == 1) || (ss.fail() == 1))
			throw std::runtime_error(std::string("premature end-of-file encountered assuming the Andor format on the file at ") + this->filePath);
		if (singleLine.find("Pixel number65") != std::string::npos) {
			// this is the first line containing info required to read the data
			break;
		}
	}
	
	result = sscanf(singleLine.c_str(), "Pixel number%d 1 %d %d 1 %d", &temp, &frameYSize, &frameXSize, &nImagesHeader);
	if (result != 4)
		throw std::runtime_error(std::string("an error occured parsing the file assuming the Andor format"));
	
	ss.getline(singleLineBuffer.get(), 4096);
	singleLine = singleLineBuffer.get();
	result = sscanf(singleLine.c_str(), "%d %d %d %d %d %d %d", &temp, &xStart, &yEnd, &xEnd, &yStart, &xBinning, &yBinning);
	if (result != 7)
		throw std::runtime_error(std::string("an error occured parsing the file assuming the Andor format"));
	
	this->xSize = (xEnd - xStart + 1) / xBinning;
	this->ySize = (yEnd - yStart + 1) / yBinning;	// integer division
    this->nImages = nImagesHeader;
	
    // apparently some Andor files have a bunch of empty lines next, and after that lines that may contain timestamps or trigger information.
    // The number of lines with timestamps is equal to the number of images in the file. In all the files I've seen the timestamps always
    // contain simply '0' and some whitespace. So we will now read lines until we find line containing '0', and assume that this is the first
    // line containing the timestamp information.
    for (;;) {
        ss.getline(singleLineBuffer.get(), 4096);
        if ((ss.eof() == 1) || (ss.fail() == 1))
			throw std::runtime_error(std::string("premature end-of-file encountered assuming the Andor format on the file at ") + this->filePath);
        singleLine = singleLineBuffer.get();
        if (singleLine.find('0') != std::string::npos)
            break;
    }
    
    // previously, this code assumed that there would be lines containing '0' next, with one line for each image in the file. But I've now seen
    // files where there is more than one such line for each image. I don't know what these lines are for, so now we just skip over any lines
    // that contain only whitespace and a single number.
    for (int i = 0; ; ++i) {
        ss.getline(singleLineBuffer.get(), 4096);
        if ((ss.eof() == 1) || (ss.fail() == 1))
			break;
        singleLine = singleLineBuffer.get();
        size_t firstNonWhiteSpace = singleLine.find_first_not_of(" \t");
        size_t lastNonWhiteSpace = singleLine.find_last_not_of(" \t");
        if ((firstNonWhiteSpace == std::string::npos) || (lastNonWhiteSpace == std::string::npos) || (firstNonWhiteSpace != lastNonWhiteSpace))
            break;
        char firstNonWhite = singleLine[firstNonWhiteSpace];
        if ((firstNonWhite < '0') || (firstNonWhite > '9'))
            break;
        this->header_length = ss.tellg();
    }
	
	// did some error happen while reading the file?
	if (ss.bad() == 1) {
		throw ERROR_READING_FILE_DATA(std::string("Error parsing the header information in \"") + this->filePath + "\" assuming the Andor format");
	}
	
	this->checkForReasonableValues();
}

ImagePtr ImageLoaderAndor::readNextImage(int &indexOfImageThatWasRead) {
	uint64_t offset;
	uint64_t nBytesInImage = xSize * ySize * sizeof(float);
	
    std::vector<char> imageBuffer(nBytesInImage);
	ImagePtr image (GetRecycledMatrix((int)xSize, (int)ySize), FreeRecycledMatrix);
	
	{
		boost::lock_guard<boost::mutex> locker(loadImagesMutex);
		if (this->nextImageToRead >= nImages)
			throw IMAGE_INDEX_BEYOND_N_IMAGES(std::string("Requested more images than there are in the file"));
		
		offset = this->header_length + this->nextImageToRead * nBytesInImage;
		
		file.seekg(offset);
		file.read(imageBuffer.data(), nBytesInImage);
		if (file.fail() != 0) {
			std::string error;
			error = "Error trying to read image data from \"";
			error += this->filePath;
			error += "\" assuming the Andor format";
			throw ERROR_READING_FILE_DATA(error);
		}
		indexOfImageThatWasRead = this->nextImageToRead;
		this->nextImageToRead += 1;
	}
	
	CopyBufferToImage<float>(imageBuffer, image);
	
	return image;
}

std::map<std::string, ImageLoaderHamamatsu::ImageOffsets> ImageLoaderHamamatsu::_offsetsMap;

ImageLoaderHamamatsu::ImageLoaderHamamatsu(std::string rhs) {
	this->filePath = rhs;
	
	file.open(this->filePath.c_str(), std::ios::binary | std::ios::in);
	if (file.fail() == 1) {
		std::string error ("Unable to open the file at ");
		error += this->filePath;
		throw CANNOT_OPEN_FILE(error);
	}
	
	parse_header_information();
}

ImageLoaderHamamatsu::~ImageLoaderHamamatsu() {
	if (file.is_open() == 1)
		file.close();
}

class ImageLoaderHamamatsu_HeaderStructure {
public:
	uint16 magic;
	uint16 commentLength;
	uint16 xSize;
	uint16 ySize;
	uint16 xBinning;	// uncertain
	uint16 yBinning;	// uncertain
	uint16 storageFormat;
	uint32 nImages;
	uint16 nChannels;
	uint16 channel;		// uncertain
	double timeStamp;
	uint32 marker;
	uint32 misc;		// function unknown
};

void ImageLoaderHamamatsu::parse_header_information() {
	
	// do we already have information on the offsets of the images for this file?
	if (_offsetsMap.count(this->filePath) != 0) {
		// we have information, now check if the timestamp is still okay
		if (GetLastModificationTime(this->filePath) != _offsetsMap[this->filePath].modificationTime) {
			// looks like the file was modified
			// delete the offset information, it will be recreated below
			_offsetsMap.erase(this->filePath);
		}
	}
	
	// now create the offsets map if needed (if it did not exist or was deleted above)
	if (_offsetsMap.count(this->filePath) == 0) {
		ImageLoaderHamamatsu::ImageOffsets imageOffsets;
		ImageLoaderHamamatsu_HeaderStructure header;
		
		file.seekg(0);
		
		int nImagesRead = 0;
		int nImagesInFile = -1;
		int thisImageXSize, thisImageYSize;
		uint64_t nextHeaderOffset = 0;
		do {
			file.seekg(nextHeaderOffset);
			
			file.read((char *)&(header.magic), sizeof(header.magic));
			file.read((char *)&(header.commentLength), sizeof(header.commentLength));
			file.read((char *)&(header.xSize), sizeof(header.xSize));
			file.read((char *)&(header.ySize), sizeof(header.ySize));
			file.read((char *)&(header.xBinning), sizeof(header.xBinning));
			file.read((char *)&(header.yBinning), sizeof(header.yBinning));
			file.read((char *)&(header.storageFormat), sizeof(header.storageFormat));
			file.read((char *)&(header.nImages), sizeof(header.nImages));
			file.read((char *)&(header.nChannels), sizeof(header.nChannels));
			file.read((char *)&(header.channel), sizeof(header.channel));
			file.read((char *)&(header.timeStamp), sizeof(header.timeStamp));
			file.read((char *)&(header.marker), sizeof(header.marker));
			file.read((char *)&(header.misc), sizeof(header.misc));
			
			if (nImagesInFile == -1)
				nImagesInFile = header.nImages;
			thisImageXSize = header.xSize;
			thisImageYSize = header.ySize;
			
			if (header.storageFormat != 2) {	// not UINT16
				std::string error;
				error = "The file at \"";
				error += this->filePath;
				error += "\" specifies that it doesn't use UINT16 for storage. This usually means that the manufacturer's software corrupted the file.";
				throw ERROR_READING_FILE_DATA(error);
			}
			
			// was there an error reading the file?
			if (file.fail() != 0) {
				std::string error;
				error = "Error parsing the header information in \"";
				error += this->filePath;
				error += "\" assuming the Hamamatsu format";
				throw ERROR_READING_FILE_DATA(error);
			}
			
			imageOffsets.offsets.push_back(nextHeaderOffset + header.commentLength + 64);
			nextHeaderOffset += header.commentLength + 64 + thisImageXSize * thisImageYSize * sizeof(uint16_t);
			
			nImagesRead += 1;
			
			#ifdef WITH_IGOR
			if ((nImagesRead % 20) == 0) {
				int abort = SpinProcess();
				if (abort)
					throw USER_ABORTED("user abort");
			}
			#endif
		} while (nImagesRead < nImagesInFile);
		
		// store the obtained offsets and information
		imageOffsets.xSize = thisImageXSize;
		imageOffsets.ySize = thisImageYSize;
		imageOffsets.modificationTime = GetLastModificationTime(this->filePath);
		
		_offsetsMap[this->filePath] = imageOffsets;
		
		file.seekg(0);
	}
	
	this->nImages = _offsetsMap[this->filePath].offsets.size();
	this->xSize = _offsetsMap[this->filePath].xSize;
	this->ySize = _offsetsMap[this->filePath].ySize;
	this->_offsets = _offsetsMap[this->filePath].offsets;
	this->storage_type = STORAGE_TYPE_UINT16;
	
	this->checkForReasonableValues();
}

ImagePtr ImageLoaderHamamatsu::readNextImage(int &indexOfImageThatWasRead) {
	uint64_t offset;
	uint64_t imageSize = xSize * ySize * 2; // assume a 16-bit format
    
    std::vector<char> imageBuffer(imageSize);
	ImagePtr image (GetRecycledMatrix((int)xSize, (int)ySize), FreeRecycledMatrix);
	
	{
		boost::lock_guard<boost::mutex> locker(loadImagesMutex);
		if (this->nextImageToRead >= nImages)
			throw IMAGE_INDEX_BEYOND_N_IMAGES(std::string("Requested more images than there are in the file"));
		
		offset = _offsets.at(nextImageToRead);
		
		file.seekg(offset);
		file.read(imageBuffer.data(), imageSize);
		if (file.fail() != 0) {
			std::string error;
			error = "Error reading image data from \"";
			error += this->filePath;
			error += "\" assuming the Hamamatsu format";
			throw ERROR_READING_FILE_DATA(error);
		}
		indexOfImageThatWasRead = this->nextImageToRead;
		this->nextImageToRead += 1;
	}
	
	CopyBufferToImage<uint16_t>(imageBuffer, image);
	
	return image;
}

ImageLoaderPDE::ImageLoaderPDE(std::string rhs) {
	this->filePath = rhs;
	
	file.open(this->filePath.c_str(), std::ios::binary | std::ios::in);
	if (file.fail() == 1) {
		std::string error ("Unable to open the file at ");
		error += this->filePath;
		throw CANNOT_OPEN_FILE(error);
	}
	
	parse_header_information();
}

ImageLoaderPDE::~ImageLoaderPDE() {
	if (file.is_open() == 1) {
		file.close();
	}
}

void ImageLoaderPDE::parse_header_information() {
	
	PDEFormatHeader header;
	
	file.seekg(0);
	file.read((char *)&(header.magic), sizeof(header.magic));
	file.read((char *)&(header.version), sizeof(header.version));
	file.read((char *)&(header.nImages), sizeof(header.nImages));
	file.read((char *)&(header.xSize), sizeof(header.xSize));
	file.read((char *)&(header.ySize), sizeof(header.ySize));
	file.read((char *)&(header.storageFormat), sizeof(header.storageFormat));
	
	if (file.fail() != 0) {
		std::string error;
		error = "Error parsing the header information in \"";
		error += this->filePath;
		error += "\" assuming the simple image format";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	if (header.magic != 27)
		throw std::runtime_error("The data is either not in the PDE format or of a different endianness");
	if (header.version != 1)
		throw std::runtime_error("Unsupported version of the PDE format.");
	
	this->nImages = header.nImages;
	this->xSize = header.xSize;
	this->ySize = header.ySize;
	this->storage_type = header.storageFormat;
	this->header_length = 24;
	
	this->checkForReasonableValues();
}

ImagePtr ImageLoaderPDE::readNextImage(int &indexOfImageThatWasRead) {
	uint64_t offset;
	
	size_t imageSize = NBytesInImage(this->xSize, this->ySize, this->storage_type);
	
    std::vector<char> imageBuffer(imageSize);
	{
		boost::lock_guard<boost::mutex> locker(loadImagesMutex);
		if (this->nextImageToRead >= this->nImages)
			throw IMAGE_INDEX_BEYOND_N_IMAGES(std::string("Requested more images than there are in the file"));
		
		offset = this->header_length + this->nextImageToRead * imageSize;
		
		file.seekg(offset);
		file.read(imageBuffer.data(), imageSize);
		if (file.fail() != 0) {
			std::string error;
			error = "Error reading image data from \"";
			error += this->filePath;
			error += "\" assuming the simple image format";
			throw ERROR_READING_FILE_DATA(error);
		}
		indexOfImageThatWasRead = this->nextImageToRead;
		this->nextImageToRead += 1;
	}
	
	
	ImagePtr image = BufferWithFormatToImage(imageBuffer, this->xSize, this->ySize, this->storage_type, 1);
	return image;
}

std::map<std::string, ImageLoaderTIFF::ImageOffsets> ImageLoaderTIFF::_offsetsMap;

ImageLoaderTIFF::ImageLoaderTIFF(std::string rhs) {
	TIFFSetWarningHandler(NULL);
	this->filePath = rhs;
	
	tiff_file = NULL;
	tiff_file = TIFFOpen(this->filePath.c_str(), "rm");
	if (tiff_file == NULL) {
		std::string error ("Unable to open the file at ");
		error += this->filePath;
		throw CANNOT_OPEN_FILE(error);
	}
	
	parse_header_information();
}


ImageLoaderTIFF::~ImageLoaderTIFF() {
	
	if (tiff_file != NULL) {
		TIFFClose(tiff_file);
	}
	
}


void ImageLoaderTIFF::parse_header_information() {
	int result;
	uint32_t result_uint32;
    
    result = TIFFSetDirectory(tiff_file, 0);
    if (result != 1)
        throw std::runtime_error("Error reading from the file");
    
    _extractSampleFormat();
    
    // do we already have information on the offsets of the images for this file?
	if (_offsetsMap.count(this->filePath) != 0) {
		// we have information, now check if the timestamp is still okay
		if (GetLastModificationTime(this->filePath) != _offsetsMap[this->filePath].modificationTime) {
			// looks like the file was modified
			// delete the offset information, it will be recreated below
			_offsetsMap.erase(this->filePath);
		} else {
            // we still have valid offsets
            _directoryOffsets = _offsetsMap[this->filePath].offsets;
        }
	}
    
    // now create the offsets map if needed (if it did not exist or was deleted above)
	if (_offsetsMap.count(this->filePath) == 0) {
		ImageLoaderTIFF::ImageOffsets imageOffsets;
        
        // how many images are there in the file?
        // because there is a possibility of a file with different
        // subfile types (e.g. Zeiss LSM files), we need to go over all of them
        // and look at which subfiles are the ones we're interested in
        imageOffsets.modificationTime = GetLastModificationTime(this->filePath);
        int nImagesChecked = 0;
        do {
            result = TIFFGetField(tiff_file, TIFFTAG_SUBFILETYPE, &result_uint32);

            if (((result_uint32 == FILETYPE_REDUCEDIMAGE) || (result_uint32 == FILETYPE_MASK)) && (result == 1)) {
                // this is one of the subtypes that we don't support
                // do nothing with it
            } else {
                // if we're here then the image is appropriate, store the offset to its IFD
                int64_t ifdOffset = TIFFCurrentDirOffset(tiff_file);
                _directoryOffsets.push_back(ifdOffset);
            }
            
            nImagesChecked += 1;
            #ifdef WITH_IGOR
			if ((nImagesChecked % 20) == 0) {
				int abort = SpinProcess();
				if (abort)
					throw USER_ABORTED("user abort");
			}
            #endif
            
        } while (TIFFReadDirectory(tiff_file) == 1);
        
        imageOffsets.offsets = _directoryOffsets;
        
        _offsetsMap[this->filePath] = imageOffsets;
    }
	
	this->nImages = _directoryOffsets.size();
    
	this->checkForReasonableValues();
}

void ImageLoaderTIFF::_extractSampleFormat() {
    int result;
    uint16_t result_uint16;
	uint32_t result_uint32;
	int isInt;
	size_t bitsPerPixel;
    
    // is the image in grayscale format?
    result = TIFFGetField(tiff_file, TIFFTAG_PHOTOMETRIC, &result_uint16);
    if (result != 1) {
        std::string error;
        error = "The image at\"";
        error += this->filePath;
        error += "\" is not a grayscale image";
        throw ERROR_READING_FILE_DATA(error);
    }
    if ((result_uint16 != 0) && (result_uint16 != 1)) {	// not a grayscale image
        std::string error;
        error = "The image at\"";
        error += this->filePath;
        error += "\" is not a grayscale image";
        throw ERROR_READING_FILE_DATA(error);
    }
    
    // is it a binary image?
    result = TIFFGetField(tiff_file, TIFFTAG_BITSPERSAMPLE, &result_uint16);
    if (result != 1) {
        std::string error;
        error = "The image at\"";
        error += this->filePath;
        error += "\" is not a grayscale image";
        throw ERROR_READING_FILE_DATA(error);
    }
    if (result_uint16 < 4) {	// 4 is the minimum number of bits allowed for grayscale images in the tiff specification, so this is a bilevel image
        std::string error;
        error = "The image at\"";
        error += this->filePath;
        error += "\" is not a grayscale image";
        throw ERROR_READING_FILE_DATA(error);
    }
    bitsPerPixel = (unsigned int)result_uint16;
	
	// is the data in unsigned integer or floating point format?
	result = TIFFGetField(tiff_file, TIFFTAG_SAMPLEFORMAT, &result_uint16);
	if (result != 1) {	// if the field does not exist then we assume that it is integer format
		result_uint16 = 1;
	}
	switch (result_uint16) {
		case 1:
			isInt = 1;
			break;
		case 3:
			isInt = 0;
			break;
		default:
			std::string error;
			error = "The SampleFormat of the image at\"";
			error += this->filePath;
			error += "\" is unknown";
			throw ERROR_READING_FILE_DATA(error);
			break;
	}
	
	if (isInt == 1) {
		switch (bitsPerPixel) {
			case 4:
				storage_type = STORAGE_TYPE_UINT4;
				break;
			case 8:
				storage_type = STORAGE_TYPE_UINT8;
				break;
			case 16:
				storage_type = STORAGE_TYPE_UINT16;
				break;
			case 32:
				storage_type = STORAGE_TYPE_UINT32;
				break;
			default:
				std::string error;
				error = "The SampleFormat of the image at\"";
				error += this->filePath;
				error += "\" is unknown";
				throw ERROR_READING_FILE_DATA(error);
				break;
		}
	} else {	// the image is not in integer format but is a floating point
		switch (bitsPerPixel) {
			case 32:
				storage_type = STORAGE_TYPE_FP32;
				break;
			case 64:
				storage_type = STORAGE_TYPE_FP64;
				break;
			default:
				std::string error;
				error = "The SampleFormat of the image at\"";
				error += this->filePath;
				error += "\" is unknown";
				throw ERROR_READING_FILE_DATA(error);
				break;
		}
	}
	
	// what is the x size?
	result = TIFFGetField(tiff_file, TIFFTAG_IMAGEWIDTH, &result_uint32);
	if (result != 1) {
		std::string error;
		error = "The image at\"";
		error += this->filePath;
		error += "\" does not specify a width";
		throw ERROR_READING_FILE_DATA(error);
	}
	xSize = (size_t)(result_uint32);
	
	// what is the y size?
	result = TIFFGetField(tiff_file, TIFFTAG_IMAGELENGTH, &result_uint32);
	if (result != 1) {
		std::string error;
		error = "The image at\"";
		error += this->filePath;
		error += "\" does not specify a height";
		throw ERROR_READING_FILE_DATA(error);
	}
	ySize = (size_t)(result_uint32);
}

ImagePtr ImageLoaderTIFF::readNextImage(int &indexOfImageThatWasRead) {
	int result;
	
	boost::lock_guard<boost::mutex> lock(loadImagesMutex);
	
	std::shared_ptr<void> single_scanline_buffer (_TIFFmalloc(TIFFScanlineSize(tiff_file)), _TIFFfree);
	if (single_scanline_buffer.get() == NULL) {
		throw std::bad_alloc();
	}
	
	if (this->nextImageToRead >= this->nImages)
		throw IMAGE_INDEX_BEYOND_N_IMAGES(std::string("Requested more images than there are in the file"));
	
	// advance the active TIFF directory to that needed to read the next image
    result = TIFFSetSubDirectory(tiff_file, _directoryOffsets.at(this->nextImageToRead));
    if (result != 1) {
        std::string error;
        error = "Unable to set the directory for the image at\"";
        error += this->filePath;
        error += "\"";
        throw ERROR_READING_FILE_DATA(error);
    }
	
	ImagePtr image(GetRecycledMatrix((int)xSize, (int)ySize), FreeRecycledMatrix);
	
	for (int j = 0; j < ySize; ++j) {
		result = TIFFReadScanline(tiff_file, single_scanline_buffer.get(), j, 0);	// sample is ignored
		if (result != 1) {
			std::string error;
			error = "Unable to read a scanline from the image at\"";
			error += this->filePath;
			error += "\"";
			throw ERROR_READING_FILE_DATA(error);
		}
		
		switch (storage_type) {	// handle the different possibilities (floating, integer) and variable sizes
			case STORAGE_TYPE_UINT4:
			{
				uint8_t *uint8Ptr = (uint8_t *)single_scanline_buffer.get();
				for (int k = 0; k < xSize; k+=2) {
					(*image)(k, j) = (double)(0x0000000F & (*uint8Ptr));
					(*image)(k + 1, j) = (double)((0x000000F0 & (*uint8Ptr)) / 16);			  
					++uint8Ptr;
				}
				break;
			}
			case STORAGE_TYPE_UINT8:
			{
				uint8_t *uint8Ptr = (uint8_t *)single_scanline_buffer.get();
				for (int k = 0; k < xSize; ++k) {
					(*image)(k, j) = (double)(*uint8Ptr);
					++uint8Ptr;
				}
				break;
			}
			case STORAGE_TYPE_UINT16:
			{
				uint16_t *uint16tPtr = (uint16_t *)single_scanline_buffer.get();
				for (int k = 0; k < xSize; ++k) {
					(*image)(k, j) = (double)(*uint16tPtr);
					++uint16tPtr;
				}
				break;
			}
			case STORAGE_TYPE_UINT32:
			{
				uint32_t *uint32tPtr = (uint32_t *)single_scanline_buffer.get();
				for (int k = 0; k < xSize; ++k) {
					(*image)(k, j) = (double)(*uint32tPtr);
					++uint32tPtr;
				}
				break;
			}
			case STORAGE_TYPE_FP32:
			{
				float *floatPtr = (float *)single_scanline_buffer.get();
				for (int k = 0; k < xSize; ++k) {
					(*image)(k, j) = (double)(*floatPtr);
					++floatPtr;
				}
				break;
			}
			case STORAGE_TYPE_FP64:
			{
				double *doublePtr = (double *)single_scanline_buffer.get();
				for (int k = 0; k < xSize; ++k) {
					(*image)(k, j) = (double)(*doublePtr);
					++doublePtr;
				}
				break;
			}
			default:
				std::string error;
				error = "Invalid floating point data size for the image at\"";
				error += this->filePath;
				error += "\"";
				throw ERROR_READING_FILE_DATA(error);
				break;
		}
	}
	
	indexOfImageThatWasRead = this->nextImageToRead;
	this->nextImageToRead += 1;
	
	return image;
}


ImageLoaderMultiFileTIFF::ImageLoaderMultiFileTIFF(std::string filePath) {
	// assume that filePath consists of something such as the following
	// /path/to/folder/XXXX123.YYY where X is character (but the last one is not a number)
	// and Y is the extension. "/path/to/folder/XXXX" is the base file path, YYY is the
	// extension, and the number of digits in between is assumed to be constant for every
	// file in the sequence
	
	size_t stringLength = filePath.length();
	size_t extensionStartIndex = filePath.rfind('.');
	size_t lengthOfExtension;
	if (extensionStartIndex == std::string::npos) {
		// no extension. Maybe this is OK but it probably is an error.
		// pretend that everything is fine.
		extension.clear();
	} else {
		extension = filePath.substr(extensionStartIndex);
	}
	lengthOfExtension = extension.length();
	
	// now determine how many digits there are that indicate the position of each file in the
	// frame sequence. Assume that all digits that precede the extension are part of this index.
	int nIndexDigits = 0;
	for (int i = extensionStartIndex - 1; i >= 0; --i) {
		if (isdigit(filePath[i])) {
			nIndexDigits += 1;
		} else {
			break;
		}
	}
	if (nIndexDigits == 0) {
		std::string error("Unable to determine the numbering format for the multi-tiff file at \"");
		error += filePath;
		error += "\"";
		throw ERROR_READING_FILE_DATA(error);
	}
	if (nIndexDigits > 9) {	// arbitrary limit
		std::string error("The numeric index of the images at \"");
		error += filePath;
		error += "\" appears to contain too many digits to be realistic";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	nDigitsInNumber = nIndexDigits;
	
	// everything coming before the index digits is the base file path
	this->baseFilePath = filePath.substr(0, stringLength - lengthOfExtension - nDigitsInNumber);
	std::string indexOfThisImageStr = filePath.substr(stringLength - lengthOfExtension - nDigitsInNumber, nDigitsInNumber);
	
	// try to load the provided image to get the metadata
	ImageLoaderTIFF imageLoaderTIFF(filePath);
	this->xSize = imageLoaderTIFF.getXSize();
	this->ySize = imageLoaderTIFF.getYSize();
	this->storage_type = imageLoaderTIFF.getStorageType();
	
	// now figure out how many images there are in the dataset
	// get the index of the current image
	int indexOfThisImage;
	int nItemsMatched = sscanf(indexOfThisImageStr.c_str(), "%d", &indexOfThisImage);
	if (nItemsMatched != 1) {
		std::string error("Error parsing the numeric index for the image at \"");
		error += filePath;
		error += "\"";
		throw ERROR_READING_FILE_DATA(error);
	}
	
	std::pair<int, int> firstAndLastIndices = findFirstAndLastValidImageIndices(indexOfThisImage);
	this->firstImageIndex = firstAndLastIndices.first;
	this->nImages = firstAndLastIndices.second - firstAndLastIndices.first + 1;
	
	this->checkForReasonableValues();
}

ImagePtr ImageLoaderMultiFileTIFF::readNextImage(int &indexOfImageThatWasRead) {
	boost::lock_guard<boost::mutex> lock(loadImagesMutex);
	
	if (this->nextImageToRead >= this->nImages)
		throw IMAGE_INDEX_BEYOND_N_IMAGES(std::string("Requested more images than there are in the file"));
	
	std::string filePath = getFilePathForImageAtIndex(this->nextImageToRead);
	ImageLoaderTIFF imageLoaderTIFF(filePath);
	int dummy;
	ImagePtr image = imageLoaderTIFF.readNextImage(dummy);
	
	indexOfImageThatWasRead = this->nextImageToRead;
	this->nextImageToRead += 1;
	
	return image;
}

std::string ImageLoaderMultiFileTIFF::getFilePathForImageAtIndex(int index) {
	std::unique_ptr<char[]> imageIndexStr(new char[this->nDigitsInNumber + 1]);
	char formatString[10];
	sprintf(formatString, "%%0%dd", nDigitsInNumber);
	sprintf(imageIndexStr.get(), formatString, index);
	
	std::string filePath = this->baseFilePath;
	filePath += imageIndexStr.get();
	filePath += this->extension;
	return filePath;
}

bool ImageLoaderMultiFileTIFF::imageFileAtIndexExists(int index) {
	std::string thisFilePath = getFilePathForImageAtIndex(index);
	
	bool exists = true;
	try {
		ImageLoaderTIFF imageLoaderTIFF(thisFilePath);
	}
	catch (...) {
		exists = false;
	}
	
	return exists;
}

std::pair<int, int> ImageLoaderMultiFileTIFF::findFirstAndLastValidImageIndices(int knownValidImageIndex) {
	std::pair<int, int> firstAndLastIndices;
	
	assert(knownValidImageIndex >= 0);
	
	int trialIndex;
	int lower, upper, mid;
	
	// find the first image index that is valid
	if (imageFileAtIndexExists(0)) {
		// check if the first image is simply at index zero
		firstAndLastIndices.first = 0;
	} else {
		// find an invalid index less than knownValidImageIndex
		trialIndex = knownValidImageIndex;
		for (int delta = -1; ; delta *= 2) {
			trialIndex = trialIndex + delta;
			if (trialIndex < 0) {
				trialIndex = -1;	// the image at -1 certainly doesn't exist
				break;
			}
			if (!imageFileAtIndexExists(trialIndex))
				break;
		}
		
		// now we have a known bad and a known good index. Find the first index that is valid.
		lower = trialIndex;
		upper = knownValidImageIndex;
		for (;;) {
			if (upper - lower == 1) {
				firstAndLastIndices.first = upper;
				break;
			}
			
			mid = (lower + upper) / 2;
			if (imageFileAtIndexExists(mid)) {
				upper = mid;
			} else {
				lower = mid;
			}
		}
	}
	
	// now find the last index that is valid. First find an index that is invalid.
	trialIndex = knownValidImageIndex + 1;
	for (int delta = 1; ; delta *= 2) {
		trialIndex = trialIndex + delta;
		if (!imageFileAtIndexExists(trialIndex))
			break;
	}
	
	lower = knownValidImageIndex; upper = trialIndex;
	for (;;) {
		if (upper - lower == 1) {
			firstAndLastIndices.second = lower;
			break;
		}
		
		mid = (lower + upper) / 2;
		if (imageFileAtIndexExists(mid)) {
			lower = mid;
		} else {
			upper = mid;
		}
	}
	
	return firstAndLastIndices;
}

#ifdef WITH_IGOR
ImageLoaderIgor::ImageLoaderIgor(std::string waveName) {
	int waveType;
	
	// try to get images from an Igor wave
	
	this->igor_data_wave = FetchWaveUsingFullPath(waveName);
	
	size_t DimensionSizes[MAX_DIMENSIONS + 1];
	int numDimensions;
	int result;
	
	result = MDGetWaveDimensions(igor_data_wave, &numDimensions, (CountInt *)DimensionSizes);
	if (result != 0) {
		throw result;
	}
	if ((numDimensions != 2) && (numDimensions != 3)) {
		throw INCOMPATIBLE_DIMENSIONING;
	}
	
	xSize = (size_t)DimensionSizes[0];
	ySize = (size_t)DimensionSizes[1];
	nImages = (size_t)DimensionSizes[2];
	
	// special case: if the wave contains only a single image
    // then it is usually two-dimensional, that is, DimensionSizes[2] == 0
	// in that case nImages is still one
	if (DimensionSizes[2] == 0) {
		nImages = 1;
	}
	
	waveType = WaveType(igor_data_wave);
    // do not handle complex waves
    if (waveType & NT_CMPLX)
        throw NO_COMPLEX_WAVE;
    
	switch (waveType) {
		case NT_I8:
			storage_type = STORAGE_TYPE_INT8;
			break;
		case NT_I16:
			storage_type = STORAGE_TYPE_INT16;
			break;
		case NT_I32:
			storage_type = STORAGE_TYPE_INT32;
			break;
		case (NT_I8 | NT_UNSIGNED):
			storage_type = STORAGE_TYPE_UINT8;
			break;
		case (NT_I16 | NT_UNSIGNED):
			storage_type = STORAGE_TYPE_UINT16;
			break;
		case (NT_I32 | NT_UNSIGNED):
			storage_type = STORAGE_TYPE_UINT32;
			break;
		case NT_FP32:
			storage_type = STORAGE_TYPE_FP32;
			break;
		case NT_FP64:
			storage_type = STORAGE_TYPE_FP64;
			break;
		default:
			storage_type = STORAGE_TYPE_FP64;
	}
	
	this->checkForReasonableValues();
}

ImagePtr ImageLoaderIgor::readNextImage(int &index) {
	int err;
	int waveType = WaveType(this->igor_data_wave);
	size_t waveDataOffset;
	char* startOfWaveData;
	size_t nPixels = this->xSize * this->ySize;
	
	// allocate a new image
	ImagePtr image (GetRecycledMatrix((int)xSize, (int)ySize), FreeRecycledMatrix);
	
	// get a pointer to the data in the wave
	err = MDAccessNumericWaveData(this->igor_data_wave, kMDWaveAccessMode0, (BCInt*)&waveDataOffset);
	if (err != 0) {
		throw err;
	}
	startOfWaveData = ((char*)(*this->igor_data_wave) + waveDataOffset);
	
	{
		boost::lock_guard<boost::mutex> lock(this->loadImagesMutex);
		if (this->nextImageToRead >= nImages)
			throw IMAGE_INDEX_BEYOND_N_IMAGES(std::string("Requested more images than there are in the file"));
		
		index = this->nextImageToRead;
		this->nextImageToRead += 1;
	}
	
	// get a pointer to the image
	double* imagePtr = image->data();
	
	switch (waveType) {
		case NT_FP32:
		{
			float* floatPtr = (float*)startOfWaveData;
			floatPtr += nPixels * index;
			for (size_t i = 0; i < nPixels; ++i) {
				imagePtr[i] = floatPtr[i];
			}
			break;
		}
		case NT_FP64:
		{
			double* doublePtr = (double*)startOfWaveData;
			doublePtr += nPixels * index;
			memcpy(imagePtr, doublePtr, nPixels * sizeof(double));
			break;
		}
		case NT_I8:
		{
			int8_t* int8Ptr = (int8_t*)startOfWaveData;
			int8Ptr += nPixels * index;
			for (size_t i = 0; i < nPixels; ++i) {
				imagePtr[i] = int8Ptr[i];
			}
			break;
		}
		case NT_I16:
		{
			int16_t* int16Ptr = (int16_t*)startOfWaveData;
			int16Ptr += nPixels * index;
			for (size_t i = 0; i < nPixels; ++i) {
				imagePtr[i] = int16Ptr[i];
			}
			break;
		}
		case NT_I32:
		{
			int32_t* int32Ptr = (int32_t*)startOfWaveData;
			int32Ptr += nPixels * index;
			for (size_t i = 0; i < nPixels; ++i) {
				imagePtr[i] = int32Ptr[i];
			}
			break;
		}
		case NT_I8 | NT_UNSIGNED:
		{
			uint8_t* uint8Ptr = (uint8_t*)startOfWaveData;
			uint8Ptr += nPixels * index;
			for (size_t i = 0; i < nPixels; ++i) {
				imagePtr[i] = uint8Ptr[i];
			}
			break;
		}
		case NT_I16 | NT_UNSIGNED:
		{
			uint16_t* uint16Ptr = (uint16_t*)startOfWaveData;
			uint16Ptr += nPixels * index;
			for (size_t i = 0; i < nPixels; ++i) {
				imagePtr[i] = uint16Ptr[i];
			}
			break;
		}
		case NT_I32 | NT_UNSIGNED:
		{
			uint32_t* uint32Ptr = (uint32_t*)startOfWaveData;
			uint32Ptr += nPixels * index;
			for (size_t i = 0; i < nPixels; ++i) {
				imagePtr[i] = uint32Ptr[i];
			}
			break;
		}
		default:
			throw std::runtime_error("Unknown or unsupported wavetype");
			break;
	}
	
	return image;
}
#endif // WITH_IGOR

#ifdef WITH_MATLAB
ImageLoaderMatlab::ImageLoaderMatlab(mxArray* matlabArray) {
	_matlabArray = matlabArray;
	
	mwSize nDims = mxGetNumberOfDimensions(matlabArray);
	if ((nDims < 2) || (nDims > 3))
		throw std::runtime_error("Expected 2D or 3D matrix");
	
	const mwSize* dimensionSizes = mxGetDimensions(matlabArray);
	xSize = static_cast<size_t>(dimensionSizes[0]);
	ySize = static_cast<size_t>(dimensionSizes[1]);
	if (nDims == 3)
		nImages = static_cast<size_t>(dimensionSizes[2]);
	else
		nImages = 1;
	
	this->storage_type = static_cast<int>(mxGetClassID(matlabArray));
	this->checkForReasonableValues();
}

ImagePtr ImageLoaderMatlab::readNextImage(int &index) {
	size_t nPixels = this->xSize * this->ySize;
	
	{
		boost::lock_guard<boost::mutex> lock(this->loadImagesMutex);
		if (this->nextImageToRead >= nImages)
			throw IMAGE_INDEX_BEYOND_N_IMAGES(std::string("Requested more images than there are in the file"));
		
		index = this->nextImageToRead;
		this->nextImageToRead += 1;
	}
	
	// allocate a new image
	ImagePtr image (GetRecycledMatrix((int)xSize, (int)ySize), FreeRecycledMatrix);
	
	char* dataPtr = reinterpret_cast<char*>(mxGetData(_matlabArray));
	
	switch (mxGetClassID(_matlabArray)) {
		case mxINT8_CLASS:
		{
			size_t offset = index * nPixels * sizeof(int8_t);
			CopyBufferToImage<int8_t>(dataPtr + offset, image);
			break;
		}
		case mxUINT8_CLASS:
		{
			size_t offset = index * nPixels * sizeof(uint8_t);
			CopyBufferToImage<uint8_t>(dataPtr + offset, image);
			break;
		}
		case mxINT16_CLASS:
		{
			size_t offset = index * nPixels * sizeof(int16_t);
			CopyBufferToImage<int16_t>(dataPtr + offset, image);
			break;
		}
		case mxUINT16_CLASS:
		{
			size_t offset = index * nPixels * sizeof(uint16_t);
			CopyBufferToImage<uint16_t>(dataPtr + offset, image);
			break;
		}
		case mxINT32_CLASS:
		{
			size_t offset = index * nPixels * sizeof(int32_t);
			CopyBufferToImage<int32_t>(dataPtr + offset, image);
			break;
		}
		case mxUINT32_CLASS:
		{
			size_t offset = index * nPixels * sizeof(uint32_t);
			CopyBufferToImage<uint32_t>(dataPtr + offset, image);
			break;
		}
		case mxINT64_CLASS:
		{
			size_t offset = index * nPixels * sizeof(int64_t);
			CopyBufferToImage<int64_t>(dataPtr + offset, image);
			break;
		}
		case mxUINT64_CLASS:
		{
			size_t offset = index * nPixels * sizeof(uint64_t);
			CopyBufferToImage<uint64_t>(dataPtr + offset, image);
			break;
		}
		case mxSINGLE_CLASS:
		{
			size_t offset = index * nPixels * sizeof(float);
			CopyBufferToImage<float>(dataPtr + offset, image);
			break;
		}
		case mxDOUBLE_CLASS:
		{
			size_t offset = index * nPixels * sizeof(double);
			CopyBufferToImage<double>(dataPtr + offset, image);
			break;
		}
		default:
			throw std::runtime_error("Unknown or unsupported matrix class ID");
			break;
	}
	
	return image;
}
#endif // WITH_MATLAB

ImageOutputWriter::ImageOutputWriter() {
	outputFilePath.assign("");
	nImagesWritten = 0;
}

PDEImageOutputWriter::PDEImageOutputWriter(const std::string &rhs,int overwrite, uint32_t storageType_rhs) {
	// if overwrite is non-zero then we overwrite any file that exists at the output path
	// if it is set to zero then we throw an error and abort instead of overwriting
	outputFilePath = rhs;
	
	if (overwrite == 0) {
		std::ifstream input_test;
		input_test.open(outputFilePath.c_str(), std::ios::in | std::ios::binary);
		if (input_test.good() == 1) {
			std::string error("The output file at ");
			error += this->outputFilePath;
			error += " already exists.";
			throw OUTPUT_FILE_ALREADY_EXISTS(error);	// escape without overwriting
		}
	}
	
	file.open(outputFilePath.c_str(), std::ios::binary | std::ios::out | std::ios::trunc);	// DANGER: OVERWRITING THE FILE
	if (file.fail() != 0) {
		std::string error("Cannot create an output file at ");
		error += this->outputFilePath;
		throw CANNOT_OPEN_OUTPUT_FILE(error);
	}
	
	this->xSize = 0;
	this->ySize = 0;
	this->nImagesWritten = 0;
	this->storageType = storageType_rhs;
	
	// write the header out in advance
	// before closing the file, after all images have been written,
	// we will come back to this and overwrite it with the correct
	// values
	this->WriteHeader();
}


PDEImageOutputWriter::~PDEImageOutputWriter() {
	if (file.is_open() != 0) {
		WriteHeader();
		file.close();
	}
}

void PDEImageOutputWriter::WriteHeader() {
	PDEFormatHeader header;
	
	header.magic = 27;
	header.version = 1;
	header.nImages = this->nImagesWritten;
	header.xSize = this->xSize;
	header.ySize = this->ySize;
	header.storageFormat = this->storageType;
	
	if (file.is_open()) {
		file.seekp(0, std::ios_base::beg);
		file.write((char *)&(header.magic), sizeof(header.magic));
		file.write((char *)&(header.version), sizeof(header.version));
		file.write((char *)&(header.nImages), sizeof(header.nImages));
		file.write((char *)&(header.xSize), sizeof(header.xSize));
		file.write((char *)&(header.ySize), sizeof(header.ySize));
		file.write((char *)&(header.storageFormat), sizeof(header.storageFormat));
		if (file.fail())
			throw std::runtime_error("Error trying to write to the PDE file");
	}
}

void PDEImageOutputWriter::write_image(ImagePtr imageToWrite) {
	
	// determine the size of the frames
	size_t currentXSize = imageToWrite->rows();
	size_t currentYSize = imageToWrite->cols();
	
	if (this->nImagesWritten == 0) {
		this->xSize = currentXSize;
		this->ySize = currentYSize;
	} else {
		if ((currentXSize != this->xSize) || (currentYSize != this->ySize)) {
			throw std::runtime_error("Tried to write an image with different dimensions to an SPE file");
		}
	}
	
    size_t nBytesToWrite = NBytesInImage(this->xSize, this->ySize, this->storageType);
    std::vector<char> imageBuffer;
    ImageToBufferWithFormat(imageToWrite, this->storageType, imageBuffer, 1);
	
	this->file.write(imageBuffer.data(), nBytesToWrite);
	++this->nImagesWritten;
}

void TIFFSampleFormatAndBitsPerSampleForFormat(const int dataFormat, int& sampleFormat, int& bitsPerSample) {
    // determine the output storage type
	switch (dataFormat) {
		case STORAGE_TYPE_INT4:
		case STORAGE_TYPE_UINT4:
		case STORAGE_TYPE_INT8:
			bitsPerSample = 8;
			sampleFormat = SAMPLEFORMAT_INT;
			break;
		case STORAGE_TYPE_UINT8:
			bitsPerSample = 8;
			sampleFormat = SAMPLEFORMAT_UINT;
			break;
		case STORAGE_TYPE_INT16:
			bitsPerSample = 16;
			sampleFormat = SAMPLEFORMAT_INT;
			break;
		case STORAGE_TYPE_UINT16:
			bitsPerSample = 16;
			sampleFormat = SAMPLEFORMAT_UINT;
			break;
		case STORAGE_TYPE_INT32:
			bitsPerSample = 32;
			sampleFormat = SAMPLEFORMAT_INT;
			break;
		case STORAGE_TYPE_UINT32:
			bitsPerSample = 32;
			sampleFormat = SAMPLEFORMAT_UINT;
			break;
		case STORAGE_TYPE_INT64:
			bitsPerSample = 64;
			sampleFormat = SAMPLEFORMAT_INT;
			break;
		case STORAGE_TYPE_UINT64:
			bitsPerSample = 64;
			sampleFormat = SAMPLEFORMAT_UINT;
			break;
		case STORAGE_TYPE_FP32:
			bitsPerSample = 32;
			sampleFormat = SAMPLEFORMAT_IEEEFP;
			break;
		case STORAGE_TYPE_FP64:
			bitsPerSample = 64;
			sampleFormat = SAMPLEFORMAT_IEEEFP;
			break;
		default:
			throw std::runtime_error("Unknown storage type requested for TIFF output");
			break;
	}
}

TIFFImageOutputWriter::TIFFImageOutputWriter(const std::string &rhs,int overwrite, int compression_rhs, int storageType_rhs) {
	// if overwrite is non-zero then we overwrite any file that exists at the output path
	// if it is set to zero then we throw an error and abort instead of overwriting
	
	TIFFSetWarningHandler(NULL);
	this->outputFilePath = rhs;
	this->compression = compression_rhs;
	this->storageType = storageType_rhs;
	
	if (overwrite == 0) {
		std::ifstream input_test;
		input_test.open(outputFilePath.c_str(), std::ios::in | std::ios::binary);
		if (input_test.good() == 1) {
			std::string error("The output file at ");
			error += this->outputFilePath;
			error += " already exists.";
			throw OUTPUT_FILE_ALREADY_EXISTS(error);	// escape without overwriting
		}
		input_test.close();
	}
	
	tiff_file = TIFFOpen(outputFilePath.c_str(), "w");
	if (tiff_file == NULL) {
		std::string error("Cannot create an output file at ");
		error += this->outputFilePath;
		throw CANNOT_OPEN_OUTPUT_FILE(error);
	}
	
	nImagesWritten = 0;
}


TIFFImageOutputWriter::~TIFFImageOutputWriter() {
	if (tiff_file != NULL) {
		TIFFClose(tiff_file);
		tiff_file = NULL;
	}
}



void TIFFImageOutputWriter::write_image(ImagePtr imageToWrite) {
	
	size_t xSize = imageToWrite->rows();
	size_t ySize = imageToWrite->cols();
	
	int result, sampleFormat, bitsPerSample;
	uint16_t current_uint16;
	uint32_t current_uint32;
	
	// determine the output storage type
    TIFFSampleFormatAndBitsPerSampleForFormat(this->storageType, sampleFormat, bitsPerSample);
	
	// make a scoped_array that will act as a single scanline buffer
	// make it a buffer of chars equal to the total number of bytes required
	std::unique_ptr<char[]> scanLine(new char[xSize * bitsPerSample / 8]);
	
	// make sure that all the image tags have the correct values
	result = TIFFSetField(tiff_file, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
	if (result != 1) {
		std::string error;
		error = "Unable to set the photometric type for the image at\"";
		error += outputFilePath;
		error += "\"";
		throw ERROR_WRITING_FILE_DATA(error);
	}
	
	current_uint32 = xSize;
	result = TIFFSetField(tiff_file, TIFFTAG_IMAGEWIDTH, current_uint32);
	if (result != 1) {
		std::string error;
		error = "Unable to set the image width for the image at\"";
		error += outputFilePath;
		error += "\"";
		throw ERROR_WRITING_FILE_DATA(error);
	}
	
	current_uint32 = ySize;
	result = TIFFSetField(tiff_file, TIFFTAG_IMAGELENGTH, current_uint32);
	if (result != 1) {
		std::string error;
		error = "Unable to set the image height for the image at\"";
		error += outputFilePath;
		error += "\"";
		throw ERROR_WRITING_FILE_DATA(error);
	}
	
	result = TIFFSetField(tiff_file, TIFFTAG_SAMPLEFORMAT, sampleFormat);
	if (result != 1) {
		std::string error;
		error = "Unable to set the SampleFormat for the image at\"";
		error += outputFilePath;
		error += "\"";
		throw ERROR_WRITING_FILE_DATA(error);
	}
	
	result = TIFFSetField(tiff_file, TIFFTAG_BITSPERSAMPLE, bitsPerSample);
	if (result != 1) {
		std::string error;
		error = "Unable to set the BitsPerSample for the image at\"";
		error += outputFilePath;
		error += "\"";
		throw ERROR_WRITING_FILE_DATA(error);
	}
	
	result = TIFFSetField(tiff_file, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
	if (result != 1) {
		std::string error;
		error = "Unable to set the SubFileType for the image at\"";
		error += outputFilePath;
		error += "\"";
		throw ERROR_WRITING_FILE_DATA(error);
	}
	
	current_uint16 = (uint16_t)nImagesWritten;
	result = TIFFSetField(tiff_file, TIFFTAG_PAGENUMBER, current_uint16);
	if (result != 1) {
		std::string error;
		error = "Unable to set the PageNumber for the image at\"";
		error += outputFilePath;
		error += "\"";
		throw ERROR_WRITING_FILE_DATA(error);
	}
	
	result = TIFFSetField(tiff_file, TIFFTAG_COMPRESSION, compression);
	if (result != 1) {
		std::string error;
		error = "Unable to set the compression method for the image at\"";
		error += outputFilePath;
		error += "\"";
		throw ERROR_WRITING_FILE_DATA(error);
	}
	
	size_t offset = 0;
	for (size_t j = 0; j < ySize; j++) {
		offset = 0;
		
		switch (this->storageType) {
			case STORAGE_TYPE_INT4:
			case STORAGE_TYPE_UINT4:
			case STORAGE_TYPE_INT8:
			{
				int8_t* charPtr = (int8_t *)scanLine.get();
				for (size_t i = 0; i < xSize; i++) {
					charPtr[i] = (int8_t)(*imageToWrite)(i, j);
				}
				break;
			}
			case STORAGE_TYPE_UINT8:
			{
				uint8_t* uCharPtr = (uint8_t *)scanLine.get();
				for (size_t i = 0; i < xSize; i++) {
					uCharPtr[i] = (uint8_t)(*imageToWrite)(i, j);
				}
				break;
			}
			case STORAGE_TYPE_INT16:
			{
				int16_t* int16Ptr = (int16_t *)scanLine.get();
				for (size_t i = 0; i < xSize; i++) {
					int16Ptr[i] = (int16_t)(*imageToWrite)(i, j);
				}
				break;
			}
			case STORAGE_TYPE_UINT16:
			{
				uint16_t* uInt16Ptr = (uint16_t *)scanLine.get();
				for (size_t i = 0; i < xSize; i++) {
					uInt16Ptr[i] = (uint16_t)(*imageToWrite)(i, j);
				}
				break;
			}
			case STORAGE_TYPE_INT32:
			{
				int32_t* int32Ptr = (int32_t *)scanLine.get();
				for (size_t i = 0; i < xSize; i++) {
					int32Ptr[i] = (int32_t)(*imageToWrite)(i, j);
				}
				break;
			}
			case STORAGE_TYPE_UINT32:
			{
				uint32_t* uInt32Ptr = (uint32_t *)scanLine.get();
				for (size_t i = 0; i < xSize; i++) {
					uInt32Ptr[i] = (uint32_t)(*imageToWrite)(i, j);
				}
				break;
			};
			case STORAGE_TYPE_INT64:
			{
				int64_t* int64Ptr = (int64_t *)scanLine.get();
				for (size_t i = 0; i < xSize; i++) {
					int64Ptr[i] = (int64_t)(*imageToWrite)(i, j);
				}
				break;
			}
			case STORAGE_TYPE_UINT64:
			{
				uint64_t* uInt64Ptr = (uint64_t *)scanLine.get();
				for (size_t i = 0; i < xSize; i++) {
					uInt64Ptr[i] = (uint64_t)(*imageToWrite)(i, j);
				}
				break;
			}
			case STORAGE_TYPE_FP32:
			{
				float* floatPtr = (float *)scanLine.get();
				for (size_t i = 0; i < xSize; i++) {
					floatPtr[i] = (float)(*imageToWrite)(i, j);
				}
				break;
			}
			case STORAGE_TYPE_FP64:
			{
				double* doublePtr = (double *)scanLine.get();
				for (size_t i = 0; i < xSize; i++) {
					doublePtr[i] = (double)(*imageToWrite)(i, j);
				}
				break;
			}
		}
		
		result = TIFFWriteScanline(tiff_file, (char *)scanLine.get(), j);
		if (result != 1) {
			std::string error;
			error = "There was an error writing a scanline for the image at\"";
			error += outputFilePath;
			error += "\"";
			throw ERROR_WRITING_FILE_DATA(error);
		}
		
	}
	
	result = TIFFWriteDirectory(tiff_file);
	if (result != 1) {
		std::string error;
		error = "Unable to write a directory for the image at\"";
		error += outputFilePath;
		error += "\"";
		throw ERROR_WRITING_FILE_DATA(error);
	}
	
	++this->nImagesWritten;
}

MultiFileTIFFImageOutputWriter::MultiFileTIFFImageOutputWriter(const std::string &baseOutputFilePath_rhs, int overwrite_rhs, bool compress, int storageType_rhs) :
    overwrite(overwrite_rhs),
    _compress(compress),
    storageType(storageType_rhs)
{
	baseOutputFilePath = baseOutputFilePath_rhs;
}

void MultiFileTIFFImageOutputWriter::write_image(ImagePtr imageToWrite) {
	char imageIndexStr[50];
	sprintf(imageIndexStr, "%06d", static_cast<int>(this->nImagesWritten));
	
	std::string thisFileName = this->baseOutputFilePath + std::string(imageIndexStr) + std::string(".tif");
	
    LocalizerTIFFImageOutputWriter singleFileOutputWriter(thisFileName, this->overwrite, _compress, this->storageType);
	//TIFFImageOutputWriter singleFileOutputWriter(thisFileName, this->overwrite, this->compression, this->storageType);
	singleFileOutputWriter.write_image(imageToWrite);
	
	++this->nImagesWritten;
}

LocalizerTIFFImageOutputWriter::LocalizerTIFFImageOutputWriter(const std::string &rhs, int overwrite, bool compress, int storageType) :
    _isBigTiff(false),
    _storageType(storageType),
    _compress(compress)
{
    outputFilePath = rhs;
    
    if (overwrite == 0) {
        std::ifstream input_test;
        input_test.open(outputFilePath.c_str(), std::ios::in | std::ios::binary);
        if (input_test.good() == 1) {
            std::string error("The output file at ");
            error += this->outputFilePath;
            error += " already exists.";
            throw OUTPUT_FILE_ALREADY_EXISTS(error);	// escape without overwriting
        }
    }
    
    file.open(outputFilePath.c_str(), std::ios::binary | std::ios::out | std::ios::trunc);	// DANGER: OVERWRITING THE FILE
    if (file.fail() != 0) {
        std::string error("Cannot create an output file at ");
        error += this->outputFilePath;
        throw CANNOT_OPEN_OUTPUT_FILE(error);
    }
    
    _writeTiffHeader();
}

LocalizerTIFFImageOutputWriter::~LocalizerTIFFImageOutputWriter() {
    if (file.is_open()) {
        _touchupOffsets();
        if (file.fail())
            throw std::runtime_error("Error trying to write to the TIFF file");
    }
}

void LocalizerTIFFImageOutputWriter::write_image(ImagePtr image) {
    if ((file.tellp() % (uint64_t)2) != 0) {
        // ensure file starts on word boundary
        WriteBinaryValue<uint8_t>(file, 0);
    }
    
    uint64_t currentOffset = file.tellp();
    std::pair<LocalizerTIFFImageOutputWriter::TIFFIFDOnDisk, std::vector<char> > newIFDParams;
    newIFDParams = _constructIFD(image, currentOffset, _isBigTiff);
    const std::vector<char>& ifdData = newIFDParams.second;
    
    // is this a regular tiff, and would the new IFD push us over the edge?
    if (!_isBigTiff && (currentOffset + (uint64_t)ifdData.size() > (uint64_t)std::numeric_limits<uint32_t>::max())) {
        _convertToBigTiff();
        write_image(image);
        return;
    }
    
    file.write(ifdData.data(), ifdData.size());
    if (file.fail())
        throw std::runtime_error("Error trying to write to the TIFF file");
    _writtenIFDs.push_back(newIFDParams.first);
}

std::pair<LocalizerTIFFImageOutputWriter::TIFFIFDOnDisk, std::vector<char> > LocalizerTIFFImageOutputWriter::_constructIFD(ImagePtr image, uint64_t ifdWillBeAtThisOffset, bool isBigTiff, bool reuseExistingData, LocalizerTIFFImageOutputWriter::TIFFIFDOnDisk existingIFDOnDisk) const {
    
    int nRows, nCols;
    if (reuseExistingData) {
        nRows = existingIFDOnDisk.nRows;
        nCols = existingIFDOnDisk.nCols;
    } else {
        nRows = image->rows();
        nCols = image->cols();
    }
    
    // determine how much storage we need.
    int sampleFormat, bitsPerSample;
    TIFFSampleFormatAndBitsPerSampleForFormat(_storageType, sampleFormat, bitsPerSample);
    uint64_t dataLength;
    std::vector<char> imageBuffer;
    if (!reuseExistingData) {
        if (!_compress) {
            ImageToBufferWithFormat(image, _storageType, imageBuffer);
        } else {
            std::vector<char> uncompressedData;
            ImageToBufferWithFormat(image, _storageType, uncompressedData);
            imageBuffer = Deflate(uncompressedData);
        }
        dataLength = imageBuffer.size();
    } else {
        dataLength = existingIFDOnDisk.dataLength;
    }
    
    int nTIFFTags = 9;
    uint64_t IFDLength = (isBigTiff) ? (8 + nTIFFTags * 20 + 8) : (2 + nTIFFTags * 12 + 4);
    uint64_t fullIFDLength = (reuseExistingData) ? IFDLength : (IFDLength + dataLength);
    uint64_t nextIFDFieldOffset = ifdWillBeAtThisOffset + IFDLength - ((isBigTiff) ? 8 : 4);
    std::vector<char> outputBuffer(fullIFDLength);
    char* bufferPtr = outputBuffer.data();
	uint64_t dataOffset;
	if (!reuseExistingData) {
		dataOffset = ifdWillBeAtThisOffset + IFDLength;
	} else {
		dataOffset = existingIFDOnDisk.dataOffset;
	}
    
    // write all of the TIFF tags
    // number of entries
    if (isBigTiff) {
        _storeInBuffer<uint64_t>(bufferPtr, nTIFFTags);
    } else {
        _storeInBuffer<uint16_t>(bufferPtr, nTIFFTags);
    }
    
    _writeTag(bufferPtr, TIFFTAG_IMAGEWIDTH, 1, nRows, isBigTiff);
    _writeTag(bufferPtr, TIFFTAG_IMAGELENGTH, 1, nCols, isBigTiff);
    _writeTag(bufferPtr, TIFFTAG_BITSPERSAMPLE, 1, bitsPerSample, isBigTiff);
    _writeTag(bufferPtr, TIFFTAG_COMPRESSION, 1, (_compress ? COMPRESSION_DEFLATE : COMPRESSION_NONE), isBigTiff);
    _writeTag(bufferPtr, TIFFTAG_PHOTOMETRIC, 1, PHOTOMETRIC_MINISBLACK, isBigTiff);
    _writeTag(bufferPtr, TIFFTAG_STRIPOFFSETS, 1, dataOffset, isBigTiff);
    _writeTag(bufferPtr, TIFFTAG_ROWSPERSTRIP, 1, nRows * nCols, isBigTiff);
    _writeTag(bufferPtr, TIFFTAG_STRIPBYTECOUNTS, 1, dataLength, isBigTiff);
    _writeTag(bufferPtr, TIFFTAG_SAMPLEFORMAT, 1, sampleFormat, isBigTiff);
    
    // offset to next IFD (0 for now)
    if (isBigTiff) {
        _storeInBuffer<uint64_t>(bufferPtr, 0);
    } else {
        _storeInBuffer<uint32_t>(bufferPtr, 0);
    }
    
    // safety check
    if ((uint64_t)(bufferPtr - outputBuffer.data()) != IFDLength) {
        throw std::logic_error("buffer length mismatch");
    }
    
    // store actual image data if needed
    if (!reuseExistingData) {
        memcpy(bufferPtr, imageBuffer.data(),imageBuffer.size());
        bufferPtr += imageBuffer.size();
    }
    
    // safety check
    if ((uint64_t)(bufferPtr - outputBuffer.data()) != fullIFDLength) {
        throw std::logic_error("buffer full length mismatch");
    }
    
    // we should be all done. compile the necessary info to piece this together
    LocalizerTIFFImageOutputWriter::TIFFIFDOnDisk ifdOnDisk;
    ifdOnDisk.ifdOffset = ifdWillBeAtThisOffset;
    ifdOnDisk.nRows = nRows;
    ifdOnDisk.nCols = nCols;
    ifdOnDisk.dataOffset = dataOffset;
    ifdOnDisk.dataLength = dataLength;
    ifdOnDisk.nextIFDFieldOffset = nextIFDFieldOffset;
    
    // safety check
    if (reuseExistingData && (ifdOnDisk.dataOffset != existingIFDOnDisk.dataOffset)) {
        throw std::logic_error("existing offset and new offset mismatch");
    }
    
    return std::pair<LocalizerTIFFImageOutputWriter::TIFFIFDOnDisk, std::vector<char> >(ifdOnDisk, outputBuffer);
}

// all tags are written with TIFF_LONG as the type for a regular TIFF and TIFF_LONG8 for BigTIFF
void LocalizerTIFFImageOutputWriter::_writeTag(char*& bufferPtr, int tagID, uint64_t count, uint64_t value, bool isBigTiff) const {
    _storeInBuffer<uint16_t>(bufferPtr, tagID);
    _storeInBuffer<uint16_t>(bufferPtr, (isBigTiff) ? TIFF_LONG8 : TIFF_LONG);
    if (isBigTiff) {
        _storeInBuffer<uint64_t>(bufferPtr, count);
        _storeInBuffer<uint64_t>(bufferPtr, value);
    } else {
        _storeInBuffer<uint32_t>(bufferPtr, count);
        _storeInBuffer<uint32_t>(bufferPtr, value);
    }
    
}

void LocalizerTIFFImageOutputWriter::_convertToBigTiff() {
    if (_isBigTiff)
        throw std::logic_error("converting tiff that is already bigtiff");
    _isBigTiff = true;
    
    _writeBigTiffHeader();
    
    // we need to loop over all of the IFDs that have been written thus far, and write new BigTIFF IFDs for all of them.
    // we don't change the existing file contents - we simply write a bunch of new ones.
    std::vector<LocalizerTIFFImageOutputWriter::TIFFIFDOnDisk> newIFDs;
    file.seekp(0, std::ios_base::end);
    for (auto it = _writtenIFDs.cbegin(); it != _writtenIFDs.cend(); ++it) {
        if ((file.tellp() % 2) != 0) {
            // enforce word boundary
            WriteBinaryValue<uint8_t>(file, 0);
        }
        std::pair<LocalizerTIFFImageOutputWriter::TIFFIFDOnDisk, std::vector<char> > newIFDParams;
        newIFDParams = _constructIFD(ImagePtr(), file.tellp(), _isBigTiff, true, *it);
        file.write(newIFDParams.second.data(), newIFDParams.second.size());
        newIFDs.push_back(newIFDParams.first);
    }
    
    _writtenIFDs = newIFDs;
}

void LocalizerTIFFImageOutputWriter::_writeTiffHeader() {
    file.seekp(0);
    WriteBinaryValue<uint16_t>(file, 0x4949);
    WriteBinaryValue<uint16_t>(file, 0x002A);
    WriteBinaryValue<uint32_t>(file, 0);
    WriteBinaryValue<uint64_t>(file, 0);
    if (file.fail())
        throw std::runtime_error("Error trying to write to the TIFF file");
}

void LocalizerTIFFImageOutputWriter::_writeBigTiffHeader() {
    file.seekp(0);
    WriteBinaryValue<uint16_t>(file, 0x4949);
    WriteBinaryValue<uint16_t>(file, 0x002B);
    WriteBinaryValue<uint16_t>(file, 0x0008);
    WriteBinaryValue<uint16_t>(file, 0);
    WriteBinaryValue<uint64_t>(file, 0);
    if (file.fail())
        throw std::runtime_error("Error trying to write to the TIFF file");
}

void LocalizerTIFFImageOutputWriter::_touchupOffsets() {
    if (!_writtenIFDs.empty()) {
        // first IFD is special since it needs to be referred to in the header
        uint64_t firstIFDOffset = _writtenIFDs.at(0).ifdOffset;
        if (_isBigTiff) {
            file.seekp(8);
            WriteBinaryValue<uint64_t>(file, firstIFDOffset);
        } else {
            file.seekp(4);
            WriteBinaryValue<uint32_t>(file, firstIFDOffset);
        }
    }
    
    for (size_t i = 1; i < _writtenIFDs.size(); ++i) {
        uint64_t previousIFDOffsetOffset = _writtenIFDs.at(i - 1).nextIFDFieldOffset;
        uint64_t thisIFDOffset = _writtenIFDs.at(i).ifdOffset;
        file.seekp(previousIFDOffsetOffset);
        if (_isBigTiff) {
            WriteBinaryValue<uint64_t>(file, thisIFDOffset);
        } else {
            WriteBinaryValue<uint32_t>(file, thisIFDOffset);
        }
    }
    
    if (file.fail())
        throw std::runtime_error("Error trying to write to the TIFF file");
    
    file.seekp(0, std::ios_base::end);
    
}

#ifdef WITH_IGOR
IgorImageOutputWriter::IgorImageOutputWriter(std::string waveName_rhs, size_t nImages_rhs, int overwrite_rhs, int storageType_rhs) {
	
	this->overwrite = overwrite_rhs;
	this->outputWave = NULL;
	this->fullPathToWave = waveName_rhs;
	this->nImagesTotal = nImages_rhs;
	this->nImagesWritten = 0;
    
    // not all storage types are supported in Igor
    // silently replace these with supported types
    switch (storageType_rhs) {
        case STORAGE_TYPE_INT4:
        case STORAGE_TYPE_UINT4:
            this->storageType = STORAGE_TYPE_INT8;
            break;
        case STORAGE_TYPE_INT64:
            this->storageType = STORAGE_TYPE_INT32;
            break;
        case STORAGE_TYPE_UINT64:
            this->storageType = STORAGE_TYPE_UINT32;
            break;
        default:
            this->storageType = storageType_rhs;
    }
	
	// check if the output wave already exists
	if (this->overwrite != 1) {
		waveHndl wave = FetchWaveUsingFullPath(waveName_rhs);
		if (wave != NULL) {
			std::string error ("the output wave ");
			error += waveName_rhs;
			error += " already exists (and the /O flag was not specified)";
		}
	}
}

IgorImageOutputWriter::IgorImageOutputWriter(DataFolderAndName outputDataFolderAndName_rhs, size_t nImages_rhs, int overwrite_rhs, int storageType_rhs) {
	
	this->overwrite = overwrite_rhs;
	this->outputWave = NULL;
	this->fullPathToWave.assign("");
	this->nImagesTotal = nImages_rhs;
	this->nImagesWritten = 0;
	this->storageType = storageType_rhs;
	this->waveDataFolderAndName = outputDataFolderAndName_rhs;
	
	// check if the output wave already exists
	if (this->overwrite != 1) {
		waveHndl wave = FetchWaveFromDataFolder(this->waveDataFolderAndName.dfH, this->waveDataFolderAndName.name);
		if (wave != NULL) {
			throw std::runtime_error("the requested output wave already exists (and the /O flag was not specified)");
		}
	}
}

void IgorImageOutputWriter::write_image(ImagePtr imageToWrite) {
	
	int result;
	size_t xSize = imageToWrite->rows();
	size_t ySize = imageToWrite->cols();
	size_t nPixels = xSize * ySize;
	int storage = this->GetIgorStorageType();
	
	if (this->outputWave == NULL) {
		// the outputwave has not been created yet, do it now
		CountInt dimensionSizes[MAX_DIMENSIONS + 1];
		
		dimensionSizes[0] = xSize;
		dimensionSizes[1] = ySize;
        if (this->nImagesTotal == 1) {
            dimensionSizes[2] = 0;  // some Igor operations insist on an mxn wave and complain about mxnx1
        } else {
            dimensionSizes[2] = this->nImagesTotal;
        }
		dimensionSizes[3] = 0;
		
		// the way to make the wave depends on whether this object was constructed with a full path
		// or with a DataFolderAndName argument
		if (this->fullPathToWave.length() != 0) {
			this->outputWave = MakeWaveUsingFullPath(this->fullPathToWave, dimensionSizes, storage, this->overwrite);
		} else {
			result = MDMakeWave(&(this->outputWave), this->waveDataFolderAndName.name, this->waveDataFolderAndName.dfH, dimensionSizes, storage, this->overwrite);
			if (result != 0)
				throw result;
		}
	}
	
	// check that we are not trying to write too many images
	// which would otherwise trigger an out-of-bounds memory
	// access
    if (this->nImagesWritten >= this->nImagesTotal)
		throw std::runtime_error("Writing too many images to the IgorImageOutputWriter");
	
	// the strategy for writing the data depends on the storage type
	size_t waveDataOffset;
	char *waveDataPtr;
	
	result = MDAccessNumericWaveData(this->outputWave, kMDWaveAccessMode0, (BCInt*)&waveDataOffset);
	if (result != 0)
		throw result;
	
	waveDataPtr = (char *)((char*)(*this->outputWave) + waveDataOffset);
	double* imagePtr = imageToWrite->data();
	
	switch (storage) {
		case NT_I8:
		{
			int8_t* int8Ptr = (int8_t*)waveDataPtr;
			int8Ptr += nPixels * this->nImagesWritten;
			for (size_t i = 0; i < nPixels; ++i) {
				int8Ptr[i] = imagePtr[i];
			}
			break;
		}
		case NT_I16:
		{
			int16_t* int16Ptr = (int16_t*)waveDataPtr;
			int16Ptr += nPixels * this->nImagesWritten;
			for (size_t i = 0; i < nPixels; ++i) {
				int16Ptr[i] = imagePtr[i];
			}
			break;
		}
		case NT_I32:
		{
			int32_t* int32Ptr = (int32_t*)waveDataPtr;
			int32Ptr += nPixels * this->nImagesWritten;
			for (size_t i = 0; i < nPixels; ++i) {
				int32Ptr[i] = imagePtr[i];
			}
			break;
		}
		case NT_I8 | NT_UNSIGNED:
		{
			uint8_t* uint8Ptr = (uint8_t*)waveDataPtr;
			uint8Ptr += nPixels * this->nImagesWritten;
			for (size_t i = 0; i < nPixels; ++i) {
				uint8Ptr[i] = imagePtr[i];
			}
			break;
		}
		case NT_I16 | NT_UNSIGNED:
		{
			uint16_t* uint16Ptr = (uint16_t*)waveDataPtr;
			uint16Ptr += nPixels * this->nImagesWritten;
			for (size_t i = 0; i < nPixels; ++i) {
				uint16Ptr[i] = imagePtr[i];
			}
			break;
		}
		case NT_I32 | NT_UNSIGNED:
		{
			uint32_t* uint32Ptr = (uint32_t*)waveDataPtr;
			uint32Ptr += nPixels * this->nImagesWritten;
			for (size_t i = 0; i < nPixels; ++i) {
				uint32Ptr[i] = imagePtr[i];
			}
			break;
		}
		case NT_FP32:
		{
			float* floatPtr = (float*)waveDataPtr;
			floatPtr += nPixels * this->nImagesWritten;
			for (size_t i = 0; i < nPixels; ++i) {
				floatPtr[i] = imagePtr[i];
			}
			break;
		}
		case NT_FP64:
		{
			double* doublePtr = (double*)waveDataPtr;
			doublePtr += nPixels * this->nImagesWritten;
			memcpy(doublePtr, imagePtr, nPixels * sizeof(double));
			break;
		}
	}
	
	++nImagesWritten;
}

int IgorImageOutputWriter::GetIgorStorageType() {
	int storage;
	
	switch (this->storageType) {
		case STORAGE_TYPE_INT4:
		case STORAGE_TYPE_UINT4:
		case STORAGE_TYPE_INT8:
			storage = NT_I8;
			break;
		case STORAGE_TYPE_UINT8:
			storage = NT_I8 | NT_UNSIGNED;
			break;
		case STORAGE_TYPE_INT16:
			storage = NT_I16;
			break;
		case STORAGE_TYPE_UINT16:
			storage = NT_I16 | NT_UNSIGNED;
			break;
		case STORAGE_TYPE_INT32:
			storage = NT_I32;
			break;
		case STORAGE_TYPE_UINT32:
			storage = NT_I32 | NT_UNSIGNED;
			break;
        case STORAGE_TYPE_INT64:
            storage = NT_I32;
			break;
		case STORAGE_TYPE_UINT64:
			storage = NT_I32 | NT_UNSIGNED;
			break;
		case STORAGE_TYPE_FP32:
			storage = NT_FP32;
			break;
		case STORAGE_TYPE_FP64:
			storage = NT_FP64;
			break;
	}
	
	return storage;
}

#endif // WITH_IGOR
