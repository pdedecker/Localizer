#include "PALM_analysis.h"


boost::shared_ptr<PALMMatrix <unsigned char> > do_processing_and_thresholding(boost::shared_ptr<PALMMatrix<double> > image, boost::shared_ptr<ThresholdImage_Preprocessor>preprocessor, 
																		 boost::shared_ptr<ThresholdImage> thresholder, boost::shared_ptr<ThresholdImage_Postprocessor> postprocessor) {
	// this function takes care of the thresholding and the associated pre- and postprocessing
	
	boost::shared_ptr<PALMMatrix<double> > preprocessed_image;
	boost::shared_ptr<PALMMatrix <unsigned char> > thresholded_image;
	boost::shared_ptr<PALMMatrix <unsigned char> > postprocessed_image;
	
	if (preprocessor.get() != NULL) {
		preprocessed_image = preprocessor->do_preprocessing(image);
		thresholded_image = thresholder->do_thresholding(preprocessed_image);
		
	} else {
		thresholded_image = thresholder->do_thresholding(image);
	}
	
	if (postprocessor.get() != NULL) {
		postprocessed_image = postprocessor->do_postprocessing(thresholded_image, image);
		thresholded_image = postprocessed_image;
	}
	
	return thresholded_image;
}

#ifdef WITH_IGOR
boost::shared_ptr<LocalizedPositionsContainer> LocalizedPositionsContainer::GetPositionsFromWave(waveHndl positionsWave) {
	int err;
	size_t findPosition;
	
	// determine the type of positions being passed
	Handle waveNoteHandle = WaveNote(positionsWave);
	size_t waveNoteSize = GetHandleSize(positionsWave);
	boost::scoped_array<char> CStringWaveNote(new char[waveNoteSize + 1]);
	
	err = GetCStringFromHandle(waveNoteHandle, CStringWaveNote.get(), waveNoteSize);
	if (err != 0)
		throw std::runtime_error("GetCStringFromHandle() returned a nonzero code");
	
	// save the wavenote as a std::string
	std::string waveNote(CStringWaveNote.get());
	
	// see if the wave note contains info on the kind of localization used
	// if not then fail
	
	findPosition = waveNote.find("LOCALIZATION METHOD:");
	if (findPosition == (size_t)-1)	// not found
		throw std::runtime_error("The positions wave does not specify a localization method");
	
	findPosition = waveNote.find("LOCALIZATION METHOD:0;");
	if (findPosition != (size_t)-1) {
		return boost::shared_ptr<LocalizedPositionsContainer> (new LocalizedPositionsContainer_2DGauss(positionsWave));
	}
	
	findPosition = waveNote.find("LOCALIZATION METHOD:1;");
	if (findPosition != (size_t)-1) {
		return boost::shared_ptr<LocalizedPositionsContainer> (new LocalizedPositionsContainer_2DGaussFixedWidth(positionsWave));
	}
	
	findPosition = waveNote.find("LOCALIZATION METHOD:2");
	if (findPosition != (size_t)-1) {
		return boost::shared_ptr<LocalizedPositionsContainer> (new LocalizedPositionsContainer_Multiplication(positionsWave));
	}
	
	findPosition = waveNote.find("LOCALIZATION METHOD:3");
	if (findPosition != (size_t)-1) {
		return boost::shared_ptr<LocalizedPositionsContainer> (new LocalizedPositionsContainer_Centroid(positionsWave));
	}
	
	findPosition = waveNote.find("LOCALIZATION METHOD:4");
	if (findPosition != (size_t)-1) {
		return boost::shared_ptr<LocalizedPositionsContainer> (new LocalizedPositionsContainer_ZeissPALM(positionsWave));
	}
	
	// if we are still here then we don't recognize the type of localization used
	throw std::runtime_error("Unknown localization method (check the wave note of the wave containing the positions)");
}
#endif // #ifdef WITH_IGOR

void LocalizedPositionsContainer_2DGauss::addPosition(boost::shared_ptr<LocalizedPosition> newPosition) {
	// check if the type of positions that we are adding is suitable
	if (newPosition->getPositionType() != LOCALIZED_POSITIONS_TYPE_2DGAUSS)
		throw std::runtime_error("Trying to append a position of a different type to a LocalizedPositionsContainer_2DGauss");
	
	// cast the pointer to the more specific type
	boost::shared_ptr<LocalizedPosition_2DGauss> newPosition_2DGauss(boost::static_pointer_cast<LocalizedPosition_2DGauss> (newPosition));
	
	this->positionsVector.push_back(*newPosition_2DGauss);
}

void LocalizedPositionsContainer_2DGauss::addPositions(boost::shared_ptr<LocalizedPositionsContainer> newPositionsContainer) {
	// are we trying to add the same container to itself?
	if (this == newPositionsContainer.get()) {
		throw std::runtime_error("Trying to append a LocalizedPositionsContainer_2DGauss to itself");
	}
	
	// check if the positions container is of the right type
	if (newPositionsContainer->getPositionsType() != LOCALIZED_POSITIONS_TYPE_2DGAUSS) {
		throw std::runtime_error("Trying to append a position of a different type to a LocalizedPositionsContainer_2DGauss");
	}
	
	// cast the pointer to the more specific type
	boost::shared_ptr<LocalizedPositionsContainer_2DGauss> newPositionsContainer_2DGauss(boost::static_pointer_cast<LocalizedPositionsContainer_2DGauss> (newPositionsContainer));
	
	for (std::vector<LocalizedPosition_2DGauss>::iterator it = newPositionsContainer_2DGauss->positionsVector.begin(); it != newPositionsContainer_2DGauss->positionsVector.end(); ++it) {
		this->positionsVector.push_back(*it);
	}
}

#ifdef WITH_IGOR
LocalizedPositionsContainer_2DGauss::LocalizedPositionsContainer_2DGauss(waveHndl positionsWave) {
	// initialize a new PositionsContainer from a wave that contains positions of the correct type
	long numDimensions;
	long dimensionSizes[MAX_DIMENSIONS+1];
	int err;
	
	err = MDGetWaveDimensions(positionsWave, &numDimensions, dimensionSizes);
	
	if ((numDimensions != 2) || (dimensionSizes[1] != 12)) {	// invalid dimensions (warning: magic numbers)
		throw (std::runtime_error("Invalid positions wave"));
	}
	
	LocalizedPosition_2DGauss singlePosition;
	size_t nPositions = dimensionSizes[0];
	this->positionsVector.reserve(nPositions);
	long indices[MAX_DIMENSIONS];
	double value[2];
	
	for (size_t i = 0; i < nPositions; ++i) {
		// get all the relevant data out of the wave and into a position object
		indices[0] = i;
		indices[1] = 0;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.frameNumber = value[0];
		indices[1] = 1;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.integral = value[0];
		indices[1] = 2;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.width = value[0];
		indices[1] = 3;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.xPosition = value[0];
		indices[1] = 4;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.yPosition = value[0];
		indices[1] = 5;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.background = value[0];
		indices[1] = 6;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.integralDeviation = value[0];
		indices[1] = 7;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.widthDeviation = value[0];
		indices[1] = 8;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.xPositionDeviation = value[0];
		indices[1] = 9;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.yPositionDeviation = value[0];
		indices[1] = 10;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.backgroundDeviation = value[0];
		indices[1] = 11;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.nFramesPresent = value[0];
		
		this->positionsVector.push_back(singlePosition);
	}
}

waveHndl LocalizedPositionsContainer_2DGauss::writePositionsToWave(std::string waveName, std::string waveNote) const {
	long dimensionSizes[MAX_DIMENSIONS+1];
	int err;
	waveHndl outputWave;
	size_t nPositions = this->positionsVector.size();
	dimensionSizes[0] = nPositions;
	dimensionSizes[1] = 12;	// magic number
	dimensionSizes[2] = 0;
	err = MDMakeWave(&outputWave, waveName.c_str(), NULL, dimensionSizes, NT_FP64, 1);
	if (err != 0)
		throw err;
	
	long indices[MAX_DIMENSIONS];
	double value[2];
	
	for (size_t i = 0; i < nPositions; ++i) {
		indices[0] = i;
		indices[1] = 0;
		value[0] = this->positionsVector.at(i).frameNumber;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 1;
		value[0] = this->positionsVector.at(i).integral;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 2;
		value[0] = this->positionsVector.at(i).width;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 3;
		value[0] = this->positionsVector.at(i).xPosition;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 4;
		value[0] = this->positionsVector.at(i).yPosition;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 5;
		value[0] = this->positionsVector.at(i).background;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 6;
		value[0] = this->positionsVector.at(i).integralDeviation;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 7;
		value[0] = this->positionsVector.at(i).widthDeviation;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 8;
		value[0] = this->positionsVector.at(i).xPositionDeviation;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 9;
		value[0] = this->positionsVector.at(i).yPositionDeviation;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 10;
		value[0] = this->positionsVector.at(i).backgroundDeviation;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 11;
		value[0] = this->positionsVector.at(i).nFramesPresent;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
	}
	
	if (waveNote.size() != 0) {
		// set the waveNote to the string passed in
		Handle waveNoteHandle = NewHandle(waveNote.length());
		if (waveNoteHandle == NULL)
			throw std::bad_alloc();
		
		PutCStringInHandle(waveNote.c_str(), waveNoteHandle);
		SetWaveNote(outputWave, waveNoteHandle);
	}
	
	return outputWave;
}
#endif // #ifdef WITH_IGOR

void LocalizedPositionsContainer_2DGauss::writePositionsToFile(std::string filePath, std::string header) const {
	boost::filesystem::ofstream outputFile;
	size_t nPositions = this->positionsVector.size();
	if (nPositions == 0) {
		throw std::runtime_error("No positions localized!");
	}
	
	outputFile.open(filePath, ios::trunc);
	if (outputFile.fail()) {
		throw std::runtime_error("Unable to create the output file");
	}
	
	// write the header. Add an extra newline to be safe.
	outputFile << "Localized positions using symmetric 2D Gauss fitting" << "\n";
	outputFile << header << "\n";
	outputFile << "First frame" << "\t";
	outputFile << "Integrated intensity" << "\t";
	outputFile << "Fitted PSF standard deviation" << "\t";
	outputFile << "X position (pixel)" << "\t";
	outputFile << "Y position (pixel)" << "\t";
	outputFile << "Background" << "Intensity deviation" << "\t";
	outputFile << "PSF deviation" << "\t";
	outputFile << "X position deviation" << "\t";
	outputFile << "Y position deviation" << "\t";
	outputFile << "Background deviation" << "\t";
	outputFile << "Number of frames where this emitter is present" << "\n";
	
	// write the actual positions
	for (std::vector<LocalizedPosition_2DGauss>::const_iterator it = this->positionsVector.begin(); it != this->positionsVector.end(); ++it) {
		outputFile << (*it).frameNumber << "\t";
		outputFile << (*it).integral << "\t";
		outputFile << (*it).width << "\t";
		outputFile << (*it).xPosition << "\t";
		outputFile << (*it).yPosition << "\t";
		outputFile << (*it).background << "\t";
		outputFile << (*it).integralDeviation << "\t";
		outputFile << (*it).widthDeviation << "\t";
		outputFile << (*it).xPositionDeviation << "\t";
		outputFile << (*it).yPositionDeviation << "\t";
		outputFile << (*it).backgroundDeviation << "\t";
		outputFile << (*it).nFramesPresent << "\t";
		
		if (outputFile.fail()) {
			if (outputFile.is_open()) {
				outputFile.close();
			}
			throw std::runtime_error("Error writing to the localized positions file");
		}
	}
	
	// write a terminating newline
	outputFile << std::endl;
	outputFile.close();
}

#ifdef WITH_IGOR
LocalizedPositionsContainer_2DGaussFixedWidth::LocalizedPositionsContainer_2DGaussFixedWidth(waveHndl positionsWave) {
	// initialize a new PositionsContainer from a wave that contains positions of the correct type
	long numDimensions;
	long dimensionSizes[MAX_DIMENSIONS+1];
	int err;
	
	err = MDGetWaveDimensions(positionsWave, &numDimensions, dimensionSizes);
	
	if ((numDimensions != 2) || (dimensionSizes[1] != 11)) {	// invalid dimensions (warning: magic numbers)
		throw (std::runtime_error("Invalid positions wave"));
	}
	
	LocalizedPosition_2DGaussFixedWidth singlePosition;
	size_t nPositions = dimensionSizes[0];
	this->positionsVector.reserve(nPositions);
	long indices[MAX_DIMENSIONS];
	double value[2];
	
	for (size_t i = 0; i < nPositions; ++i) {
		// get all the relevant data out of the wave and into a position object
		indices[0] = i;
		indices[1] = 0;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.frameNumber = value[0];
		indices[1] = 1;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.integral = value[0];
		indices[1] = 2;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.width = value[0];
		indices[1] = 3;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.xPosition = value[0];
		indices[1] = 4;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.yPosition = value[0];
		indices[1] = 5;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.background = value[0];
		indices[1] = 6;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.integralDeviation = value[0];
		indices[1] = 7;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.xPositionDeviation = value[0];
		indices[1] = 8;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.yPositionDeviation = value[0];
		indices[1] = 9;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.backgroundDeviation = value[0];
		indices[1] = 10;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.nFramesPresent = value[0];
		
		this->positionsVector.push_back(singlePosition);
	}
}

waveHndl LocalizedPositionsContainer_2DGaussFixedWidth::writePositionsToWave(std::string waveName, std::string waveNote) const {
	long dimensionSizes[MAX_DIMENSIONS+1];
	int err;
	waveHndl outputWave;
	size_t nPositions = this->positionsVector.size();
	dimensionSizes[0] = nPositions;
	dimensionSizes[1] = 11;	// magic number
	dimensionSizes[2] = 0;
	err = MDMakeWave(&outputWave, waveName.c_str(), NULL, dimensionSizes, NT_FP64, 1);
	if (err != 0)
		throw err;
	
	long indices[MAX_DIMENSIONS];
	double value[2];
	
	for (size_t i = 0; i < nPositions; ++i) {
		indices[0] = i;
		indices[1] = 0;
		value[0] = this->positionsVector.at(i).frameNumber;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 1;
		value[0] = this->positionsVector.at(i).integral;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 2;
		value[0] = this->positionsVector.at(i).width;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 3;
		value[0] = this->positionsVector.at(i).xPosition;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 4;
		value[0] = this->positionsVector.at(i).yPosition;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 5;
		value[0] = this->positionsVector.at(i).background;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 6;
		value[0] = this->positionsVector.at(i).integralDeviation;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 7;
		value[0] = this->positionsVector.at(i).xPositionDeviation;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 8;
		value[0] = this->positionsVector.at(i).yPositionDeviation;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 9;
		value[0] = this->positionsVector.at(i).backgroundDeviation;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 10;
		value[0] = this->positionsVector.at(i).nFramesPresent;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
	}
	
	if (waveNote.size() != 0) {
		// set the waveNote to the string passed in
		Handle waveNoteHandle = NewHandle(waveNote.length());
		if (waveNoteHandle == NULL)
			throw std::bad_alloc();
		
		PutCStringInHandle(waveNote.c_str(), waveNoteHandle);
		SetWaveNote(outputWave, waveNoteHandle);
	}
	
	return outputWave;
}
#endif // #ifdef WITH_IGOR

void LocalizedPositionsContainer_2DGaussFixedWidth::writePositionsToFile(std::string filePath, std::string header) const {
	boost::filesystem::ofstream outputFile;
	size_t nPositions = this->positionsVector.size();
	if (nPositions == 0) {
		throw std::runtime_error("No positions localized!");
	}
	
	outputFile.open(filePath, ios::trunc);
	if (outputFile.fail()) {
		throw std::runtime_error("Unable to create the output file");
	}
	
	// write the header. Add an extra newline to be safe.
	outputFile << "Localized positions using symmetric 2D Gauss fitting with fixed PSF width" << "\n";
	outputFile << header << "\n";
	outputFile << "First frame" << "\t";
	outputFile << "Integrated intensity" << "\t";
	outputFile << "Fitted PSF standard deviation" << "\t";
	outputFile << "X position (pixel)" << "\t";
	outputFile << "Y position (pixel)" << "\t";
	outputFile << "Background" << "Intensity deviation" << "\t";
	outputFile << "X position deviation" << "\t";
	outputFile << "Y position deviation" << "\t";
	outputFile << "Background deviation" << "\t";
	outputFile << "Number of frames where this emitter is present" << "\n";
	
	// write the actual positions
	for (std::vector<LocalizedPosition_2DGaussFixedWidth>::const_iterator it = this->positionsVector.begin(); it != this->positionsVector.end(); ++it) {
		outputFile << (*it).frameNumber << "\t";
		outputFile << (*it).integral << "\t";
		outputFile << (*it).width << "\t";
		outputFile << (*it).xPosition << "\t";
		outputFile << (*it).yPosition << "\t";
		outputFile << (*it).background << "\t";
		outputFile << (*it).integralDeviation << "\t";
		outputFile << (*it).xPositionDeviation << "\t";
		outputFile << (*it).yPositionDeviation << "\t";
		outputFile << (*it).backgroundDeviation << "\t";
		outputFile << (*it).nFramesPresent << "\t";
		
		if (outputFile.fail()) {
			if (outputFile.is_open()) {
				outputFile.close();
			}
			throw std::runtime_error("Error writing to the localized positions file");
		}
	}
	
	// write a terminating newline
	outputFile << std::endl;
	outputFile.close();
}

void LocalizedPositionsContainer_2DGaussFixedWidth::addPosition(boost::shared_ptr<LocalizedPosition> newPosition) {
	// check if the type of positions that we are adding is suitable
	if (newPosition->getPositionType() != LOCALIZED_POSITIONS_TYPE_2DGAUSS_FIXED_WIDTH)
		throw std::runtime_error("Trying to append a position of a different type to a LocalizedPositionsContainer_2DGauss");
	
	// cast the pointer to the more specific type
	boost::shared_ptr<LocalizedPosition_2DGaussFixedWidth> newPosition_2DGaussFixedWidth(boost::static_pointer_cast<LocalizedPosition_2DGaussFixedWidth> (newPosition));
	
	this->positionsVector.push_back(*newPosition_2DGaussFixedWidth);
}

void LocalizedPositionsContainer_2DGaussFixedWidth::addPositions(boost::shared_ptr<LocalizedPositionsContainer> newPositionsContainer) {
	// are we trying to add the same container to itself?
	if (this == newPositionsContainer.get()) {
		throw std::runtime_error("Trying to append a LocalizedPositionsContainer_2DGaussFixedWidth to itself");
	}
	
	// check if the positions container is of the right type
	if (newPositionsContainer->getPositionsType() != LOCALIZED_POSITIONS_TYPE_2DGAUSS_FIXED_WIDTH) {
		throw std::runtime_error("Trying to append a position of a different type to a LocalizedPositionsContainer_2DGauss");
	}
	
	// cast the pointer to the more specific type
	boost::shared_ptr<LocalizedPositionsContainer_2DGaussFixedWidth> newPositionsContainer_2DGauss(boost::static_pointer_cast<LocalizedPositionsContainer_2DGaussFixedWidth> (newPositionsContainer));
	
	for (std::vector<LocalizedPosition_2DGaussFixedWidth>::iterator it = newPositionsContainer_2DGauss->positionsVector.begin(); it != newPositionsContainer_2DGauss->positionsVector.end(); ++it) {
		this->positionsVector.push_back(*it);
	}
}

#ifdef WITH_IGOR
LocalizedPositionsContainer_Centroid::LocalizedPositionsContainer_Centroid(waveHndl positionsWave) {
	// initialize a new PositionsContainer from a wave that contains positions of the correct type
	long numDimensions;
	long dimensionSizes[MAX_DIMENSIONS+1];
	int err;
	
	err = MDGetWaveDimensions(positionsWave, &numDimensions, dimensionSizes);
	
	if ((numDimensions != 2) || (dimensionSizes[1] != 4)) {	// invalid dimensions (warning: magic numbers)
		throw (std::runtime_error("Invalid positions wave"));
	}
	
	LocalizedPosition_Centroid singlePosition;
	size_t nPositions = dimensionSizes[0];
	this->positionsVector.reserve(nPositions);
	long indices[MAX_DIMENSIONS];
	double value[2];
	
	for (size_t i = 0; i < nPositions; ++i) {
		// get all the relevant data out of the wave and into a position object
		indices[0] = i;
		indices[1] = 0;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.frameNumber = value[0];
		indices[1] = 1;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.xPosition = value[0];
		indices[1] = 2;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.yPosition = value[0];
		indices[1] = 3;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.nFramesPresent = value[0];
		
		this->positionsVector.push_back(singlePosition);
	}
}

waveHndl LocalizedPositionsContainer_Centroid::writePositionsToWave(std::string waveName, std::string waveNote) const {
	long dimensionSizes[MAX_DIMENSIONS+1];
	int err;
	waveHndl outputWave;
	size_t nPositions = this->positionsVector.size();
	dimensionSizes[0] = nPositions;
	dimensionSizes[1] = 4;	// magic number
	dimensionSizes[2] = 0;
	err = MDMakeWave(&outputWave, waveName.c_str(), NULL, dimensionSizes, NT_FP64, 1);
	if (err != 0)
		throw err;
	
	long indices[MAX_DIMENSIONS];
	double value[2];
	
	for (size_t i = 0; i < nPositions; ++i) {
		indices[0] = i;
		indices[1] = 0;
		value[0] = this->positionsVector.at(i).frameNumber;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 1;
		value[0] = this->positionsVector.at(i).xPosition;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 2;
		value[0] = this->positionsVector.at(i).yPosition;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 3;
		value[0] = this->positionsVector.at(i).nFramesPresent;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
	}
	
	if (waveNote.size() != 0) {
		// set the waveNote to the string passed in
		Handle waveNoteHandle = NewHandle(waveNote.length());
		if (waveNoteHandle == NULL)
			throw std::bad_alloc();
		
		PutCStringInHandle(waveNote.c_str(), waveNoteHandle);
		SetWaveNote(outputWave, waveNoteHandle);
	}
	
	return outputWave;
}
#endif // WITH_IGOR

void LocalizedPositionsContainer_Centroid::writePositionsToFile(std::string filePath, std::string header) const {
	boost::filesystem::ofstream outputFile;
	size_t nPositions = this->positionsVector.size();
	if (nPositions == 0) {
		throw std::runtime_error("No positions localized!");
	}
	
	outputFile.open(filePath, ios::trunc);
	if (outputFile.fail()) {
		throw std::runtime_error("Unable to create the output file");
	}
	
	// write the header. Add an extra newline to be safe.
	outputFile << "Localized positions using symmetric 2D Gauss fitting" << "\n";
	outputFile << header << "\n";
	outputFile << "First frame" << "\t";
	outputFile << "X position (pixel)" << "\t";
	outputFile << "Y position (pixel)" << "\t";
	outputFile << "Number of frames where this emitter is present" << "\n";
	
	// write the actual positions
	for (std::vector<LocalizedPosition_Centroid>::const_iterator it = this->positionsVector.begin(); it != this->positionsVector.end(); ++it) {
		outputFile << (*it).frameNumber << "\t";
		outputFile << (*it).xPosition << "\t";
		outputFile << (*it).yPosition << "\t";
		outputFile << (*it).nFramesPresent << "\t";
		
		if (outputFile.fail()) {
			if (outputFile.is_open()) {
				outputFile.close();
			}
			throw std::runtime_error("Error writing to the localized positions file");
		}
	}
	
	// write a terminating newline
	outputFile << std::endl;
	outputFile.close();
}

void LocalizedPositionsContainer_Centroid::addPosition(boost::shared_ptr<LocalizedPosition> newPosition) {
	// check if the type of positions that we are adding is suitable
	if (newPosition->getPositionType() != LOCALIZED_POSITIONS_TYPE_CENTROID)
		throw std::runtime_error("Trying to append a position of a different type to a LocalizedPositionsContainer_Centroid");
	
	// cast the pointer to the more specific type
	boost::shared_ptr<LocalizedPosition_Centroid> newPosition_Centroid(boost::static_pointer_cast<LocalizedPosition_Centroid> (newPosition));
	
	this->positionsVector.push_back(*newPosition_Centroid);
}

void LocalizedPositionsContainer_Centroid::addPositions(boost::shared_ptr<LocalizedPositionsContainer> newPositionsContainer) {
	// are we trying to add the same container to itself?
	if (this == newPositionsContainer.get()) {
		throw std::runtime_error("Trying to append a LocalizedPositionsContainer_Centroid to itself");
	}
	
	// check if the positions container is of the right type
	if (newPositionsContainer->getPositionsType() != LOCALIZED_POSITIONS_TYPE_CENTROID) {
		throw std::runtime_error("Trying to append a position of a different type to a LocalizedPositionsContainer_Centroid");
	}
	
	// cast the pointer to the more specific type
	boost::shared_ptr<LocalizedPositionsContainer_Centroid> newPositionsContainer_Centroid(boost::static_pointer_cast<LocalizedPositionsContainer_Centroid> (newPositionsContainer));
	
	for (std::vector<LocalizedPosition_Centroid>::iterator it = newPositionsContainer_Centroid->positionsVector.begin(); it != newPositionsContainer_Centroid->positionsVector.end(); ++it) {
		this->positionsVector.push_back(*it);
	}
}

#ifdef WITH_IGOR
LocalizedPositionsContainer_Multiplication::LocalizedPositionsContainer_Multiplication(waveHndl positionsWave) {
	// initialize a new PositionsContainer from a wave that contains positions of the correct type
	long numDimensions;
	long dimensionSizes[MAX_DIMENSIONS+1];
	int err;
	
	err = MDGetWaveDimensions(positionsWave, &numDimensions, dimensionSizes);
	
	if ((numDimensions != 2) || (dimensionSizes[1] != 5)) {	// invalid dimensions (warning: magic numbers)
		throw (std::runtime_error("Invalid positions wave"));
	}
	
	LocalizedPosition_Multiplication singlePosition;
	size_t nPositions = dimensionSizes[0];
	this->positionsVector.reserve(nPositions);
	long indices[MAX_DIMENSIONS];
	double value[2];
	
	for (size_t i = 0; i < nPositions; ++i) {
		// get all the relevant data out of the wave and into a position object
		indices[0] = i;
		indices[1] = 0;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.frameNumber = value[0];
		indices[1] = 1;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.width = value[0];
		indices[1] = 2;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.xPosition = value[0];
		indices[1] = 3;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.yPosition = value[0];
		indices[1] = 4;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.nFramesPresent = value[0];
		
		this->positionsVector.push_back(singlePosition);
	}
}

waveHndl LocalizedPositionsContainer_Multiplication::writePositionsToWave(std::string waveName, std::string waveNote) const {
	long dimensionSizes[MAX_DIMENSIONS+1];
	int err;
	waveHndl outputWave;
	size_t nPositions = this->positionsVector.size();
	dimensionSizes[0] = nPositions;
	dimensionSizes[1] = 5;	// magic number
	dimensionSizes[2] = 0;
	err = MDMakeWave(&outputWave, waveName.c_str(), NULL, dimensionSizes, NT_FP64, 1);
	if (err != 0)
		throw err;
	
	long indices[MAX_DIMENSIONS];
	double value[2];
	
	for (size_t i = 0; i < nPositions; ++i) {
		indices[0] = i;
		indices[1] = 0;
		value[0] = this->positionsVector.at(i).frameNumber;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 1;
		value[0] = this->positionsVector.at(i).width;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 2;
		value[0] = this->positionsVector.at(i).xPosition;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 3;
		value[0] = this->positionsVector.at(i).yPosition;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 4;
		value[0] = this->positionsVector.at(i).nFramesPresent;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
	}
	
	if (waveNote.size() != 0) {
		// set the waveNote to the string passed in
		Handle waveNoteHandle = NewHandle(waveNote.length());
		if (waveNoteHandle == NULL)
			throw std::bad_alloc();
		
		PutCStringInHandle(waveNote.c_str(), waveNoteHandle);
		SetWaveNote(outputWave, waveNoteHandle);
	}
	
	return outputWave;
}
#endif // WITH_IGOR

void LocalizedPositionsContainer_Multiplication::writePositionsToFile(std::string filePath, std::string header) const {
	boost::filesystem::ofstream outputFile;
	size_t nPositions = this->positionsVector.size();
	if (nPositions == 0) {
		throw std::runtime_error("No positions localized!");
	}
	
	outputFile.open(filePath, ios::trunc);
	if (outputFile.fail()) {
		throw std::runtime_error("Unable to create the output file");
	}
	
	// write the header. Add an extra newline to be safe.
	outputFile << "Localized positions using symmetric 2D Gauss fitting" << "\n";
	outputFile << header << "\n";
	outputFile << "First frame" << "\t";
	outputFile << "Used PSF standard deviation" << "\t";
	outputFile << "X position (pixel)" << "\t";
	outputFile << "Y position (pixel)" << "\t";
	outputFile << "Number of frames where this emitter is present" << "\n";
	
	// write the actual positions
	for (std::vector<LocalizedPosition_Multiplication>::const_iterator it = this->positionsVector.begin(); it != this->positionsVector.end(); ++it) {
		outputFile << (*it).frameNumber << "\t";
		outputFile << (*it).width << "\t";
		outputFile << (*it).xPosition << "\t";
		outputFile << (*it).yPosition << "\t";
		outputFile << (*it).nFramesPresent << "\t";
		
		if (outputFile.fail()) {
			if (outputFile.is_open()) {
				outputFile.close();
			}
			throw std::runtime_error("Error writing to the localized positions file");
		}
	}
	
	// write a terminating newline
	outputFile << std::endl;
	outputFile.close();
}

void LocalizedPositionsContainer_Multiplication::addPosition(boost::shared_ptr<LocalizedPosition> newPosition) {
	// check if the type of positions that we are adding is suitable
	if (newPosition->getPositionType() != LOCALIZED_POSITIONS_TYPE_CENTROID)
		throw std::runtime_error("Trying to append a position of a different type to a LocalizedPositionsContainer_Multiplication");
	
	// cast the pointer to the more specific type
	boost::shared_ptr<LocalizedPosition_Multiplication> newPosition_Multiplication(boost::static_pointer_cast<LocalizedPosition_Multiplication> (newPosition));
	
	this->positionsVector.push_back(*newPosition_Multiplication);
}

void LocalizedPositionsContainer_Multiplication::addPositions(boost::shared_ptr<LocalizedPositionsContainer> newPositionsContainer) {
	// are we trying to add the same container to itself?
	if (this == newPositionsContainer.get()) {
		throw std::runtime_error("Trying to append a LocalizedPositionsContainer_Multiplication to itself");
	}
	
	// check if the positions container is of the right type
	if (newPositionsContainer->getPositionsType() != LOCALIZED_POSITIONS_TYPE_MULTIPLICATION) {
		throw std::runtime_error("Trying to append a position of a different type to a LocalizedPositionsContainer_Multiplication");
	}
	
	// cast the pointer to the more specific type
	boost::shared_ptr<LocalizedPositionsContainer_Multiplication> newPositionsContainer_2DGauss(boost::static_pointer_cast<LocalizedPositionsContainer_Multiplication> (newPositionsContainer));
	
	for (std::vector<LocalizedPosition_Multiplication>::iterator it = newPositionsContainer_2DGauss->positionsVector.begin(); it != newPositionsContainer_2DGauss->positionsVector.end(); ++it) {
		this->positionsVector.push_back(*it);
	}
}

#ifdef WITH_IGOR
LocalizedPositionsContainer_ZeissPALM::LocalizedPositionsContainer_ZeissPALM(waveHndl positionsWave) {
	// initialize a new PositionsContainer from a wave that contains positions of the correct type
	long numDimensions;
	long dimensionSizes[MAX_DIMENSIONS+1];
	int err;
	
	err = MDGetWaveDimensions(positionsWave, &numDimensions, dimensionSizes);
	
	if ((numDimensions != 2) || (dimensionSizes[1] != 6)) {	// invalid dimensions (warning: magic numbers)
		throw (std::runtime_error("Invalid positions wave"));
	}
	
	LocalizedPosition_ZeissPALM singlePosition;
	size_t nPositions = dimensionSizes[0];
	this->positionsVector.reserve(nPositions);
	long indices[MAX_DIMENSIONS];
	double value[2];
	
	for (size_t i = 0; i < nPositions; ++i) {
		// get all the relevant data out of the wave and into a position object
		indices[0] = i;
		indices[1] = 0;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.frameNumber = value[0];
		indices[1] = 1;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.integral = value[0];
		indices[1] = 2;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.xPosition = value[0];
		indices[1] = 3;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.yPosition = value[0];
		indices[1] = 4;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.positionDeviation = value[0];
		indices[1] = 5;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.nFramesPresent = value[0];
		
		this->positionsVector.push_back(singlePosition);
	}
}

waveHndl LocalizedPositionsContainer_ZeissPALM::writePositionsToWave(std::string waveName, std::string waveNote) const {
	long dimensionSizes[MAX_DIMENSIONS+1];
	int err;
	waveHndl outputWave;
	size_t nPositions = this->positionsVector.size();
	dimensionSizes[0] = nPositions;
	dimensionSizes[1] = 6;	// magic number
	dimensionSizes[2] = 0;
	err = MDMakeWave(&outputWave, waveName.c_str(), NULL, dimensionSizes, NT_FP64, 1);
	if (err != 0)
		throw err;
	
	long indices[MAX_DIMENSIONS];
	double value[2];
	
	for (size_t i = 0; i < nPositions; ++i) {
		indices[0] = i;
		indices[1] = 0;
		value[0] = this->positionsVector.at(i).frameNumber;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 1;
		value[0] = this->positionsVector.at(i).integral;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 2;
		value[0] = this->positionsVector.at(i).xPosition;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 3;
		value[0] = this->positionsVector.at(i).yPosition;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 4;
		value[0] = this->positionsVector.at(i).positionDeviation;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 5;
		value[0] = this->positionsVector.at(i).nFramesPresent;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
	}
	
	if (waveNote.size() != 0) {
		// set the waveNote to the string passed in
		Handle waveNoteHandle = NewHandle(waveNote.length());
		if (waveNoteHandle == NULL)
			throw std::bad_alloc();
		
		PutCStringInHandle(waveNote.c_str(), waveNoteHandle);
		SetWaveNote(outputWave, waveNoteHandle);
	}
	
	return outputWave;
}
#endif // WITH_IGOR

void LocalizedPositionsContainer_ZeissPALM::addPosition(boost::shared_ptr<LocalizedPosition> newPosition) {
	// check if the type of positions that we are adding is suitable
	if (newPosition->getPositionType() != LOCALIZED_POSITIONS_TYPE_ZEISSPALM)
		throw std::runtime_error("Trying to append a position of a different type to a LocalizedPositionsContainer_Multiplication");
	
	// cast the pointer to the more specific type
	boost::shared_ptr<LocalizedPosition_ZeissPALM> newPosition_ZeissPALM(boost::static_pointer_cast<LocalizedPosition_ZeissPALM> (newPosition));
	
	this->positionsVector.push_back(*newPosition_ZeissPALM);
}

void LocalizedPositionsContainer_ZeissPALM::addPositions(boost::shared_ptr<LocalizedPositionsContainer> newPositionsContainer) {
	// are we trying to add the same container to itself?
	if (this == newPositionsContainer.get()) {
		throw std::runtime_error("Trying to append a LocalizedPositionsContainer_ZeissPALM to itself");
	}
	
	// check if the positions container is of the right type
	if (newPositionsContainer->getPositionsType() != LOCALIZED_POSITIONS_TYPE_ZEISSPALM) {
		throw std::runtime_error("Trying to append a position of a different type to a LocalizedPositionsContainer_ZeissPALM");
	}
	
	// cast the pointer to the more specific type
	boost::shared_ptr<LocalizedPositionsContainer_ZeissPALM> newPositionsContainer_ZeissPALM(boost::static_pointer_cast<LocalizedPositionsContainer_ZeissPALM> (newPositionsContainer));
	
	for (std::vector<LocalizedPosition_ZeissPALM>::iterator it = newPositionsContainer_ZeissPALM->positionsVector.begin(); it != newPositionsContainer_ZeissPALM->positionsVector.end(); ++it) {
		this->positionsVector.push_back(*it);
	}
}

PALMAnalysisController::PALMAnalysisController(boost::shared_ptr<ImageLoader> imageLoader_rhs, boost::shared_ptr<ThresholdImage> thresholder_rhs,
											   boost::shared_ptr<ThresholdImage_Preprocessor> thresholdImagePreprocessor_rhs,
											   boost::shared_ptr<ThresholdImage_Postprocessor> thresholdImagePostprocessor_rhs,
											   boost::shared_ptr<ParticleFinder> particleFinder_rhs, boost::shared_ptr<FitPositions> fitPositions_rhs,
											   boost::shared_ptr<PALMAnalysisProgressReporter> progressReporter_rhs) {
	imageLoader = imageLoader_rhs;
	thresholder = thresholder_rhs;
	thresholdImagePreprocessor = thresholdImagePreprocessor_rhs;
	thresholdImagePostprocessor = thresholdImagePostprocessor_rhs;
	particleFinder = particleFinder_rhs;
	fitPositions = fitPositions_rhs;
	progressReporter = progressReporter_rhs;
	
	nImages = imageLoader->get_total_number_of_images();
	
	this->errorMessage.assign("");
}

boost::shared_ptr<LocalizedPositionsContainer> PALMAnalysisController::DoPALMAnalysis() {
	size_t numberOfProcessors = boost::thread::hardware_concurrency();
	size_t numberOfThreads;
	vector<boost::shared_ptr<boost::thread> > threads;
	boost::shared_ptr<boost::thread> singleThreadPtr;
	int firstThreadHasFinished, status;
	double percentDone;
	
	if (this->nImages == 0) {	// if there are no images to load, do not do any processing
		return boost::shared_ptr<LocalizedPositionsContainer_2DGauss> (new LocalizedPositionsContainer_2DGauss());
	}
	
	numberOfThreads = numberOfProcessors * 2;	// take two threads for every processor since every thread will be blocked on I/O sooner or later
	if (numberOfThreads == 0) {
		numberOfThreads = 1;
	}
	if (numberOfThreads > this->nImages) {
		numberOfThreads = nImages;
	}
	
	// fill the queue holding the frames to be processed with the frames in the sequence
	for (size_t i = 0; i < this->nImages; ++i) {
		this->framesToBeProcessed.push(i);
	}
	
	progressReporter->CalculationStarted();
	
	// start the thread pool
	threads.clear();
	for (size_t j = 0; j < numberOfThreads; ++j) {
		singleThreadPtr = boost::shared_ptr<boost::thread>(new boost::thread(boost::bind(ThreadPoolWorker, this)));
		threads.push_back(singleThreadPtr);
	}
	
	// test if the threads have finished
	for (;;) {
		firstThreadHasFinished = threads.at(0)->timed_join(boost::posix_time::milliseconds(500));
		if (firstThreadHasFinished == 0) {	// the thread is not done yet, we're just waiting
											// while we wait we check for various things and give the interface the chance to update
			
			// did one of the threads run into an error?
			this->errorReportingMutex.lock();
			if (this->errorMessage.length() != 0) {
				errorReportingMutex.unlock();
				// an error occured, time to abort this analysis
				for (size_t j = 0; j < numberOfThreads; ++j) {
					threads.at(j)->interrupt();
				}
				// wait until the threads have completed
				for (size_t j = 0; j < numberOfThreads; ++j) {
				threads.at(j)->join();
				}
				
				throw ERROR_RUNNING_THREADED_ANALYSIS(this->errorMessage);
			}
			this->errorReportingMutex.unlock();
			
			// does the user want to abort?
#ifdef WITH_IGOR
			status = CheckAbort(0);
			if (status == -1) {
				for (size_t j = 0; j < numberOfThreads; ++j) {
					threads.at(j)->interrupt();
				}
				// wait until the threads have completed
				for (size_t j = 0; j < numberOfThreads; ++j) {
					threads.at(j)->join();
				}
				progressReporter->CalculationAborted();
				throw USER_ABORTED("Analysis aborted on user request");
			}
#endif // WITH_IGOR
			
			// allow the reporter to update with new progress
			this->acquireFrameForProcessingMutex.lock();
			percentDone = (double)(nImages - framesToBeProcessed.size()) / (double)(this->nImages) * 100.0;
			progressReporter->UpdateCalculationProgress(percentDone);
			this->acquireFrameForProcessingMutex.unlock();
			continue;
		} else {
			break;
		}
	}
	
	for (size_t j = 1; j < numberOfThreads; ++j) {
		threads.at(j)->join();
	}
	
	// the processing itself is now done, but the results will not have been returned in the correct order
	this->localizedPositions->sortPositionsByFrameNumber();
	
	progressReporter->CalculationDone();
	
	// return the results
	return this->localizedPositions;
}

void ThreadPoolWorker(PALMAnalysisController* controller) {
	size_t currentImageToProcess;
	boost::shared_ptr<PALMMatrix<double> > currentImage;
	boost::shared_ptr<PALMMatrix <unsigned char> > thresholdedImage;
	boost::shared_ptr<PALMMatrix<double> > locatedParticles;
	boost::shared_ptr<LocalizedPositionsContainer> localizedPositions;
	
	try {
		for (;;) {	// loop continuously looking for more images until there are none left
			// if the main thread wants us to interrupt, then give it an opportunity to do so
			if (boost::this_thread::interruption_requested()) {
				return;
			}
			
			// start by a acquiring an image to process
			controller->acquireFrameForProcessingMutex.lock();
			if (controller->framesToBeProcessed.size() == 0) {
				// no more frames to be processed
				controller->acquireFrameForProcessingMutex.unlock();
				return;
			}
			currentImageToProcess = controller->framesToBeProcessed.front();
			controller->framesToBeProcessed.pop();
			controller->acquireFrameForProcessingMutex.unlock();
			
			// we need to process the image with index currentImageToProcess
			currentImage = controller->imageLoader->get_nth_image(currentImageToProcess);
			thresholdedImage = do_processing_and_thresholding(currentImage, controller->thresholdImagePreprocessor, controller->thresholder,
															  controller->thresholdImagePostprocessor);
			locatedParticles = controller->particleFinder->findPositions(currentImage, thresholdedImage);
			localizedPositions = controller->fitPositions->fit_positions(currentImage, locatedParticles);
			
			// the localization routines have no idea what frame number the processed frame was
			// so we need to add the appropriate frame number to the fitted positions here
			localizedPositions->setFrameNumbers(currentImageToProcess);
			
			// pass the result to the output queue
			controller->addLocalizedPositionsMutex.lock();
			// if this is the first time that positions are being returned then controller will contain a NULL pointer
			// set it to the positions we are now returning
			// this way ThreadPoolWorker does not need to know what the type of positions is
			if (controller->localizedPositions.get() == NULL) {
				controller->localizedPositions = localizedPositions;
			} else {
				controller->localizedPositions->addPositions(localizedPositions);
			}
			
			controller->addLocalizedPositionsMutex.unlock();
		}
	}
	catch (runtime_error &e) {
		// an error has appeared somewhere in this thread during the calculation
		// since there seems to be no easy way to communicate the exception to the
		// main thread, set an error message in the analysis controller.
		// The controller will periodically check this message and handle the error
		controller->errorReportingMutex.lock();
		controller->errorMessage.assign(e.what());
		controller->errorReportingMutex.unlock();
		
		// no point in continuing this thread
		return;
	}
	catch (boost::thread_interrupted) {
		// the main thread wants us to stop
		return;
	}
	catch (...) {
		// catch any other exception not handled by the above block
		controller->errorReportingMutex.lock();
		controller->errorMessage.assign("Encountered an unspecified exception while doing the PALM analysis");
		controller->errorReportingMutex.unlock();
		
		// no point in continuing this thread
		return;
	}
	
}

#ifdef WITH_IGOR
void PALMAnalysisProgressReporter_IgorCommandLine::UpdateCalculationProgress(double percentDone) {
	char XOPOut[10];
	if (percentDone - previousPercentage > 10.0) {
		previousPercentage = floor(percentDone / 10.0) * 10.0;
		sprintf(XOPOut, "%.0lf%% ", previousPercentage);
		XOPNotice(XOPOut);
	}
}
#endif // WITH_IGOR

boost::shared_ptr<vector<double> > RipleysKFunctionCalculator::CalculateRipleysKFunction(boost::shared_ptr<PALMMatrix<double> > positions, double startBin, double endBin, double binWidth) {
	size_t nBins = (size_t)ceil((endBin - startBin) / binWidth + 0.5);
	size_t nPositions = positions->getXSize();
	
	size_t minX = 0, maxX = (size_t)-1, minY = 0, maxY = (size_t)-1;	// the outermost positions
	
	// determine the coordinates of the box enclosing the points
	for (size_t i = 0; i < nPositions; ++i) {
		if ((*positions)(i, 3) < minX)	// MAGIC NUMBERS
			minX = (*positions)(i, 3);
		if ((*positions)(i, 3) > maxX)
			maxX = (*positions)(i, 3);
		if ((*positions)(i, 4) < minY)
			minY = (*positions)(i, 3);
		if ((*positions)(i, 4) > maxY)
			maxY = (*positions)(i, 3);
	}
	
	// allocate the output memory
	boost::shared_ptr<vector <double> > kFunction(new vector<double>());
	kFunction->resize(nBins, double());
	
	// the main and slow calculation loop
	// calculate the distance between all points
	double distance;
	double lowerLimit = startBin - binWidth / 2.0;
	double upperLimit = endBin + binWidth / 2.0;
	size_t currentBin;
	for (size_t i = 0; i < nPositions; ++i) {
		for (size_t j = i + 1; j < nPositions; ++j) {
			distance = sqrt(((*positions)(i, 3) - (*positions)(j, 3)) * ((*positions)(i, 3) - (*positions)(j, 3)) + ((*positions)(i, 4) - (*positions)(j, 4)) * ((*positions)(i, 4) - (*positions)(j, 4)));
			// check if the distance between these points is within the limit that we want to support
			if ((distance < lowerLimit) || (distance > upperLimit)) {
				continue;
			}
			
			currentBin = ceil((distance - lowerLimit) / binWidth);
			// we need to increment kFunction at currentBin and all higher bins
			// the base value is 2 (since the same pair occurs twice in the equations in literature
			// but we only calculate each pair once (summation conditions are different)
			
			// however, we need to apply edge correction. We will choose the correction due to Ripley,
			// which states that we weigh with the proportion of the circumreference of a circle
			// centered at the point
		}
	}
}

