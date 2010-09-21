/*
 *  PALM_analysis_storage.cpp
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 18/04/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "PALM_analysis_storage.h"

#ifdef WITH_IGOR

boost::shared_ptr<LocalizedPositionsContainer> LocalizedPositionsContainer::GetPositionsFromWave(waveHndl positionsWave) {
	int err;
	size_t findPosition;
	
	// determine the type of positions being passed
	Handle waveNoteHandle = WaveNote(positionsWave);
	size_t waveNoteSize = GetHandleSize(positionsWave);
	if (waveNoteSize == 0) {	// no wavenote
		throw std::runtime_error("The localized positions wave does not contain a wavenote, cannot determine the storage type");
	}
	
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
	
	findPosition = waveNote.find("LOCALIZATION METHOD:5");
	if (findPosition != (size_t)-1) {
		return boost::shared_ptr<LocalizedPositionsContainer> (new LocalizedPositionsContainer_Ellipsoidal2DGaussian(positionsWave));
	}
	
	findPosition = waveNote.find("LOCALIZATION METHOD:6");
	if (findPosition != (size_t)-1) {
		return boost::shared_ptr<LocalizedPositionsContainer> (new LocalizedPositionsContainer_MLEwG(positionsWave));
	}
	
	// if we are still here then we don't recognize the type of localization used
	throw std::runtime_error("Unknown localization method (check the wave note of the wave containing the positions)");
}
#endif // #ifdef WITH_IGOR

void LocalizedPositionsContainer_2DGauss::addPosition(boost::shared_ptr<LocalizedPosition> newPosition) {
	// check if the type of positions that we are adding is suitable
	if (newPosition->getPositionType() != this->getPositionsType())
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
	if (newPositionsContainer->getPositionsType() != this->getPositionsType()) {
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
	int numDimensions;
	CountInt dimensionSizes[MAX_DIMENSIONS+1];
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

waveHndl LocalizedPositionsContainer_2DGauss::writePositionsToWave(DataFolderAndName outputWaveParams, std::string waveNote) const {
	long dimensionSizes[MAX_DIMENSIONS+1];
	int err;
	waveHndl outputWave;
	size_t nPositions = this->positionsVector.size();
	dimensionSizes[0] = nPositions;
	dimensionSizes[1] = 12;	// magic number
	dimensionSizes[2] = 0;
	
	err = MDMakeWave(&outputWave, outputWaveParams.name, outputWaveParams.dfH, dimensionSizes, NT_FP64, 1);
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
	
	outputFile.open(filePath, std::ios::trunc);
	if (outputFile.fail()) {
		throw std::runtime_error("Unable to create the output file");
	}
	
	// write the header. Add an extra newline to be safe.
	outputFile << "Localized positions using symmetric 2D Gauss fitting" << "\n";
	outputFile << header << "\n";
	outputFile << "DATA FOLLOWS" << "\n";
	outputFile << "First frame" << "\t";
	outputFile << "Integrated intensity" << "\t";
	outputFile << "Fitted PSF standard deviation (pixel)" << "\t";
	outputFile << "X position (pixel)" << "\t";
	outputFile << "Y position (pixel)" << "\t";
	outputFile << "Background" << "\t";
	outputFile << "Intensity deviation" << "\t";
	outputFile << "PSF width deviation (pixel)" << "\t";
	outputFile << "X position deviation (pixel)" << "\t";
	outputFile << "Y position deviation (pixel)" << "\t";
	outputFile << "Background deviation (pixel)" << "\t";
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
		outputFile << (*it).nFramesPresent << "\n";
		
		if (outputFile.fail()) {
			if (outputFile.is_open()) {
				outputFile.close();
			}
			throw std::runtime_error("Error writing to the localized positions file");
		}
	}
	
	outputFile.close();
}

#ifdef WITH_IGOR
LocalizedPositionsContainer_2DGaussFixedWidth::LocalizedPositionsContainer_2DGaussFixedWidth(waveHndl positionsWave) {
	// initialize a new PositionsContainer from a wave that contains positions of the correct type
	int numDimensions;
	CountInt dimensionSizes[MAX_DIMENSIONS+1];
	int err;
	
	err = MDGetWaveDimensions(positionsWave, &numDimensions, dimensionSizes);
	
	if ((numDimensions != 2) || (dimensionSizes[1] != 10)) {	// invalid dimensions (warning: magic numbers)
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
		singlePosition.xPosition = value[0];
		indices[1] = 3;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.yPosition = value[0];
		indices[1] = 4;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.background = value[0];
		indices[1] = 5;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.integralDeviation = value[0];
		indices[1] = 6;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.xPositionDeviation = value[0];
		indices[1] = 7;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.yPositionDeviation = value[0];
		indices[1] = 8;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.backgroundDeviation = value[0];
		indices[1] = 9;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.nFramesPresent = value[0];
		
		this->positionsVector.push_back(singlePosition);
	}
}

waveHndl LocalizedPositionsContainer_2DGaussFixedWidth::writePositionsToWave(DataFolderAndName outputWaveParams, std::string waveNote) const {
	long dimensionSizes[MAX_DIMENSIONS+1];
	int err;
	waveHndl outputWave;
	size_t nPositions = this->positionsVector.size();
	dimensionSizes[0] = nPositions;
	dimensionSizes[1] = 10;	// magic number
	dimensionSizes[2] = 0;
	
	err = MDMakeWave(&outputWave, outputWaveParams.name, outputWaveParams.dfH, dimensionSizes, NT_FP64, 1);
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
		value[0] = this->positionsVector.at(i).background;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 5;
		value[0] = this->positionsVector.at(i).integralDeviation;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 6;
		value[0] = this->positionsVector.at(i).xPositionDeviation;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 7;
		value[0] = this->positionsVector.at(i).yPositionDeviation;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 8;
		value[0] = this->positionsVector.at(i).backgroundDeviation;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 9;
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
	
	outputFile.open(filePath, std::ios::trunc);
	if (outputFile.fail()) {
		throw std::runtime_error("Unable to create the output file");
	}
	
	// write the header. Add an extra newline to be safe.
	outputFile << "Localized positions using symmetric 2D Gauss fitting with fixed PSF width" << "\n";
	outputFile << header << "\n";
	outputFile << "DATA FOLLOWS" << "\n";
	outputFile << "First frame" << "\t";
	outputFile << "Integrated intensity" << "\t";
	outputFile << "X position (pixel)" << "\t";
	outputFile << "Y position (pixel)" << "\t";
	outputFile << "Background" << "\t";
	outputFile << "Intensity deviation" << "\t";
	outputFile << "X position deviation (pixel)" << "\t";
	outputFile << "Y position deviation (pixel)" << "\t";
	outputFile << "Background deviation" << "\t";
	outputFile << "Number of frames where this emitter is present" << "\n";
	
	// write the actual positions
	for (std::vector<LocalizedPosition_2DGaussFixedWidth>::const_iterator it = this->positionsVector.begin(); it != this->positionsVector.end(); ++it) {
		outputFile << (*it).frameNumber << "\t";
		outputFile << (*it).integral << "\t";
		outputFile << (*it).xPosition << "\t";
		outputFile << (*it).yPosition << "\t";
		outputFile << (*it).background << "\t";
		outputFile << (*it).integralDeviation << "\t";
		outputFile << (*it).xPositionDeviation << "\t";
		outputFile << (*it).yPositionDeviation << "\t";
		outputFile << (*it).backgroundDeviation << "\t";
		outputFile << (*it).nFramesPresent << "\n";
		
		if (outputFile.fail()) {
			if (outputFile.is_open()) {
				outputFile.close();
			}
			throw std::runtime_error("Error writing to the localized positions file");
		}
	}
	
	outputFile.close();
}

void LocalizedPositionsContainer_2DGaussFixedWidth::addPosition(boost::shared_ptr<LocalizedPosition> newPosition) {
	// check if the type of positions that we are adding is suitable
	if (newPosition->getPositionType() != this->getPositionsType())
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
	if (newPositionsContainer->getPositionsType() != this->getPositionsType()) {
		throw std::runtime_error("Trying to append a position of a different type to a LocalizedPositionsContainer_2DGauss");
	}
	
	// cast the pointer to the more specific type
	boost::shared_ptr<LocalizedPositionsContainer_2DGaussFixedWidth> newPositionsContainer_2DGaussFixedWidth(boost::static_pointer_cast<LocalizedPositionsContainer_2DGaussFixedWidth> (newPositionsContainer));
	
	for (std::vector<LocalizedPosition_2DGaussFixedWidth>::iterator it = newPositionsContainer_2DGaussFixedWidth->positionsVector.begin(); it != newPositionsContainer_2DGaussFixedWidth->positionsVector.end(); ++it) {
		this->positionsVector.push_back(*it);
	}
}

void LocalizedPositionsContainer_Ellipsoidal2DGaussian::addPosition(boost::shared_ptr<LocalizedPosition> newPosition) {
	// check if the type of positions that we are adding is suitable
	if (newPosition->getPositionType() != this->getPositionsType())
		throw std::runtime_error("Trying to append a position of a different type to a LocalizedPositionsContainer_2DGauss");
	
	// cast the pointer to the more specific type
	boost::shared_ptr<LocalizedPosition_Ellipsoidal2DGauss> newPosition_Ellipsoidal2DGauss(boost::static_pointer_cast<LocalizedPosition_Ellipsoidal2DGauss> (newPosition));
	
	this->positionsVector.push_back(*newPosition_Ellipsoidal2DGauss);
}

void LocalizedPositionsContainer_Ellipsoidal2DGaussian::addPositions(boost::shared_ptr<LocalizedPositionsContainer> newPositionsContainer) {
	// are we trying to add the same container to itself?
	if (this == newPositionsContainer.get()) {
		throw std::runtime_error("Trying to append a LocalizedPositionsContainer_2DGauss to itself");
	}
	
	// check if the positions container is of the right type
	if (newPositionsContainer->getPositionsType() != this->getPositionsType()) {
		throw std::runtime_error("Trying to append a position of a different type to a LocalizedPositionsContainer_2DGauss");
	}
	
	// cast the pointer to the more specific type
	boost::shared_ptr<LocalizedPositionsContainer_Ellipsoidal2DGaussian> newPositionsContainer_Ellipsoidal2DGauss(boost::static_pointer_cast<LocalizedPositionsContainer_Ellipsoidal2DGaussian> (newPositionsContainer));
	
	for (std::vector<LocalizedPosition_Ellipsoidal2DGauss>::iterator it = newPositionsContainer_Ellipsoidal2DGauss->positionsVector.begin(); it != newPositionsContainer_Ellipsoidal2DGauss->positionsVector.end(); ++it) {
		this->positionsVector.push_back(*it);
	}
}

#ifdef WITH_IGOR
LocalizedPositionsContainer_Ellipsoidal2DGaussian::LocalizedPositionsContainer_Ellipsoidal2DGaussian(waveHndl positionsWave) {
	// initialize a new PositionsContainer from a wave that contains positions of the correct type
	int numDimensions;
	CountInt dimensionSizes[MAX_DIMENSIONS+1];
	int err;
	
	err = MDGetWaveDimensions(positionsWave, &numDimensions, dimensionSizes);
	
	if ((numDimensions != 2) || (dimensionSizes[1] != 16)) {	// invalid dimensions (warning: magic numbers)
		throw (std::runtime_error("Invalid positions wave"));
	}
	
	LocalizedPosition_Ellipsoidal2DGauss singlePosition;
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
		singlePosition.xWidth = value[0];
		indices[1] = 3;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.yWidth = value[0];
		indices[1] = 4;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.xPosition = value[0];
		indices[1] = 5;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.yPosition = value[0];
		indices[1] = 6;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.correlation = value[0];
		indices[1] = 7;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.background = value[0];
		indices[1] = 8;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.integralDeviation = value[0];
		indices[1] = 9;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.xWidthDeviation = value[0];
		indices[1] = 10;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.yWidthDeviation = value[0];
		indices[1] = 11;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.xPositionDeviation = value[0];
		indices[1] = 12;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.yPositionDeviation = value[0];
		indices[1] = 13;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.correlationDeviation = value[0];
		indices[1] = 14;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.backgroundDeviation = value[0];
		indices[1] = 15;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.nFramesPresent = value[0];
		
		this->positionsVector.push_back(singlePosition);
	}
}

waveHndl LocalizedPositionsContainer_Ellipsoidal2DGaussian::writePositionsToWave(DataFolderAndName outputWaveParams, std::string waveNote) const {
	long dimensionSizes[MAX_DIMENSIONS+1];
	int err;
	waveHndl outputWave;
	size_t nPositions = this->positionsVector.size();
	dimensionSizes[0] = nPositions;
	dimensionSizes[1] = 16;	// magic number
	dimensionSizes[2] = 0;
	
	err = MDMakeWave(&outputWave, outputWaveParams.name, outputWaveParams.dfH, dimensionSizes, NT_FP64, 1);
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
		value[0] = this->positionsVector.at(i).xWidth;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 3;
		value[0] = this->positionsVector.at(i).yWidth;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 4;
		value[0] = this->positionsVector.at(i).xPosition;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 5;
		value[0] = this->positionsVector.at(i).yPosition;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 6;
		value[0] = this->positionsVector.at(i).correlation;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 7;
		value[0] = this->positionsVector.at(i).background;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 8;
		value[0] = this->positionsVector.at(i).integralDeviation;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 9;
		value[0] = this->positionsVector.at(i).xWidthDeviation;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 10;
		value[0] = this->positionsVector.at(i).yWidthDeviation;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 11;
		value[0] = this->positionsVector.at(i).xPositionDeviation;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 12;
		value[0] = this->positionsVector.at(i).yPositionDeviation;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 13;
		value[0] = this->positionsVector.at(i).correlationDeviation;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 14;
		value[0] = this->positionsVector.at(i).backgroundDeviation;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 15;
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

void LocalizedPositionsContainer_Ellipsoidal2DGaussian::writePositionsToFile(std::string filePath, std::string header) const {
	boost::filesystem::ofstream outputFile;
	size_t nPositions = this->positionsVector.size();
	if (nPositions == 0) {
		throw std::runtime_error("No positions localized!");
	}
	
	outputFile.open(filePath, std::ios::trunc);
	if (outputFile.fail()) {
		throw std::runtime_error("Unable to create the output file");
	}
	
	// write the header. Add an extra newline to be safe.
	outputFile << "Localized positions using ellipsoidal 2D Gauss fitting" << "\n";
	outputFile << header << "\n";
	outputFile << "DATA FOLLOWS" << "\n";
	outputFile << "First frame" << "\t";
	outputFile << "Integrated intensity" << "\t";
	outputFile << "Fitted PSF standard deviation along x (pixel)" << "\t";
	outputFile << "Fitted PSF standard deviation along y (pixel)" << "\t";
	outputFile << "X position (pixel)" << "\t";
	outputFile << "Y position (pixel)" << "\t";
	outputFile << "Correlation between x and y" << "\t";
	outputFile << "Background" << "\t";
	outputFile << "Intensity deviation" << "\t";
	outputFile << "PSF width deviation along x (pixel)" << "\t";
	outputFile << "PSF width deviation along y (pixel)" << "\t";
	outputFile << "X position deviation (pixel)" << "\t";
	outputFile << "Y position deviation (pixel)" << "\t";
	outputFile << "Correlation deviation" << "\t";
	outputFile << "Background deviation" << "\t";
	outputFile << "Number of frames where this emitter is present" << "\n";
	
	// write the actual positions
	for (std::vector<LocalizedPosition_Ellipsoidal2DGauss>::const_iterator it = this->positionsVector.begin(); it != this->positionsVector.end(); ++it) {
		outputFile << (*it).frameNumber << "\t";
		outputFile << (*it).integral << "\t";
		outputFile << (*it).xWidth << "\t";
		outputFile << (*it).yWidth << "\t";
		outputFile << (*it).xPosition << "\t";
		outputFile << (*it).yPosition << "\t";
		outputFile << (*it).correlation << "\t";
		outputFile << (*it).background << "\t";
		outputFile << (*it).integralDeviation << "\t";
		outputFile << (*it).xWidthDeviation << "\t";
		outputFile << (*it).yWidthDeviation << "\t";
		outputFile << (*it).xPositionDeviation << "\t";
		outputFile << (*it).yPositionDeviation << "\t";
		outputFile << (*it).correlationDeviation << "\t";
		outputFile << (*it).backgroundDeviation << "\t";
		outputFile << (*it).nFramesPresent << "\n";
		
		if (outputFile.fail()) {
			if (outputFile.is_open()) {
				outputFile.close();
			}
			throw std::runtime_error("Error writing to the localized positions file");
		}
	}
	
	outputFile.close();
}

#ifdef WITH_IGOR
LocalizedPositionsContainer_Centroid::LocalizedPositionsContainer_Centroid(waveHndl positionsWave) {
	// initialize a new PositionsContainer from a wave that contains positions of the correct type
	int numDimensions;
	CountInt dimensionSizes[MAX_DIMENSIONS+1];
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

waveHndl LocalizedPositionsContainer_Centroid::writePositionsToWave(DataFolderAndName outputWaveParams, std::string waveNote) const {
	long dimensionSizes[MAX_DIMENSIONS+1];
	int err;
	waveHndl outputWave;
	size_t nPositions = this->positionsVector.size();
	dimensionSizes[0] = nPositions;
	dimensionSizes[1] = 4;	// magic number
	dimensionSizes[2] = 0;
	
	err = MDMakeWave(&outputWave, outputWaveParams.name, outputWaveParams.dfH, dimensionSizes, NT_FP64, 1);
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
	
	outputFile.open(filePath, std::ios::trunc);
	if (outputFile.fail()) {
		throw std::runtime_error("Unable to create the output file");
	}
	
	// write the header. Add an extra newline to be safe.
	outputFile << "Localized positions using centroid calculation" << "\n";
	outputFile << header << "\n";
	
	outputFile << "DATA FOLLOWS" << "\n";
	outputFile << "First frame" << "\t";
	outputFile << "X position (pixel)" << "\t";
	outputFile << "Y position (pixel)" << "\t";
	outputFile << "Number of frames where this emitter is present" << "\n";
	
	// write the actual positions
	for (std::vector<LocalizedPosition_Centroid>::const_iterator it = this->positionsVector.begin(); it != this->positionsVector.end(); ++it) {
		outputFile << (*it).frameNumber << "\t";
		outputFile << (*it).xPosition << "\t";
		outputFile << (*it).yPosition << "\t";
		outputFile << (*it).nFramesPresent << "\n";
		
		if (outputFile.fail()) {
			if (outputFile.is_open()) {
				outputFile.close();
			}
			throw std::runtime_error("Error writing to the localized positions file");
		}
	}
	
	outputFile.close();
}

void LocalizedPositionsContainer_Centroid::addPosition(boost::shared_ptr<LocalizedPosition> newPosition) {
	// check if the type of positions that we are adding is suitable
	if (newPosition->getPositionType() != this->getPositionsType())
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
	int numDimensions;
	CountInt dimensionSizes[MAX_DIMENSIONS+1];
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

waveHndl LocalizedPositionsContainer_Multiplication::writePositionsToWave(DataFolderAndName outputWaveParams, std::string waveNote) const {
	long dimensionSizes[MAX_DIMENSIONS+1];
	int err;
	waveHndl outputWave;
	size_t nPositions = this->positionsVector.size();
	dimensionSizes[0] = nPositions;
	dimensionSizes[1] = 5;	// magic number
	dimensionSizes[2] = 0;
	
	err = MDMakeWave(&outputWave, outputWaveParams.name, outputWaveParams.dfH, dimensionSizes, NT_FP64, 1);
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
	
	outputFile.open(filePath, std::ios::trunc);
	if (outputFile.fail()) {
		throw std::runtime_error("Unable to create the output file");
	}
	
	// write the header. Add an extra newline to be safe.
	outputFile << "Localized positions using iterative multiplication" << "\n";
	outputFile << header << "\n";
	outputFile << "DATA FOLLOWS" << "\n";
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
		outputFile << (*it).nFramesPresent << "\n";
		
		if (outputFile.fail()) {
			if (outputFile.is_open()) {
				outputFile.close();
			}
			throw std::runtime_error("Error writing to the localized positions file");
		}
	}
	
	outputFile.close();
}

void LocalizedPositionsContainer_Multiplication::addPosition(boost::shared_ptr<LocalizedPosition> newPosition) {
	// check if the type of positions that we are adding is suitable
	if (newPosition->getPositionType() != this->getPositionsType())
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
	int numDimensions;
	CountInt dimensionSizes[MAX_DIMENSIONS+1];
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

waveHndl LocalizedPositionsContainer_ZeissPALM::writePositionsToWave(DataFolderAndName outputWaveParams, std::string waveNote) const {
	long dimensionSizes[MAX_DIMENSIONS+1];
	int err;
	waveHndl outputWave;
	size_t nPositions = this->positionsVector.size();
	dimensionSizes[0] = nPositions;
	dimensionSizes[1] = 6;	// magic number
	dimensionSizes[2] = 0;
	
	err = MDMakeWave(&outputWave, outputWaveParams.name, outputWaveParams.dfH, dimensionSizes, NT_FP64, 1);
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
	if (newPosition->getPositionType() != this->getPositionsType())
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

#ifdef WITH_IGOR
LocalizedPositionsContainer_MLEwG::LocalizedPositionsContainer_MLEwG(waveHndl positionsWave) {
	// initialize a new PositionsContainer from a wave that contains positions of the correct type
	int numDimensions;
	CountInt dimensionSizes[MAX_DIMENSIONS+1];
	int err;
	
	err = MDGetWaveDimensions(positionsWave, &numDimensions, dimensionSizes);
	
	if ((numDimensions != 2) || (dimensionSizes[1] != 8)) {	// invalid dimensions (warning: magic numbers)
		throw (std::runtime_error("Invalid positions wave"));
	}
	
	LocalizedPosition_MLEwG singlePosition;
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
		singlePosition.positionDeviation = value[0];
		indices[1] = 7;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.nFramesPresent = value[0];
		
		this->positionsVector.push_back(singlePosition);
	}
}

waveHndl LocalizedPositionsContainer_MLEwG::writePositionsToWave(DataFolderAndName outputWaveParams, std::string waveNote) const {
	long dimensionSizes[MAX_DIMENSIONS+1];
	int err;
	waveHndl outputWave;
	size_t nPositions = this->positionsVector.size();
	dimensionSizes[0] = nPositions;
	dimensionSizes[1] = 8;	// magic number
	dimensionSizes[2] = 0;
	
	err = MDMakeWave(&outputWave, outputWaveParams.name, outputWaveParams.dfH, dimensionSizes, NT_FP64, 1);
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
		value[0] = this->positionsVector.at(i).positionDeviation;
		err = MDSetNumericWavePointValue(outputWave, indices, value);
		indices[1] = 7;
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

void LocalizedPositionsContainer_MLEwG::writePositionsToFile(std::string filePath, std::string header) const {
	boost::filesystem::ofstream outputFile;
	size_t nPositions = this->positionsVector.size();
	if (nPositions == 0) {
		throw std::runtime_error("No positions localized!");
	}
	
	outputFile.open(filePath, std::ios::trunc);
	if (outputFile.fail()) {
		throw std::runtime_error("Unable to create the output file");
	}
	
	// write the header. Add an extra newline to be safe.
	outputFile << "Localized positions using MLEwG estimation" << "\n";
	outputFile << header << "\n";
	outputFile << "DATA FOLLOWS" << "\n";
	outputFile << "First frame" << "\t";
	outputFile << "Integrated intensity" << "\t";
	outputFile << "Fitted PSF standard deviation (pixel)" << "\t";
	outputFile << "X position (pixel)" << "\t";
	outputFile << "Y position (pixel)" << "\t";
	outputFile << "Background" << "\t";
	outputFile << "Position deviation (pixel)" << "\t";
	outputFile << "Number of frames where this emitter is present" << "\n";
	
	// write the actual positions
	for (std::vector<LocalizedPosition_MLEwG>::const_iterator it = this->positionsVector.begin(); it != this->positionsVector.end(); ++it) {
		outputFile << (*it).frameNumber << "\t";
		outputFile << (*it).integral << "\t";
		outputFile << (*it).width << "\t";
		outputFile << (*it).xPosition << "\t";
		outputFile << (*it).yPosition << "\t";
		outputFile << (*it).background << "\t";
		outputFile << (*it).positionDeviation << "\t";
		outputFile << (*it).nFramesPresent << "\n";
		
		if (outputFile.fail()) {
			if (outputFile.is_open()) {
				outputFile.close();
			}
			throw std::runtime_error("Error writing to the localized positions file");
		}
	}
	
	outputFile.close();
}

void LocalizedPositionsContainer_MLEwG::addPosition(boost::shared_ptr<LocalizedPosition> newPosition) {
	// check if the type of positions that we are adding is suitable
	if (newPosition->getPositionType() != this->getPositionsType())
		throw std::runtime_error("Trying to append a position of a different type to a LocalizedPositionsContainer_MLEwG");
	
	// cast the pointer to the more specific type
	boost::shared_ptr<LocalizedPosition_MLEwG> newPosition_MLEwG(boost::static_pointer_cast<LocalizedPosition_MLEwG> (newPosition));
	
	this->positionsVector.push_back(*newPosition_MLEwG);
}

void LocalizedPositionsContainer_MLEwG::addPositions(boost::shared_ptr<LocalizedPositionsContainer> newPositionsContainer) {
	// are we trying to add the same container to itself?
	if (this == newPositionsContainer.get()) {
		throw std::runtime_error("Trying to append a LocalizedPositionsContainer_MLEwG to itself");
	}
	
	// check if the positions container is of the right type
	if (newPositionsContainer->getPositionsType() != this->getPositionsType()) {
		throw std::runtime_error("Trying to append a position of a different type to a LocalizedPositionsContainer_MLEwG");
	}
	
	// cast the pointer to the more specific type
	boost::shared_ptr<LocalizedPositionsContainer_MLEwG> newPositionsContainer_MLEwG(boost::static_pointer_cast<LocalizedPositionsContainer_MLEwG> (newPositionsContainer));
	
	for (std::vector<LocalizedPosition_MLEwG>::iterator it = newPositionsContainer_MLEwG->positionsVector.begin(); it != newPositionsContainer_MLEwG->positionsVector.end(); ++it) {
		this->positionsVector.push_back(*it);
	}
}
