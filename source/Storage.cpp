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

#include "PALM_analysis_storage.h"

#ifdef WITH_IGOR

std::shared_ptr<LocalizedPositionsContainer> LocalizedPositionsContainer::GetPositionsFromWave(waveHndl positionsWave) {
	int err;
	size_t findPosition;
	
	// determine the type of positions being passed
	Handle waveNoteHandle = WaveNote(positionsWave);
	size_t waveNoteSize = GetHandleSize(positionsWave);
	if (waveNoteSize == 0) {	// no wavenote
		throw std::runtime_error("The localized positions wave does not contain a wavenote, cannot determine the storage type");
	}
	
	std::unique_ptr<char[]> CStringWaveNote(new char[waveNoteSize + 1]);
	
	err = GetCStringFromHandle(waveNoteHandle, CStringWaveNote.get(), waveNoteSize);
	if (err != 0) {
		if (err == USING_NULL_STRVAR) {
			throw std::runtime_error("No wavenote found. Did you make or load these positions with this software?");
		} else {
			throw std::runtime_error("GetCStringFromHandle() returned a nonzero code");
		}
	}
	
	// save the wavenote as a std::string
	std::string waveNote(CStringWaveNote.get());
	
	// see if the wave note contains info on the kind of localization used
	// if not then fail
	
	findPosition = waveNote.find("LOCALIZATION METHOD:");
	if (findPosition == (size_t)-1)	// not found
		throw std::runtime_error("The positions wave does not specify a localization method");
	
	findPosition = waveNote.find("LOCALIZATION METHOD:0;");
	if (findPosition != (size_t)-1) {
		return std::shared_ptr<LocalizedPositionsContainer> (new LocalizedPositionsContainer_2DGauss(positionsWave));
	}
	
	findPosition = waveNote.find("LOCALIZATION METHOD:1;");
	if (findPosition != (size_t)-1) {
		return std::shared_ptr<LocalizedPositionsContainer> (new LocalizedPositionsContainer_2DGaussFixedWidth(positionsWave));
	}
	
	findPosition = waveNote.find("LOCALIZATION METHOD:2");
	if (findPosition != (size_t)-1) {
		return std::shared_ptr<LocalizedPositionsContainer> (new LocalizedPositionsContainer_Multiplication(positionsWave));
	}
	
	findPosition = waveNote.find("LOCALIZATION METHOD:3");
	if (findPosition != (size_t)-1) {
		return std::shared_ptr<LocalizedPositionsContainer> (new LocalizedPositionsContainer_Centroid(positionsWave));
	}
	
	findPosition = waveNote.find("LOCALIZATION METHOD:4");
	if (findPosition != (size_t)-1) {
		return std::shared_ptr<LocalizedPositionsContainer> (new LocalizedPositionsContainer_ZeissPALM(positionsWave));
	}
	
	findPosition = waveNote.find("LOCALIZATION METHOD:5");
	if (findPosition != (size_t)-1) {
		return std::shared_ptr<LocalizedPositionsContainer> (new LocalizedPositionsContainer_Ellipsoidal2DGaussian(positionsWave));
	}
	
	findPosition = waveNote.find("LOCALIZATION METHOD:6");
	if (findPosition != (size_t)-1) {
		return std::shared_ptr<LocalizedPositionsContainer> (new LocalizedPositionsContainer_MLEwG(positionsWave));
	}
	
	// if we are still here then we don't recognize the type of localization used
	throw std::runtime_error("Unknown localization method (check the wave note of the wave containing the positions)");
}
#endif // #ifdef WITH_IGOR

void LocalizedPositionsContainer_2DGauss::addPosition(std::shared_ptr<LocalizedPosition> newPosition) {
	// check if the type of positions that we are adding is suitable
	if (newPosition->getPositionType() != this->getPositionsType())
		throw std::runtime_error("Trying to append a position of a different type to a LocalizedPositionsContainer_2DGauss");
	
	// cast the pointer to the more specific type
	std::shared_ptr<LocalizedPosition_2DGauss> newPosition_2DGauss(std::static_pointer_cast<LocalizedPosition_2DGauss> (newPosition));
	
	this->positionsVector.push_back(*newPosition_2DGauss);
}

void LocalizedPositionsContainer_2DGauss::addPositions(std::shared_ptr<LocalizedPositionsContainer> newPositionsContainer) {
	// are we trying to add the same container to itself?
	if (this == newPositionsContainer.get()) {
		throw std::runtime_error("Trying to append a LocalizedPositionsContainer_2DGauss to itself");
	}
	
	// check if the positions container is of the right type
	if (newPositionsContainer->getPositionsType() != this->getPositionsType()) {
		throw std::runtime_error("Trying to append a container of a different type to a LocalizedPositionsContainer_2DGauss");
	}
	
	// cast the pointer to the more specific type
	std::shared_ptr<LocalizedPositionsContainer_2DGauss> newPositionsContainer_2DGauss(std::static_pointer_cast<LocalizedPositionsContainer_2DGauss> (newPositionsContainer));
	
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
	CountInt indices[MAX_DIMENSIONS];
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
#endif // #ifdef WITH_IGOR

ImagePtr LocalizedPositionsContainer_2DGauss::getLocalizedPositionsAsMatrix() const {
	size_t nPositions = this->positionsVector.size();
	
	ImagePtr matrix(new Image(static_cast<int>(nPositions), 12));	// magic number
	matrix->setConstant(0.0);
	
	for (size_t i = 0; i < nPositions; ++i) {
		(*matrix)(i, 0) = this->positionsVector.at(i).frameNumber;
		(*matrix)(i, 1) = this->positionsVector.at(i).integral;
		(*matrix)(i, 2) = this->positionsVector.at(i).width;
		(*matrix)(i, 3) = this->positionsVector.at(i).xPosition;
		(*matrix)(i, 4) = this->positionsVector.at(i).yPosition;
		(*matrix)(i, 5) = this->positionsVector.at(i).background;
		(*matrix)(i, 6) = this->positionsVector.at(i).integralDeviation;
		(*matrix)(i, 7) = this->positionsVector.at(i).widthDeviation;
		(*matrix)(i, 8) = this->positionsVector.at(i).xPositionDeviation;
		(*matrix)(i, 9) = this->positionsVector.at(i).yPositionDeviation;
		(*matrix)(i, 10) = this->positionsVector.at(i).backgroundDeviation;
		(*matrix)(i, 11) = this->positionsVector.at(i).nFramesPresent;
	}
	
	return matrix;
}

void LocalizedPositionsContainer_2DGauss::writePositionsToFile(std::string filePath, std::string header) const {
	std::ofstream outputFile;
	size_t nPositions = this->positionsVector.size();
	if (nPositions == 0) {
		throw std::runtime_error("No positions localized!");
	}
	
	outputFile.open(filePath.c_str(), std::ios::trunc);
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
	CountInt indices[MAX_DIMENSIONS];
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
#endif // #ifdef WITH_IGOR

ImagePtr LocalizedPositionsContainer_2DGaussFixedWidth::getLocalizedPositionsAsMatrix() const {
	size_t nPositions = this->positionsVector.size();
	
	ImagePtr matrix(new Image(static_cast<int>(nPositions), 10));	// magic number
	matrix->setConstant(0.0);
	
	for (size_t i = 0; i < nPositions; ++i) {
		(*matrix)(i, 0) = this->positionsVector.at(i).frameNumber;
		(*matrix)(i, 1) = this->positionsVector.at(i).integral;
		(*matrix)(i, 2) = this->positionsVector.at(i).xPosition;
		(*matrix)(i, 3) = this->positionsVector.at(i).yPosition;
		(*matrix)(i, 4) = this->positionsVector.at(i).background;
		(*matrix)(i, 5) = this->positionsVector.at(i).integralDeviation;
		(*matrix)(i, 6) = this->positionsVector.at(i).xPositionDeviation;
		(*matrix)(i, 7) = this->positionsVector.at(i).yPositionDeviation;
		(*matrix)(i, 8) = this->positionsVector.at(i).backgroundDeviation;
		(*matrix)(i, 9) = this->positionsVector.at(i).nFramesPresent;
	}
	
	return matrix;
}

void LocalizedPositionsContainer_2DGaussFixedWidth::writePositionsToFile(std::string filePath, std::string header) const {
	std::ofstream outputFile;
	size_t nPositions = this->positionsVector.size();
	if (nPositions == 0) {
		throw std::runtime_error("No positions localized!");
	}
	
	outputFile.open(filePath.c_str(), std::ios::trunc);
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

void LocalizedPositionsContainer_2DGaussFixedWidth::addPosition(std::shared_ptr<LocalizedPosition> newPosition) {
	// check if the type of positions that we are adding is suitable
	if (newPosition->getPositionType() != this->getPositionsType())
		throw std::runtime_error("Trying to append a position of a different type to a LocalizedPositionsContainer_2DGauss");
	
	// cast the pointer to the more specific type
	std::shared_ptr<LocalizedPosition_2DGaussFixedWidth> newPosition_2DGaussFixedWidth(std::static_pointer_cast<LocalizedPosition_2DGaussFixedWidth> (newPosition));
	
	this->positionsVector.push_back(*newPosition_2DGaussFixedWidth);
}

void LocalizedPositionsContainer_2DGaussFixedWidth::addPositions(std::shared_ptr<LocalizedPositionsContainer> newPositionsContainer) {
	// are we trying to add the same container to itself?
	if (this == newPositionsContainer.get()) {
		throw std::runtime_error("Trying to append a LocalizedPositionsContainer_2DGaussFixedWidth to itself");
	}
	
	// check if the positions container is of the right type
	if (newPositionsContainer->getPositionsType() != this->getPositionsType()) {
		throw std::runtime_error("Trying to append a container of a different type to a LocalizedPositionsContainer_2DGauss");
	}
	
	// cast the pointer to the more specific type
	std::shared_ptr<LocalizedPositionsContainer_2DGaussFixedWidth> newPositionsContainer_2DGaussFixedWidth(std::static_pointer_cast<LocalizedPositionsContainer_2DGaussFixedWidth> (newPositionsContainer));
	
	for (std::vector<LocalizedPosition_2DGaussFixedWidth>::iterator it = newPositionsContainer_2DGaussFixedWidth->positionsVector.begin(); it != newPositionsContainer_2DGaussFixedWidth->positionsVector.end(); ++it) {
		this->positionsVector.push_back(*it);
	}
}

void LocalizedPositionsContainer_Ellipsoidal2DGaussian::addPosition(std::shared_ptr<LocalizedPosition> newPosition) {
	// check if the type of positions that we are adding is suitable
	if (newPosition->getPositionType() != this->getPositionsType())
		throw std::runtime_error("Trying to append a position of a different type to a LocalizedPositionsContainer_2DGauss");
	
	// cast the pointer to the more specific type
	std::shared_ptr<LocalizedPosition_Ellipsoidal2DGauss> newPosition_Ellipsoidal2DGauss(std::static_pointer_cast<LocalizedPosition_Ellipsoidal2DGauss> (newPosition));
	
	this->positionsVector.push_back(*newPosition_Ellipsoidal2DGauss);
}

void LocalizedPositionsContainer_Ellipsoidal2DGaussian::addPositions(std::shared_ptr<LocalizedPositionsContainer> newPositionsContainer) {
	// are we trying to add the same container to itself?
	if (this == newPositionsContainer.get()) {
		throw std::runtime_error("Trying to append a LocalizedPositionsContainer_2DGauss to itself");
	}
	
	// check if the positions container is of the right type
	if (newPositionsContainer->getPositionsType() != this->getPositionsType()) {
		throw std::runtime_error("Trying to append a container of a different type to a LocalizedPositionsContainer_2DGauss");
	}
	
	// cast the pointer to the more specific type
	std::shared_ptr<LocalizedPositionsContainer_Ellipsoidal2DGaussian> newPositionsContainer_Ellipsoidal2DGauss(std::static_pointer_cast<LocalizedPositionsContainer_Ellipsoidal2DGaussian> (newPositionsContainer));
	
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
	CountInt indices[MAX_DIMENSIONS];
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
		singlePosition.stdDev1 = value[0];
		indices[1] = 3;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.stdDev2 = value[0];
		indices[1] = 4;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.xPosition = value[0];
		indices[1] = 5;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.yPosition = value[0];
		indices[1] = 6;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.theta = value[0];
		indices[1] = 7;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.background = value[0];
		indices[1] = 8;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.integralDeviation = value[0];
		indices[1] = 9;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.stdDev1Deviation = value[0];
		indices[1] = 10;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.stdDev2Deviation = value[0];
		indices[1] = 11;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.xPositionDeviation = value[0];
		indices[1] = 12;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.yPositionDeviation = value[0];
		indices[1] = 13;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.thetaDeviation = value[0];
		indices[1] = 14;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.backgroundDeviation = value[0];
		indices[1] = 15;
		err = MDGetNumericWavePointValue(positionsWave, indices, value);
		singlePosition.nFramesPresent = value[0];
		
		this->positionsVector.push_back(singlePosition);
	}
}
#endif // #ifdef WITH_IGOR

ImagePtr LocalizedPositionsContainer_Ellipsoidal2DGaussian::getLocalizedPositionsAsMatrix() const {
	size_t nPositions = this->positionsVector.size();
	
	ImagePtr matrix(new Image(static_cast<int>(nPositions), 16));	// magic number
	matrix->setConstant(0.0);
	
	for (size_t i = 0; i < nPositions; ++i) {
		(*matrix)(i, 0) = this->positionsVector.at(i).frameNumber;
		(*matrix)(i, 1) = this->positionsVector.at(i).integral;
		(*matrix)(i, 2) = this->positionsVector.at(i).stdDev1;
		(*matrix)(i, 3) = this->positionsVector.at(i).stdDev2;
		(*matrix)(i, 4) = this->positionsVector.at(i).xPosition;
		(*matrix)(i, 5) = this->positionsVector.at(i).yPosition;
		(*matrix)(i, 6) = this->positionsVector.at(i).theta;
		(*matrix)(i, 7) = this->positionsVector.at(i).background;
		(*matrix)(i, 8) = this->positionsVector.at(i).integralDeviation;
		(*matrix)(i, 9) = this->positionsVector.at(i).stdDev1Deviation;
		(*matrix)(i, 10) = this->positionsVector.at(i).stdDev2Deviation;
		(*matrix)(i, 11) = this->positionsVector.at(i).xPositionDeviation;
		(*matrix)(i, 12) = this->positionsVector.at(i).yPositionDeviation;
		(*matrix)(i, 13) = this->positionsVector.at(i).thetaDeviation;
		(*matrix)(i, 14) = this->positionsVector.at(i).backgroundDeviation;
		(*matrix)(i, 15) = this->positionsVector.at(i).nFramesPresent;
	}
	
	return matrix;
}

void LocalizedPositionsContainer_Ellipsoidal2DGaussian::writePositionsToFile(std::string filePath, std::string header) const {
	std::ofstream outputFile;
	size_t nPositions = this->positionsVector.size();
	if (nPositions == 0) {
		throw std::runtime_error("No positions localized!");
	}
	
	outputFile.open(filePath.c_str(), std::ios::trunc);
	if (outputFile.fail()) {
		throw std::runtime_error("Unable to create the output file");
	}
	
	// write the header. Add an extra newline to be safe.
	outputFile << "Localized positions using ellipsoidal 2D Gauss fitting" << "\n";
	outputFile << header << "\n";
	outputFile << "DATA FOLLOWS" << "\n";
	outputFile << "First frame" << "\t";
	outputFile << "Integrated intensity" << "\t";
	outputFile << "Fitted PSF standard deviation along first principal axis (pixel)" << "\t";
	outputFile << "Fitted PSF standard deviation along second principal axis (pixel)" << "\t";
	outputFile << "X position (pixel)" << "\t";
	outputFile << "Y position (pixel)" << "\t";
	outputFile << "Angular rotation (radians)" << "\t";
	outputFile << "Background" << "\t";
	outputFile << "Intensity deviation" << "\t";
	outputFile << "PSF width deviation along first principal axis (pixel)" << "\t";
	outputFile << "PSF width deviation along second principal axis (pixel)" << "\t";
	outputFile << "X position deviation (pixel)" << "\t";
	outputFile << "Y position deviation (pixel)" << "\t";
	outputFile << "Rotation deviation" << "\t";
	outputFile << "Background deviation" << "\t";
	outputFile << "Number of frames where this emitter is present" << "\n";
	
	// write the actual positions
	for (std::vector<LocalizedPosition_Ellipsoidal2DGauss>::const_iterator it = this->positionsVector.begin(); it != this->positionsVector.end(); ++it) {
		outputFile << (*it).frameNumber << "\t";
		outputFile << (*it).integral << "\t";
		outputFile << (*it).stdDev1 << "\t";
		outputFile << (*it).stdDev2 << "\t";
		outputFile << (*it).xPosition << "\t";
		outputFile << (*it).yPosition << "\t";
		outputFile << (*it).theta << "\t";
		outputFile << (*it).background << "\t";
		outputFile << (*it).integralDeviation << "\t";
		outputFile << (*it).stdDev1Deviation << "\t";
		outputFile << (*it).stdDev2Deviation << "\t";
		outputFile << (*it).xPositionDeviation << "\t";
		outputFile << (*it).yPositionDeviation << "\t";
		outputFile << (*it).thetaDeviation << "\t";
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
	CountInt indices[MAX_DIMENSIONS];
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
#endif // #ifdef WITH_IGOR

ImagePtr LocalizedPositionsContainer_Centroid::getLocalizedPositionsAsMatrix() const {
	size_t nPositions = this->positionsVector.size();
	
	ImagePtr matrix(new Image(static_cast<int>(nPositions), 4));	// magic number
	matrix->setConstant(0.0);
	
	for (size_t i = 0; i < nPositions; ++i) {
		(*matrix)(i, 0) = this->positionsVector.at(i).frameNumber;
		(*matrix)(i, 1) = this->positionsVector.at(i).xPosition;
		(*matrix)(i, 2) = this->positionsVector.at(i).yPosition;
		(*matrix)(i, 3) = this->positionsVector.at(i).nFramesPresent;
	}
	
	return matrix;
}

void LocalizedPositionsContainer_Centroid::writePositionsToFile(std::string filePath, std::string header) const {
	std::ofstream outputFile;
	size_t nPositions = this->positionsVector.size();
	if (nPositions == 0) {
		throw std::runtime_error("No positions localized!");
	}
	
	outputFile.open(filePath.c_str(), std::ios::trunc);
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

void LocalizedPositionsContainer_Centroid::addPosition(std::shared_ptr<LocalizedPosition> newPosition) {
	// check if the type of positions that we are adding is suitable
	if (newPosition->getPositionType() != this->getPositionsType())
		throw std::runtime_error("Trying to append a position of a different type to a LocalizedPositionsContainer_Centroid");
	
	// cast the pointer to the more specific type
	std::shared_ptr<LocalizedPosition_Centroid> newPosition_Centroid(std::static_pointer_cast<LocalizedPosition_Centroid> (newPosition));
	
	this->positionsVector.push_back(*newPosition_Centroid);
}

void LocalizedPositionsContainer_Centroid::addPositions(std::shared_ptr<LocalizedPositionsContainer> newPositionsContainer) {
	// are we trying to add the same container to itself?
	if (this == newPositionsContainer.get()) {
		throw std::runtime_error("Trying to append a LocalizedPositionsContainer_Centroid to itself");
	}
	
	// check if the positions container is of the right type
	if (newPositionsContainer->getPositionsType() != LOCALIZED_POSITIONS_TYPE_CENTROID) {
		throw std::runtime_error("Trying to append a container of a different type to a LocalizedPositionsContainer_Centroid");
	}
	
	// cast the pointer to the more specific type
	std::shared_ptr<LocalizedPositionsContainer_Centroid> newPositionsContainer_Centroid(std::static_pointer_cast<LocalizedPositionsContainer_Centroid> (newPositionsContainer));
	
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
	CountInt indices[MAX_DIMENSIONS];
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
#endif // #ifdef WITH_IGOR

ImagePtr LocalizedPositionsContainer_Multiplication::getLocalizedPositionsAsMatrix() const {
	size_t nPositions = this->positionsVector.size();
	
	ImagePtr matrix(new Image(static_cast<int>(nPositions), 5));	// magic number
	matrix->setConstant(0.0);
	
	for (size_t i = 0; i < nPositions; ++i) {
		(*matrix)(i, 0) = this->positionsVector.at(i).frameNumber;
		(*matrix)(i, 1) = this->positionsVector.at(i).width;
		(*matrix)(i, 2) = this->positionsVector.at(i).xPosition;
		(*matrix)(i, 3) = this->positionsVector.at(i).yPosition;
		(*matrix)(i, 4) = this->positionsVector.at(i).nFramesPresent;
	}
	
	return matrix;
}

void LocalizedPositionsContainer_Multiplication::writePositionsToFile(std::string filePath, std::string header) const {
	std::ofstream outputFile;
	size_t nPositions = this->positionsVector.size();
	if (nPositions == 0) {
		throw std::runtime_error("No positions localized!");
	}
	
	outputFile.open(filePath.c_str(), std::ios::trunc);
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

void LocalizedPositionsContainer_Multiplication::addPosition(std::shared_ptr<LocalizedPosition> newPosition) {
	// check if the type of positions that we are adding is suitable
	if (newPosition->getPositionType() != this->getPositionsType())
		throw std::runtime_error("Trying to append a position of a different type to a LocalizedPositionsContainer_Multiplication");
	
	// cast the pointer to the more specific type
	std::shared_ptr<LocalizedPosition_Multiplication> newPosition_Multiplication(std::static_pointer_cast<LocalizedPosition_Multiplication> (newPosition));
	
	this->positionsVector.push_back(*newPosition_Multiplication);
}

void LocalizedPositionsContainer_Multiplication::addPositions(std::shared_ptr<LocalizedPositionsContainer> newPositionsContainer) {
	// are we trying to add the same container to itself?
	if (this == newPositionsContainer.get()) {
		throw std::runtime_error("Trying to append a LocalizedPositionsContainer_Multiplication to itself");
	}
	
	// check if the positions container is of the right type
	if (newPositionsContainer->getPositionsType() != LOCALIZED_POSITIONS_TYPE_MULTIPLICATION) {
		throw std::runtime_error("Trying to append a container of a different type to a LocalizedPositionsContainer_Multiplication");
	}
	
	// cast the pointer to the more specific type
	std::shared_ptr<LocalizedPositionsContainer_Multiplication> newPositionsContainer_2DGauss(std::static_pointer_cast<LocalizedPositionsContainer_Multiplication> (newPositionsContainer));
	
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
	CountInt indices[MAX_DIMENSIONS];
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
#endif // #ifdef WITH_IGOR

ImagePtr LocalizedPositionsContainer_ZeissPALM::getLocalizedPositionsAsMatrix() const {
	size_t nPositions = this->positionsVector.size();
	
	ImagePtr matrix(new Image(static_cast<int>(nPositions), 6));	// magic number
	matrix->setConstant(0.0);
	
	for (size_t i = 0; i < nPositions; ++i) {
		(*matrix)(i, 0) = this->positionsVector.at(i).frameNumber;
		(*matrix)(i, 1) = this->positionsVector.at(i).integral;
		(*matrix)(i, 2) = this->positionsVector.at(i).xPosition;
		(*matrix)(i, 3) = this->positionsVector.at(i).yPosition;
		(*matrix)(i, 4) = this->positionsVector.at(i).positionDeviation;
		(*matrix)(i, 5) = this->positionsVector.at(i).nFramesPresent;
	}
	
	return matrix;
}

void LocalizedPositionsContainer_ZeissPALM::addPosition(std::shared_ptr<LocalizedPosition> newPosition) {
	// check if the type of positions that we are adding is suitable
	if (newPosition->getPositionType() != this->getPositionsType())
		throw std::runtime_error("Trying to append a position of a different type to a LocalizedPositionsContainer_Multiplication");
	
	// cast the pointer to the more specific type
	std::shared_ptr<LocalizedPosition_ZeissPALM> newPosition_ZeissPALM(std::static_pointer_cast<LocalizedPosition_ZeissPALM> (newPosition));
	
	this->positionsVector.push_back(*newPosition_ZeissPALM);
}

void LocalizedPositionsContainer_ZeissPALM::addPositions(std::shared_ptr<LocalizedPositionsContainer> newPositionsContainer) {
	// are we trying to add the same container to itself?
	if (this == newPositionsContainer.get()) {
		throw std::runtime_error("Trying to append a LocalizedPositionsContainer_ZeissPALM to itself");
	}
	
	// check if the positions container is of the right type
	if (newPositionsContainer->getPositionsType() != LOCALIZED_POSITIONS_TYPE_ZEISSPALM) {
		throw std::runtime_error("Trying to append a container of a different type to a LocalizedPositionsContainer_ZeissPALM");
	}
	
	// cast the pointer to the more specific type
	std::shared_ptr<LocalizedPositionsContainer_ZeissPALM> newPositionsContainer_ZeissPALM(std::static_pointer_cast<LocalizedPositionsContainer_ZeissPALM> (newPositionsContainer));
	
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
	CountInt indices[MAX_DIMENSIONS];
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
#endif // #ifdef WITH_IGOR

ImagePtr LocalizedPositionsContainer_MLEwG::getLocalizedPositionsAsMatrix() const {
	size_t nPositions = this->positionsVector.size();
	
	ImagePtr matrix(new Image(static_cast<int>(nPositions), 8));	// magic number
	matrix->setConstant(0.0);
	
	for (size_t i = 0; i < nPositions; ++i) {
		(*matrix)(i, 0) = this->positionsVector.at(i).frameNumber;
		(*matrix)(i, 1) = this->positionsVector.at(i).integral;
		(*matrix)(i, 2) = this->positionsVector.at(i).width;
		(*matrix)(i, 3) = this->positionsVector.at(i).xPosition;
		(*matrix)(i, 4) = this->positionsVector.at(i).yPosition;
		(*matrix)(i, 5) = this->positionsVector.at(i).background;
		(*matrix)(i, 6) = this->positionsVector.at(i).positionDeviation;
		(*matrix)(i, 7) = this->positionsVector.at(i).nFramesPresent;
	}
	
	return matrix;
}

void LocalizedPositionsContainer_MLEwG::writePositionsToFile(std::string filePath, std::string header) const {
	std::ofstream outputFile;
	size_t nPositions = this->positionsVector.size();
	if (nPositions == 0) {
		throw std::runtime_error("No positions localized!");
	}
	
	outputFile.open(filePath.c_str(), std::ios::trunc);
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

void LocalizedPositionsContainer_MLEwG::addPosition(std::shared_ptr<LocalizedPosition> newPosition) {
	// check if the type of positions that we are adding is suitable
	if (newPosition->getPositionType() != this->getPositionsType())
		throw std::runtime_error("Trying to append a position of a different type to a LocalizedPositionsContainer_MLEwG");
	
	// cast the pointer to the more specific type
	std::shared_ptr<LocalizedPosition_MLEwG> newPosition_MLEwG(std::static_pointer_cast<LocalizedPosition_MLEwG> (newPosition));
	
	this->positionsVector.push_back(*newPosition_MLEwG);
}

void LocalizedPositionsContainer_MLEwG::addPositions(std::shared_ptr<LocalizedPositionsContainer> newPositionsContainer) {
	// are we trying to add the same container to itself?
	if (this == newPositionsContainer.get()) {
		throw std::runtime_error("Trying to append a LocalizedPositionsContainer_MLEwG to itself");
	}
	
	// check if the positions container is of the right type
	if (newPositionsContainer->getPositionsType() != this->getPositionsType()) {
		throw std::runtime_error("Trying to append a container of a different type to a LocalizedPositionsContainer_MLEwG");
	}
	
	// cast the pointer to the more specific type
	std::shared_ptr<LocalizedPositionsContainer_MLEwG> newPositionsContainer_MLEwG(std::static_pointer_cast<LocalizedPositionsContainer_MLEwG> (newPositionsContainer));
	
	for (std::vector<LocalizedPosition_MLEwG>::iterator it = newPositionsContainer_MLEwG->positionsVector.begin(); it != newPositionsContainer_MLEwG->positionsVector.end(); ++it) {
		this->positionsVector.push_back(*it);
	}
}
