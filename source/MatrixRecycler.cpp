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

MatrixRecycler::~MatrixRecycler() {
	// time to clean up unused memory
	// if there is still memory being reported in use
	// then throw an error
	if (this->usedMatrixList.size() > 0)
		throw std::runtime_error("MatrixRecycler still has memory reported in use");
	
	// delete all entries in the unusedMatrixList
	for (std::list<Image*>::iterator it = this->unusedMatrixList.begin(); it != this->unusedMatrixList.end(); ++it) {
		delete (*it);
	}
}

Image* MatrixRecycler::getMatrix(int nRows, int nCols) {
	boost::lock_guard<boost::mutex> locker(this->recyclingMutex);
	
	// try to find an unused matrix that meets the requirements
	
	for (std::list<Image*>::iterator it = this->unusedMatrixList.begin(); it != this->unusedMatrixList.end(); ++it) {
		if (((*it)->rows() == nRows) && ((*it)->cols() == nCols)) {
			// this matrix is appropriate
			// copy it to the used list, remove it from this one, and return it
			this->usedMatrixList.push_front(*it);
			this->unusedMatrixList.erase(it);
			return this->usedMatrixList.front();
		}
	}
	
	// if we are still here then there was no unused matrix available
	// so allocate some new memory, add it to the used list, and return it
	
	Image* newMatrix = new Image((int)nRows, (int)nCols);
	this->usedMatrixList.push_front(newMatrix);
	return newMatrix;
}

void MatrixRecycler::freeMatrix(Image *matrixToFree) {
	boost::lock_guard<boost::mutex> locker(this->recyclingMutex);
	
	// matrixToFree is not required anymore for now
	// remove it from the used list
	// and move it to the unused one
	// but only if it can be found in the list!
	
	// first check to see how many unused matrices are available
	// if there are too many then delete the last one
	if (this->unusedMatrixList.size() >= kMaxUnusedMatrices) {
		delete this->unusedMatrixList.back();
		this->unusedMatrixList.pop_back();
	}
	
	for (std::list<Image*>::iterator it = this->usedMatrixList.begin(); it != this->usedMatrixList.end(); ++it) {
		if ((*it) == matrixToFree) {
			this->unusedMatrixList.push_front(*it);
			this->usedMatrixList.erase(it);
			return;
		}
	}
	
	// if we are still here then the matrix was not in the list
	// this is an error
	throw std::runtime_error("no matrix found in freeMatrix()");
}

void MatrixRecycler::freeAllMatrices() {
	boost::lock_guard<boost::mutex> locker(this->recyclingMutex);
	
	if (this->usedMatrixList.size() != 0) {
		throw std::runtime_error("tried to free all matrices while some were still in use.");
	}
	
	for (std::list<Image*>::iterator it = this->unusedMatrixList.begin(); it != this->unusedMatrixList.end(); ++it) {
		delete (*it);
	}
	
	unusedMatrixList.clear();
}

/**
 * A global instance of MatrixRecycler to be used in the
 * segmentation
 */
std::shared_ptr<MatrixRecycler> globalMatrixRecycler(new MatrixRecycler);

Image* GetRecycledMatrix(int nRows, int nCols) {
	return globalMatrixRecycler->getMatrix(nRows, nCols);
}

void FreeRecycledMatrix(Image* matrixToFree) {
	globalMatrixRecycler->freeMatrix(matrixToFree);
}

void FreeAllRecycledMatrices() {
	globalMatrixRecycler->freeAllMatrices();
}
