/*
 *  PALM_analysis_MatrixRecycler.cpp
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 28/02/2011.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "PALM_analysis_MatrixRecycler.h"

MatrixRecycler::~MatrixRecycler() {
	// time to clean up unused memory
	// if there is still memory being reported in use
	// then throw an error
	if (this->usedMatrixList.size() > 0)
		throw std::runtime_error("MatrixRecycler still has memory reported in use");
	
	// delete all entries in the unusedMatrixList
	for (std::list<Eigen::MatrixXd*>::iterator it = this->usedMatrixList.begin(); it != this->usedMatrixList.end(); ++it) {
		delete (*it);
	}
}

Eigen::MatrixXd* MatrixRecycler::getMatrix(size_t nRows, size_t nCols) {
	boost::lock_guard<boost::mutex> locker(this->recyclingMutex);
	
	// try to find an unused matrix that meets the requirements
	
	for (std::list<Eigen::MatrixXd*>::iterator it = this->unusedMatrixList.begin(); it != this->unusedMatrixList.end(); ++it) {
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
	
	Eigen::MatrixXd* newMatrix = new Eigen::MatrixXd((int)nRows, (int)nCols);
	this->usedMatrixList.push_front(newMatrix);
	return newMatrix;
}

void MatrixRecycler::freeMatrix(Eigen::MatrixXd *matrixToFree) {
	boost::lock_guard<boost::mutex> locker(this->recyclingMutex);
	
	// matrixToFree is not required anymore for now
	// remove it from the used list
	// and copy it to the unused one
	// but only if it can be found in the list!
	for (std::list<Eigen::MatrixXd*>::iterator it = this->usedMatrixList.begin(); it != this->usedMatrixList.end(); ++it) {
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
