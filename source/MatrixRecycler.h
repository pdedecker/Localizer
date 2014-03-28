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

#ifndef PALM_ANALYSIS_MATRIXRECYCLER
#define PALM_ANALYSIS_MATRIXRECYCLER

#include <list>
#include "boost/thread.hpp"
#include <eigen3/Eigen/Eigen>
#include "Storage.h"

const size_t kMaxUnusedMatrices = 20;

class MatrixRecycler {
public:
	MatrixRecycler() {;}
	~MatrixRecycler();
	
	Image *getMatrix(int nRows, int nCols);
	void freeMatrix(Image *matrixToFree);
	
	// free all matrices allocated by this matrix
	// if this function is called while some memory
	// allocated by this class is still in use then
	// an error will be thrown
	void freeAllMatrices();
	
protected:
	std::list<Image *> unusedMatrixList;
	std::list<Image *> usedMatrixList;
	
	boost::mutex recyclingMutex;
};

// the functions below act as proxies for the functions in the MatrixRecycler class,
// except that they act on a single, global instance

/**
 * Obtain a matrix of the requested dimensions from the globalMatrixRecycler
 */
Image* GetRecycledMatrix(int nRows, int nCols);

/**
 * Mark a matrix from the globalMatrixRecycler as no longer in use
 */
void FreeRecycledMatrix(Image* matrixToFree);

/**
 * Request that all reserved memory held in the recycler be freed
 */
void FreeAllRecycledMatrices();

#endif
