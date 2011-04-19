/*
 *  PALM_analysis_MatrixRecycler.h
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 28/02/2011.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <list>
#include "boost/thread.hpp"
#include <Eigen/Eigen>

const size_t kMaxUnusedMatrices = 20;

class MatrixRecycler {
public:
	MatrixRecycler() {;}
	~MatrixRecycler();
	
	Eigen::MatrixXd *getMatrix(size_t nRows, size_t nCols);
	void freeMatrix(Eigen::MatrixXd *matrixToFree);
	
	// free all matrices allocated by this matrix
	// if this function is called while some memory
	// allocated by this class is still in use then
	// an error will be thrown
	void freeAllMatrices();
	
protected:
	std::list<Eigen::MatrixXd *> unusedMatrixList;
	std::list<Eigen::MatrixXd *> usedMatrixList;
	
	boost::mutex recyclingMutex;
};

// the functions below act as proxies for the functions in the MatrixRecycler class,
// except that they act on a single, global instance

/**
 * Obtain a matrix of the requested dimensions from the globalMatrixRecycler
 */
Eigen::MatrixXd* GetRecycledMatrix(size_t nRows, size_t nCols);

/**
 * Mark a matrix from the globalMatrixRecycler as no longer in use
 */
void FreeRecycledMatrix(Eigen::MatrixXd* matrixToFree);

/**
 * Request that all reserved memory held in the recycler be freed
 */
void FreeAllRecycledMatrices();
