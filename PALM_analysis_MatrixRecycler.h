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
	
protected:
	std::list<Eigen::MatrixXd *> unusedMatrixList;
	std::list<Eigen::MatrixXd *> usedMatrixList;
	
	boost::mutex recyclingMutex;
};

/**
 * A function that will handle allocation of memory from globalMatrixRecycler
 */
Eigen::MatrixXd* GetRecycledMatrix(size_t nRows, size_t nCols);

/**
 * A function that will handle freeing of memory from globalMatrixRecycler
 */
void FreeRecycledMatrix(Eigen::MatrixXd* matrixToFree);
