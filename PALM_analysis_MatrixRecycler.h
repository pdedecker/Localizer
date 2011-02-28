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
