/*
 *  PALM_analysis_ProgressReporting.cpp
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 30/04/2011.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "PALM_analysis_ProgressReporting.h"

int ProgressReporter_Silent::UpdateCalculationProgress(double progress, double maxProgress) {
#ifdef WITH_IGOR
	int abortStatus = CheckAbort(0);
	if (abortStatus != 0)
		return abortStatus;
#endif
	return 0;
}

#ifdef WITH_IGOR
int ProgressReporter_IgorCommandLine::UpdateCalculationProgress(double progress, double maxProgress) {
	double percentDone = progress / maxProgress * 100.0;
	char XOPOut[10];
	
	// check if the user wants to abort
	int abortStatus = CheckAbort(0);
	if (abortStatus != 0)
		return abortStatus;
	
	if (percentDone - previousPercentage > 10.0) {
		previousPercentage = floor(percentDone / 10.0) * 10.0;
		sprintf(XOPOut, "%.0lf%% ", previousPercentage);
		XOPNotice(XOPOut);
	}
	return 0;
}

ProgressReporter_IgorUserFunction::ProgressReporter_IgorUserFunction(FUNCREF igorProgressFunction) {
	int err;
	FunctionInfo fi;
	int requiredParameterTypes[2];
	int badParameterNumber;
	
	// Make sure the function exists and get information about it
	err = GetFunctionInfoFromFuncRef(igorProgressFunction, &fi);
	if (err != 0)
		throw err;
	
	// Make sure the function has the right form
	requiredParameterTypes[0] = NT_FP64;
	requiredParameterTypes[1] = NT_FP64;
	
	err = CheckFunctionForm(&fi, 2, requiredParameterTypes, &badParameterNumber, NT_FP64);
	if (err != 0)
		throw err;
	
	// the function is valid, we're all set
	this->igorProgressFunction = fi;
}

int ProgressReporter_IgorUserFunction::UpdateCalculationProgress(double progress, double maxProgress) {
	// check if the user wants to abort
	int abortStatus = CheckAbort(0);
	if (abortStatus != 0)
		return abortStatus;
	
	// call the progress function
	double result;
	int err;
	IgorUserFunctionParams params;
	params.progress = progress;
	params.maxProgress = maxProgress;
	
	err = CallFunction(&(this->igorProgressFunction), (void *)&params, &result);
	if (err != 0)
		throw err;
	
	return (int)(result + 0.5);
}

#endif // WITH_IGOR
