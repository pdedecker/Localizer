/*
 Copyright 2008-2011 Peter Dedecker.
 
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

#include "PALM_analysis_ProgressReporting.h"
#include <math.h>

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

#if WITH_MATLAB
int ProgressReporter_MatlabCommandLine::UpdateCalculationProgress(double progress, double maxProgress) {
	double percentDone = progress / maxProgress * 100.0;
	
	// unfortunately there does not appear to be a simple way to detect a user abort in Matlab
	
	if (percentDone - previousPercentage > 5.0) {
		previousPercentage = floor(percentDone / 5.0) * 5.0;
		mexPrintf("%.0lf%% ", previousPercentage);
	}
	return 0;
}
#endif // WITH_MATLAB

