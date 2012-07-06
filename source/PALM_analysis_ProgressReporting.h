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

#ifndef PALM_ANALYSIS_PROGRESSREPORTING_H
#define PALM_ANALYSIS_PROGRESSREPORTING_H

#include <iostream>

#ifdef WITH_IGOR
#include "XOPStandardHeaders.h"
#endif

#ifdef WITH_MATLAB
#include "mex.h"
#include "WaitMex.h"
#endif


/**
 * @brief An abstract base class that provides progress updates and allows user aborts
 */
class ProgressReporter {
public:
	ProgressReporter() {;}
	virtual ~ProgressReporter() {;}
	
	virtual void CalculationStarted() = 0;
	virtual int UpdateCalculationProgress(double progress, double maxProgress) = 0;
	virtual void CalculationDone() = 0;
	virtual void CalculationAborted() = 0;
};

class ProgressReporter_Silent : public ProgressReporter {
public:
	ProgressReporter_Silent() {;}
	~ProgressReporter_Silent() {;}
	
	void CalculationStarted() {;}
	int UpdateCalculationProgress(double progress, double maxProgress);
	void CalculationDone() {;}
	void CalculationAborted() {;}
	
protected:
	double previousPercentage;
};

#ifdef WITH_IGOR
/**
 * @brief Print updates on the calculation progress to the Igor command line
 */
class ProgressReporter_IgorCommandLine : public ProgressReporter {
public:
	ProgressReporter_IgorCommandLine() {;}
	~ProgressReporter_IgorCommandLine() {;}
	
	void CalculationStarted() {previousPercentage = 0; XOPNotice("Running calculation... ");}
	int UpdateCalculationProgress(double progress, double maxProgress);
	void CalculationDone() {XOPNotice("Calculation finished!\r");}
	void CalculationAborted() {XOPNotice("Abort requested by user\r");}
	
protected:
	double previousPercentage;
};

#pragma pack(2) 
struct IgorUserFunctionParams {
	double progress;
	double maxProgress;
};
typedef struct IgorUserFunctionParams IgorUserFunctionParams;
#pragma pack()

class ProgressReporter_IgorUserFunction : public ProgressReporter {
public:
	ProgressReporter_IgorUserFunction(FUNCREF igorProgressFunction);
	~ProgressReporter_IgorUserFunction() {;}
	
	void CalculationStarted() {;}
	int UpdateCalculationProgress(double progress, double maxProgress);
	void CalculationDone() {;}
	void CalculationAborted() {;}
	
protected:
	FunctionInfo igorProgressFunction;
};
#endif // WITH_IGOR

class ProgressReporter_stdout : public ProgressReporter {
public:
	ProgressReporter_stdout() {;}
	~ProgressReporter_stdout() {;}
	
	void CalculationStarted() {std::cout << "Running calculation... "; std::cout.flush();}
	int UpdateCalculationProgress(double progress, double maxProgress) {printf("\rRunning calculation... %4.1f%%", progress / maxProgress * 100.0); std::cout.flush(); return 0;}
	void CalculationDone() {std::cout << "\rRunning calculation... Finished!\n"; std::cout.flush();}
	void CalculationAborted() {std::cout << " Abort requested by user\n"; std::cout.flush();}
};

#ifdef WITH_MATLAB
class ProgressReporter_MatlabCommandLine : public ProgressReporter {
public:
	ProgressReporter_MatlabCommandLine() {;}
	~ProgressReporter_MatlabCommandLine() {;}
	
	void CalculationStarted() {previousPercentage = 0; mexPrintf("Running calculation...");}
	int UpdateCalculationProgress(double progress, double maxProgress);
	void CalculationDone() {mexPrintf("Calculation finished!\n");}
	void CalculationAborted() {mexPrintf("Abort requested by user\r");}
	
protected:
	double previousPercentage;
};

class ProgressReporter_MatlabWaitMex : public ProgressReporter {
	public:
	ProgressReporter_MatlabWaitMex() : previousPercentage(0), waitBar(NULL) {}
	~ProgressReporter_MatlabWaitMex();
	
	void CalculationStarted();
	int UpdateCalculationProgress(double progress, double maxProgress);
	void CalculationDone();
	void CalculationAborted() {this->CalculationDone();}
	
protected:
	double previousPercentage;
	waitbar* waitBar;
};
#endif

#endif
