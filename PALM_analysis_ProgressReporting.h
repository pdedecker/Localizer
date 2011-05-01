/*
 *  PALM_analysis_ProgressReporting.h
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 30/04/2011.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef PALM_ANALYSIS_PROGRESSREPORTING_H
#define PALM_ANALYSIS_PROGRESSREPORTING_H

#include <iostream>

#ifdef WITH_IGOR
#include "XOPStandardHeaders.h"
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

#endif
