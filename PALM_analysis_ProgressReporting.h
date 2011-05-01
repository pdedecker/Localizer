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
class PALMAnalysisProgressReporter {
public:
	PALMAnalysisProgressReporter() {;}
	virtual ~PALMAnalysisProgressReporter() {;}
	
	virtual void CalculationStarted() = 0;
	virtual int UpdateCalculationProgress(double progress, double maxProgress) = 0;
	virtual void CalculationDone() = 0;
	virtual void CalculationAborted() = 0;
};

class PALMAnalysisProgressReporter_Silent : public PALMAnalysisProgressReporter {
public:
	PALMAnalysisProgressReporter_Silent() {;}
	~PALMAnalysisProgressReporter_Silent() {;}
	
	void CalculationStarted() {;}
	int UpdateCalculationProgress(double progress, double maxProgress) {return 0;}
	void CalculationDone() {;}
	void CalculationAborted() {;}
	
protected:
	double previousPercentage;
};

#ifdef WITH_IGOR
/**
 * @brief Print updates on the calculation progress to the Igor command line
 */
class PALMAnalysisProgressReporter_IgorCommandLine : public PALMAnalysisProgressReporter {
public:
	PALMAnalysisProgressReporter_IgorCommandLine() {;}
	~PALMAnalysisProgressReporter_IgorCommandLine() {;}
	
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

class PALMAnalysisProgressReporter_IgorUserFunction : public PALMAnalysisProgressReporter {
public:
	PALMAnalysisProgressReporter_IgorUserFunction(FUNCREF igorProgressFunction);
	~PALMAnalysisProgressReporter_IgorUserFunction() {;}
	
	void CalculationStarted() {;}
	int UpdateCalculationProgress(double progress, double maxProgress);
	void CalculationDone() {;}
	void CalculationAborted() {;}
	
protected:
	FunctionInfo igorProgressFunction;
};
#endif // WITH_IGOR

class PALMAnalysisProgressReporter_stdout : public PALMAnalysisProgressReporter {
public:
	PALMAnalysisProgressReporter_stdout() {;}
	~PALMAnalysisProgressReporter_stdout() {;}
	
	void CalculationStarted() {std::cout << "Running calculation... "; std::cout.flush();}
	int UpdateCalculationProgress(double progress, double maxProgress) {printf("\rRunning calculation... %4.1f%%", progress / maxProgress * 100.0); std::cout.flush(); return 0;}
	void CalculationDone() {std::cout << "\rRunning calculation... Finished!\n"; std::cout.flush();}
	void CalculationAborted() {std::cout << " Abort requested by user\n"; std::cout.flush();}
};

#endif
