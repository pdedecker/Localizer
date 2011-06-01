#include "XOPStandardHeaders.r"

resource 'vers' (1) {						// XOP version info.
	0x01, 0x20, release, 0x00, 0,			// version bytes and country integer.
	"1.20",
	"1.20, © 1993-2002 WaveMetrics, Inc., all rights reserved."
};

resource 'vers' (2) {						// Igor version info.
	0x05, 0x00, release, 0x00, 0,			// Version bytes and country integer.
	"5.00",
	"(for Igor 5.00 or later)"
};

resource 'STR#' (1100) { // Custom error messages 
	{ 
		"An error occurred while performing the requested operation. See the command window (command-J or control-J) for more information.",
		"A rolling average requires an odd number of frames to be specified (/AVG flag)",
	} 
};

// No menu item

resource 'XOPI' (1100) {
	XOP_VERSION,							// XOP protocol version.
	DEV_SYS_CODE,							// Development system information.
	0,										// Obsolete - set to zero.
	0,										// Obsolete - set to zero.
	XOP_TOOLKIT_VERSION,					// XOP Toolkit version.
};

resource 'XOPC' (1100) {
	{
		"LocalizationAnalysis",
		XOPOp+utilOp+compilableOp+threadSafeOp,
		
		"ReadCCDImages",
		XOPOp+utilOp+compilableOp+threadSafeOp,
		
		"ProcessCCDImages",
		XOPOp+utilOp+compilableOp+threadSafeOp,
		
		"AnalyzeCCDImages",
		XOPOp+utilOp+compilableOp+threadSafeOp,
		
		"EmitterSegmentation",
		XOPOp+utilOp+compilableOp+threadSafeOp,
		
		"ConvolveImages",
		XOPOp+utilOp+compilableOp,
		
		"LocalizationBitmap",
		XOPOp+utilOp+compilableOp+threadSafeOp,
		
		"RipleyLFunctionClustering",
		XOPOp+utilOp+compilableOp,
		
		"SOFIAnalysis",
		XOPOp+utilOp+compilableOp,
	}
};
