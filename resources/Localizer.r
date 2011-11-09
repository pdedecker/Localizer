// Copyright 2008-2011 Peter Dedecker.

// This file is part of Localizer.

// Localizer is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// Localizer is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with Localizer.  If not, see <http://www.gnu.org/licenses/>.


// Additional permission under GNU GPL version 3 section 7

// If you modify this Program, or any covered work, by 
// linking or combining it with libraries required for interaction 
// with analysis programs such as Igor Pro or Matlab, or to acquire
// data from or control hardware related to an experimental measurement,
// the licensors of this Program grant you additional permission 
// to convey the resulting work.


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
