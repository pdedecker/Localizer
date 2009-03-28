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
		"get_nth_image returned NULL",
		"The last image should have a larger or equal index to the starting image",
		"get_nth_image requested an image index beyond the number of images in the file",
		"get_nth_image reported that there is no file open to load from",
		"Invalid SPE storage type in get_nth_image()",
		"Cannot open the input file.",
		"Unable to create the output file.",
		"The size of a 'char' variable is not one byte, I cannot continue like this",
		"The size of a 'float' variable is not four bytes, I cannot continue like this",
		"The output file already exists. Use '/O' to overwrite",
		"The number of frames to average over should be an odd number",
		"Unknown method (/M flag)",
		"Unknown or unsupported CCD file type (/Y flag)",
		"Unknown method (/M flag)",
		"Unknown threshold preprocessing method (/G flag)",
		"Unknown threshold postprocessing method (/G flag)",
		"There was an error reading from the data file",
		"There was an error writing to the output file",
		
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
		"AnalyzePALMImages",
		XOPOp+utilOp+compilableOp,
		
		"ReadCCDImages",
		XOPOp+utilOp+compilableOp+threadSafeOp,
		
		"ProcessCCDImages",
		XOPOp+utilOp+compilableOp,
		
		"AnalyzeCCDImages",
		XOPOp+utilOp+compilableOp,
		
		"TestThreshold",
		XOPOp+utilOp+compilableOp,
		
		"ConvolveImages",
		XOPOp+utilOp+compilableOp,
		
		"MakeBitmapPALMImage",
		XOPOp+utilOp+compilableOp,
	}
};