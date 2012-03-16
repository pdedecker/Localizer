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
 with analysis programs such as Igor Pro or Matlab, 
 the licensors of this Program grant you additional permission 
 to convey the resulting work.
 */

#ifndef PALM_ANALYSIS_CONVOLVER_H
#define PALM_ANALYSIS_CONVOLVER_H

#include <fftw3.h>
#include "boost/thread/mutex.hpp"
#include "boost/thread/shared_mutex.hpp"

#include "PALM_analysis_defines.h"
#include "PALM_analysis_MatrixRecycler.h"

class ConvolveMatricesWithFFTClass {
public:
	ConvolveMatricesWithFFTClass() {;}
	~ConvolveMatricesWithFFTClass();
	
	ImagePtr ConvolveMatrixWithSmallKernel(ImagePtr image, ImagePtr kernel);
	ImagePtr ConvolveMatricesWithFFT(ImagePtr image1, ImagePtr image2);
	ImagePtr ConvolveMatrixWithGivenFFT(ImagePtr image, boost::shared_ptr<fftw_complex> array2_FFT, size_t FFT_xSize2, size_t FFT_ySize2);
	
	// http://www.leptonica.com/convolution.html
	ImagePtr ConvolveMatrixWithFlatKernel(ImagePtr image, size_t kernelXSize, size_t kernelYSize);
	
	/**
	 * Get the FFT of a single image, possibly for later use in a convolution.
	 * The x and y size of the transformed image will be returned by reference in FFT_xSize and FFT_ySize.
	 */
	boost::shared_ptr<fftw_complex> DoForwardFFT(ImagePtr image);
	
	/**
	 * Calculate the reverse FFT
	 */
	ImagePtr DoReverseFFT(boost::shared_ptr<fftw_complex> array_FFT, size_t xSize, size_t ySize);
protected:
	static boost::mutex FFTWPlannerMutex;
};

#endif
