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

#include "PALM_analysis_SOFI.h"

#include "tbb/tbb.h"
#include "tbb/spin_mutex.h"

void DoSOFIAnalysis(std::shared_ptr<ImageLoader> imageLoader, std::vector<std::shared_ptr<SOFIFrameVerifier> > frameVerifiers, 
					std::shared_ptr<ProgressReporter> progressReporter,
					int nFramesToSkip, int nFramesToInclude, std::vector<int> lagTimes, int order, int crossCorrelate, int nFramesToGroup, 
					std::vector<ImagePtr>& sofiOutputImages, std::vector<ImagePtr>& averageOutputImages) {
	if (lagTimes.empty() && (order > 0))
		lagTimes.resize(order - 1, 0);
	int largestLagTime = *(std::max_element(lagTimes.begin(), lagTimes.end()));
	if (largestLagTime > 100)
		throw std::runtime_error("Lag times larger than 100 are not supported");
	size_t blockSize = std::max(50, 2 * largestLagTime);
	size_t nImages = imageLoader->getNImages();
	sofiOutputImages.clear();
	averageOutputImages.clear();
	
	if (nFramesToSkip < 0)
		nFramesToSkip = 0;
    
    if (nFramesToInclude <= 0)
        nFramesToInclude = nImages - nFramesToSkip;
    
    if (nFramesToSkip + nFramesToInclude > nImages)
        throw std::runtime_error("Too many images requested");
	
	if ((nFramesToInclude <= largestLagTime) || (largestLagTime > blockSize))
		throw std::runtime_error("Lag time too long");
	
	if (nFramesToSkip > nImages - largestLagTime)
		throw std::runtime_error("Too many frames are being skipped (/SKIP flag)");
	
	if ((nFramesToGroup <= largestLagTime) && (nFramesToGroup != 0))
		throw std::runtime_error("Cannot group less frames than the lag time");
    
    // adjust to the range requested by the user
    nImages = nFramesToSkip + nFramesToInclude;
	
	std::shared_ptr<SOFICalculator> sofiCalculator;
	if (crossCorrelate == 0) {
		sofiCalculator = std::shared_ptr<SOFICalculator>(new SOFICalculator_AutoCorrelation(order, lagTimes, blockSize));
	} else {
		sofiCalculator = std::shared_ptr<SOFICalculator>(new SOFICalculator_CrossCorrelation(order, lagTimes, blockSize));
	}
	
	size_t nImagesToProcess = nImages - nFramesToSkip;
	size_t nGroups;
	if (nFramesToGroup != 0) {
		nGroups = ceil((double)(nImagesToProcess) / (double)(nFramesToGroup));
	} else {
		nGroups = 1;
	}
	
	ImagePtr currentImage;
	imageLoader->spoolTo(nFramesToSkip);
	
	size_t nFramesInThisGroup;
	size_t nFramesProcessedTotal = 0;
	ImagePtr sofiImage;
    ImagePtr averageImage;
	int status, spinProcessStatus;
	int isValidFrame;
	
	progressReporter->CalculationStarted();
	
	for (size_t currentGroup = 0; currentGroup < nGroups; currentGroup += 1) {
		if (nFramesToGroup > 0) {
			nFramesInThisGroup = std::min(nFramesToGroup, (int)(nImagesToProcess - (currentGroup * nFramesToGroup)));
		} else {
			nFramesInThisGroup = nImagesToProcess;
		}
		
		for (size_t i = 0; i < nFramesInThisGroup; ++i) {
			currentImage = imageLoader->readNextImage();
			
			// verify if this is an acceptable image
			isValidFrame = 1;
			for (size_t j = 0; j < frameVerifiers.size(); j+=1) {
				if (frameVerifiers[j]->isValidFrame(currentImage) != 1) {
					isValidFrame = 0;
					break;
				}
			}
			
			sofiCalculator->addNewImage(currentImage);
			nFramesProcessedTotal += 1;
			if (nFramesProcessedTotal % 25 == 0) {
				status = progressReporter->UpdateCalculationProgress(nFramesProcessedTotal, nImagesToProcess);
				
				#ifdef WITH_IGOR
					spinProcessStatus = SpinProcess();
				#else
					spinProcessStatus = 0;
				#endif // WITH_IGOR
				
				if ((status != 0) || (spinProcessStatus != 0)) {
					progressReporter->CalculationAborted();
					throw USER_ABORTED("Abort requested by user");
				}
			}
		}
		
		sofiCalculator->getResult(sofiImage, averageImage);
		// safety check
		if ((sofiImage.get() == NULL) || (averageImage.get() == NULL))
			throw std::runtime_error("No output from sofiCalculator");
		
		sofiOutputImages.push_back(sofiImage);
		averageOutputImages.push_back(averageImage);
	}
	
	progressReporter->CalculationDone();
}

void SOFICalculator::addNewImage(ImagePtr newImage) {
	// a new image is available
	// add it to the queue, but if the queue reaches the batch size
    // then a result needs to be computed and stored
    
    if (newImage.get() == NULL)
        throw std::runtime_error("NULL image provided to SOFICalculator::addNewImage");
	
	this->imageVector.push_back(newImage);
    
    if (imageVector.size() == _batchSize) {
        // it's time to calculate and store some output
        addNewBlockFromStoredImages();
		imageVector.clear();
    }
}

void SOFICalculator::addNewBlockFromStoredImages() {
    int largestLagTime = *(std::max_element(_lagTimes.begin(), _lagTimes.end()));
    double weightOfThisBlock = static_cast<double>(imageVector.size() - largestLagTime) / static_cast<double>(_batchSize - largestLagTime);
    
    ImagePtr thisBatchAverage, thisBatchSOFI;
    this->performCalculation(thisBatchSOFI, thisBatchAverage);
    
    // set up the output images if required
    if (this->sofiImage.get() == NULL) {
        this->sofiImage = thisBatchSOFI;
    } else {
        *(this->sofiImage) += *thisBatchSOFI * weightOfThisBlock;
    }
    
    if (this->averageImage.get() == NULL) {
        this->averageImage = thisBatchAverage;
    } else {
        *(this->averageImage) += *thisBatchAverage * weightOfThisBlock;
    }
    
    sumOfWeightsIncluded += weightOfThisBlock;
}

void SOFICalculator::getResult(ImagePtr &calculatedSOFIImage, ImagePtr &calculatedAverageImage) {
    // process any images that may remain in the queue
    try {
        addNewBlockFromStoredImages();
    }
    catch (SOFINoImageInCalculation e) {
        // this is no problem, it just means that there were not enough images left
    }
    
    // now check whether the output images even exist
    if ((this->averageImage.get() == NULL) || (this->sofiImage.get() == NULL))
        throw std::runtime_error("Requested a SOFI image even though there are no evaluations");
    
    // calculate the averaged SOFI and average image
    calculatedSOFIImage = ImagePtr(new Image(*sofiImage));
    calculatedAverageImage = ImagePtr(new Image(*averageImage));
    
    *calculatedSOFIImage /= sumOfWeightsIncluded;
    *calculatedAverageImage /= sumOfWeightsIncluded;
    
    performCorrection(calculatedSOFIImage);
}

SOFICalculator_AutoCorrelation::SOFICalculator_AutoCorrelation(int order, std::vector<int> lagTimes, size_t batchSize) :
    SOFICalculator(order, lagTimes, batchSize)
{
	// for an order of n there should be exactly (n - 1) lag times
	if (_lagTimes.size() != order - 1)
		throw std::runtime_error("A SOFI image of order n requires (n - 1) lag times to be specified");
}

void SOFICalculator_AutoCorrelation::performCalculation(ImagePtr &calculatedSOFIImage, ImagePtr &calculatedAverageImage) {
	int largestLagTime = *(std::max_element(_lagTimes.begin(), _lagTimes.end()));
	if (this->imageVector.size() == 0) {
		throw SOFINoImageInCalculation("Requested a SOFI image even though there are no input images available");
    }
	if (this->imageVector.size() <= largestLagTime) {
		throw SOFINoImageInCalculation("Requested a SOFI image even though there are not enough input images available");
    }
    
    ImagePtr firstImage = this->imageVector.front();
	
	// calculate the average image
    ImagePtr thisBatchAverageImage(new Image(firstImage->rows(), firstImage->cols()));
    thisBatchAverageImage->setConstant(0.0);
    for (std::vector<ImagePtr>::iterator it = this->imageVector.begin(); it != imageVector.end(); ++it) {
        *thisBatchAverageImage += *(*it);
    }
    *thisBatchAverageImage /= static_cast<double>(this->imageVector.size());
    
    // and convert all the input images to fluctuation images by subtracting the average
    for (std::vector<ImagePtr>::iterator it = this->imageVector.begin(); it != imageVector.end(); ++it) {
        *(*it) -= *thisBatchAverageImage;
    }
	
	// determine how many evaluations we can make given the lag times
	int nEvaluations = this->imageVector.size() - largestLagTime;
	// set up a lag time vector that explicitly includes the zero-lag time for the first contribution,
	// so that all contributions can be treated symetrically.
	std::vector<int> expandedLagTimes = _lagTimes;
	expandedLagTimes.insert(expandedLagTimes.begin(), 0);
	
	// allocate support images
    ImagePtr outputImage(new Image(firstImage->rows(), firstImage->cols()));
    ImagePtr subImage(new Image(firstImage->rows(), firstImage->cols()));
	
    outputImage->setConstant(0.0);
    for (size_t n = 0; n < nEvaluations; ++n) {
        subImage->setConstant(1.0);
        for (int i = 0; i < _order; ++i) {
            (*subImage).array() *= (*imageVector.at(n + expandedLagTimes.at(i))).array();
        }
        *outputImage += *subImage;
    }
    
    *outputImage /= static_cast<double>(nEvaluations);
    
    calculatedSOFIImage = ImagePtr(new Image(*outputImage));
    calculatedAverageImage = ImagePtr(new Image(*thisBatchAverageImage));
}

SOFICalculator_CrossCorrelation::SOFICalculator_CrossCorrelation(int order, std::vector<int> lagTimes, size_t batchSize) :
    SOFICalculator(order, lagTimes, batchSize)
{
	// only allow order 2 or 3 for now
	if ((order != 2) && (order != 3))
		throw std::runtime_error("SOFI cross correlation currently supports only orders 2 or 3");
	
	// for an order of n there should be exactly (n - 1) lag times
	if (_lagTimes.size() != order - 1)
		throw std::runtime_error("A SOFI image of order n requires (n - 1) lag times to be specified");
}

void SOFICalculator_CrossCorrelation::performCalculation(ImagePtr &calculatedSOFIImage, ImagePtr &calculatedAverageImage) {
	int largestLagTime = *(std::max_element(_lagTimes.begin(), _lagTimes.end()));
    // verify that there are enough images in the input buffer
    if (this->imageVector.size() == 0)
        throw SOFINoImageInCalculation("no images for SOFICalculator_CrossCorrelation::getResult()");
	if (this->imageVector.size() <= largestLagTime) {
		throw SOFINoImageInCalculation("Requested a SOFI image even though there are not enough input images available");
    }
    
    ImagePtr firstImage = this->imageVector.front();
	
	// calculate the average image
    ImagePtr thisBatchAverageImage(new Image(firstImage->rows(), firstImage->cols()));
    thisBatchAverageImage->setConstant(0.0);
    for (std::vector<ImagePtr>::iterator it = this->imageVector.begin(); it != imageVector.end(); ++it) {
        *thisBatchAverageImage += *(*it);
    }
    *thisBatchAverageImage /= static_cast<double>(this->imageVector.size());
    
    // and convert all the input images to fluctuation images by subtracting the average
    for (std::vector<ImagePtr>::iterator it = this->imageVector.begin(); it != imageVector.end(); ++it) {
        *(*it) -= *thisBatchAverageImage;
    }
	
	// determine how many evaluations we can make given the lag times
	size_t nEvaluations = this->imageVector.size() - largestLagTime;
	// set up a lag time vector that explicitly includes the zero-lag time for the first contribution,
	// so that all contributions can be treated symetrically.
	std::vector<int> expandedLagTimes = _lagTimes;
	expandedLagTimes.insert(expandedLagTimes.begin(), 0);
    
    // set everything up
    XCSOFIKernelProvider kernelProvider;
    int firstRow, firstCol, lastRow, lastCol;
    int nRowsOutput, nColsOutput;
    int nKernelRows, nKernelCols;
    
    kernelProvider.getSizeOfOutputImage(_order, firstImage->rows(), firstImage->cols(), nRowsOutput, nColsOutput);
    kernelProvider.getCoordinatesOfFirstUsableInputPixel(_order, firstImage->rows(), firstImage->cols(), firstRow, firstCol);
    kernelProvider.getCoordinatesOfLastUsableInputPixel(_order, firstImage->rows(), firstImage->cols(), lastRow, lastCol);
    boost::shared_array<std::vector<SOFIPixelCalculation> > kernel = kernelProvider.getKernel(_order, 0, nKernelRows, nKernelCols);
    int nPixelsInKernel = nKernelRows * nKernelCols;
    
    // allocate the output image
    ImagePtr outputImage(new Image(nRowsOutput, nColsOutput));
	outputImage->setConstant(0.0);
    
    tbb::spin_mutex spinMutex;
    // main calculation
    // loop over all images
    tbb::parallel_for(tbb::blocked_range<size_t>(0, nEvaluations), [&](const tbb::blocked_range<size_t> &thisRange) {
        int outputRow, outputCol;
        for (size_t n = thisRange.begin(); n != thisRange.end(); ++n) {
            // loop over all pixels in the input
            for (int i = firstRow; i <= lastRow; ++i) {
                for (int j = firstCol; j <= lastCol; ++j) {
                    // loop over the kernel
                    for (int k = 0; k < nPixelsInKernel; ++k) {
                        // loop over all calculations within this kernel pixel
                        const std::vector<SOFIPixelCalculation> *calculationsForThisPixel = &(kernel[k]);
                        double summedVal = 0.0;
                        for (size_t calculationIndex = 0; calculationIndex < calculationsForThisPixel->size(); ++calculationIndex) {
                            // and finally, loop over the different operations that each input requires
                            const SOFIPixelCalculation *thisCalculation = &((*calculationsForThisPixel)[calculationIndex]);
                            double currentVal = 1.0;
                            for (size_t multiplicationIndex = 0; multiplicationIndex < thisCalculation->inputRowDeltas.size(); ++multiplicationIndex) {
                                currentVal *= (*imageVector.at(n + expandedLagTimes.at(multiplicationIndex)))(i + thisCalculation->inputRowDeltas[multiplicationIndex], j + thisCalculation->inputColDeltas[multiplicationIndex]);
                            }
                            summedVal += currentVal;
                        }
                        summedVal /= static_cast<double>(calculationsForThisPixel->size());
                        (*calculationsForThisPixel)[0].getOutputPixelCoordinates(_order, i, j, outputRow, outputCol);
                        {
                            tbb::spin_mutex::scoped_lock(spinMutex);
                            (*outputImage)(outputRow + (*calculationsForThisPixel)[0].outputRowDelta, outputCol + (*calculationsForThisPixel)[0].outputColDelta) += summedVal;
                        }
                    }
                }
            }
        }
    });
    
    // take the average
    (*outputImage) /= static_cast<double>(nEvaluations);
    
	calculatedSOFIImage = outputImage;
    calculatedAverageImage = thisBatchAverageImage;
    
    return;
}

void SOFICalculator_CrossCorrelation::performCorrection(ImagePtr imageToCorrect) {
    // Each type of cross-cumulant has a different scaling factor depending on the distance between the contributing
    // pixels. Originally we corrected for this by scaling the values so that every kind of pixel has the same mean
    // value. However, this did not fix the checkerboarding completely since there appears to be some kind of residual
    // correlation as an artifact of camera readout. So now we normalize by requiring that the mean and standard deviation
    // of each type of pixel are the same. This is achieved by transforming each type of pixel according to y = ax + b.
    int nRowsOutput = imageToCorrect->rows();
    int nColsOutput = imageToCorrect->cols();
    
    XCSOFIKernelProvider kernelProvider;
    int nKernelRows, nKernelCols;
    kernelProvider.getKernel(_order, 0, nKernelRows, nKernelCols);
    int nKindsOfPixels = nKernelRows * nKernelCols;
    
    int kindOfRow, kindOfCol;
    Eigen::MatrixXd pixelAverages(nKernelRows, nKernelCols);
    Eigen::MatrixXd pixelVariances(nKernelRows, nKernelCols);
    Eigen::MatrixXd aFactor(nKernelRows, nKernelCols), bTerm(nKernelRows, nKernelCols);
    pixelAverages.setConstant(0.0);
    pixelVariances.setConstant(0.0);
    int nPixelsOfEachKind = nRowsOutput * nColsOutput / nKindsOfPixels;
    
    // calculate averages
    for (int i = 0; i < nRowsOutput; ++i) {
        for (int j = 0; j < nColsOutput; ++j) {
            kindOfRow = i % nKernelRows;
            kindOfCol = j % nKernelCols;
            pixelAverages(kindOfRow, kindOfCol) += (*imageToCorrect)(i, j);
        }
    }
    pixelAverages /= static_cast<double>(nPixelsOfEachKind);
    
    // calculate variances
    for (int i = 0; i < nRowsOutput; ++i) {
        for (int j = 0; j < nColsOutput; ++j) {
            kindOfRow = i % nKernelRows;
            kindOfCol = j % nKernelCols;
            pixelVariances(kindOfRow, kindOfCol) += square<double>((*imageToCorrect)(i, j) - pixelAverages(kindOfRow, kindOfCol));
        }
    }
    pixelVariances /= static_cast<double>(nPixelsOfEachKind - 1);
    
    // now determine correction factors, arbitrary to the first element
    for (int i = 0; i < nKernelRows; ++i) {
        for (int j = 0; j < nKernelCols; ++j) {
            if (i == 0 && j == 0) {
                aFactor(i, j) = 1.0;
                bTerm(i, j) = 0.0;
                continue;
            }
            aFactor(i, j) = std::sqrt(pixelVariances(0, 0) / pixelVariances(i, j));
            bTerm(i, j) = pixelAverages(0, 0) - aFactor(i, j) * pixelAverages(i, j);
        }
    }
    
    // perform the correction
    for (int i = 0; i < nRowsOutput; ++i) {
        for (int j = 0; j < nColsOutput; ++j) {
            kindOfRow = i % nKernelRows;
            kindOfCol = j % nKernelCols;
            (*imageToCorrect)(i, j) = (*imageToCorrect)(i, j) * aFactor(kindOfRow, kindOfCol) + bTerm(kindOfRow, kindOfCol);
        }
    }
}

void SOFIPixelCalculation::getOutputPixelCoordinates(int order, int inputRow, int inputCol, int &outputRow, int &outputCol) const {
    switch (order) {
        case 2:
            outputRow = 2 * inputRow - 2;
            outputCol = 2 * inputCol - 2;
            break;
        case 3:
            outputRow = 3 * inputRow - 3;
            outputCol = 3 * inputCol - 3;
            break;
        default:
            throw std::runtime_error("Invalid case in getOutputPixelCoordinates");
            break;
    }
}

void XCSOFIKernelProvider::getCoordinatesOfFirstUsableInputPixel(int order, int nRowsInput, int nColsInput, int &firstRow, int &firstCol) const {
    switch (order) {
        case 2:
            firstRow = 1;
            firstCol = 1;
            break;
        case 3:
            firstRow = 1;
            firstCol = 1;
            break;
        default:
            throw std::runtime_error("Invalid case in getCoordinatesOfFirstUsableInputPixel");
            break;
    }
}

void XCSOFIKernelProvider::getCoordinatesOfLastUsableInputPixel(int order, int nRowsInput, int nColsInput, int &lastRow, int &lastCol) const {
    switch (order) {
        case 2:
            lastRow = (nRowsInput - 1) - 1;
            lastCol = (nColsInput - 1) - 1;
            break;
        case 3:
            lastRow = (nRowsInput - 1) - 1;
            lastCol = (nColsInput - 1) - 1;
            break;
        default:
            throw std::runtime_error("Invalid case in getCoordinatesOfLastUsableInputPixel");
            break;
    }
}

void XCSOFIKernelProvider::getSizeOfOutputImage(int order, int nRowsInput, int nColsInput, int &nRows, int &nCols) const {
    switch (order) {
        case 2:
            nRows = 2 * nRowsInput - 4;
            nCols = 2 * nColsInput - 4;
            break;
        case 3:
            nRows = 3 * nRowsInput - 6;
            nCols = 3 * nColsInput - 6;
            break;
        default:
            throw std::runtime_error("Invalid case in getSizeOfOutputImage");
            break;
    }
}

boost::shared_array<std::vector<SOFIPixelCalculation> > XCSOFIKernelProvider::getKernel(int order, int lagTime, int &nRows, int &nCols) const {
    boost::shared_array<std::vector<SOFIPixelCalculation> > kernel;
    
    switch (order) {
        case 2:
        {
            // kernel size is 2x2
            nRows = 2;
            nCols = 2;
            kernel = boost::shared_array<std::vector<SOFIPixelCalculation> >(new std::vector<SOFIPixelCalculation>[4]);
            SOFIPixelCalculation sofiPixelCalculation;
            // always two inputs for 2nd order
            sofiPixelCalculation.inputRowDeltas.resize(2);
            sofiPixelCalculation.inputColDeltas.resize(2);
            
            // the auto pixel at relative coordinates (0,0), expressed as a crosscorrelation
            sofiPixelCalculation.inputRowDeltas[0] = -1;
            sofiPixelCalculation.inputColDeltas[0] = 0;
            sofiPixelCalculation.inputRowDeltas[1] = +1;
            sofiPixelCalculation.inputColDeltas[1] = 0;
            sofiPixelCalculation.outputRowDelta = 0;
            sofiPixelCalculation.outputColDelta = 0;
            kernel[0].push_back(sofiPixelCalculation);
            sofiPixelCalculation.inputRowDeltas[0] = 0;
            sofiPixelCalculation.inputColDeltas[0] = -1;
            sofiPixelCalculation.inputRowDeltas[1] = 0;
            sofiPixelCalculation.inputColDeltas[1] = +1;
            sofiPixelCalculation.outputRowDelta = 0;
            sofiPixelCalculation.outputColDelta = 0;
            kernel[0].push_back(sofiPixelCalculation);
            
            // vertical (along cols) cross pixel
            sofiPixelCalculation.inputRowDeltas[0] = 0;
            sofiPixelCalculation.inputColDeltas[0] = 0;
            sofiPixelCalculation.inputRowDeltas[1] = +1;
            sofiPixelCalculation.inputColDeltas[1] = 0;
            sofiPixelCalculation.outputRowDelta = +1;
            sofiPixelCalculation.outputColDelta = 0;
            kernel[1].push_back(sofiPixelCalculation);
            
            // horizontal (along rows) cross pixel
            sofiPixelCalculation.inputRowDeltas[0] = 0;
            sofiPixelCalculation.inputColDeltas[0] = 0;
            sofiPixelCalculation.inputRowDeltas[1] = 0;
            sofiPixelCalculation.inputColDeltas[1] = +1;
            sofiPixelCalculation.outputRowDelta = 0;
            sofiPixelCalculation.outputColDelta = +1;
            kernel[2].push_back(sofiPixelCalculation);
            
            // diagonal cross pixel
            sofiPixelCalculation.inputRowDeltas[0] = +1;
            sofiPixelCalculation.inputColDeltas[0] = 0;
            sofiPixelCalculation.inputRowDeltas[1] = 0;
            sofiPixelCalculation.inputColDeltas[1] = +1;
            sofiPixelCalculation.outputRowDelta = +1;
            sofiPixelCalculation.outputColDelta = +1;
            kernel[3].push_back(sofiPixelCalculation);
            sofiPixelCalculation.inputRowDeltas[0] = 0;
            sofiPixelCalculation.inputColDeltas[0] = 0;
            sofiPixelCalculation.inputRowDeltas[1] = +1;
            sofiPixelCalculation.inputColDeltas[1] = +1;
            sofiPixelCalculation.outputRowDelta = +1;
            sofiPixelCalculation.outputColDelta = +1;
            kernel[3].push_back(sofiPixelCalculation);
            
            break;
        }
        case 3:
        {
            // kernel size is 3x3
            nRows = 3;
            nCols = 3;
            kernel = boost::shared_array<std::vector<SOFIPixelCalculation> >(new std::vector<SOFIPixelCalculation>[9]);
            SOFIPixelCalculation sofiPixelCalculation;
            // always 3 inputs for 3rd order
            sofiPixelCalculation.inputRowDeltas.resize(3);
            sofiPixelCalculation.inputColDeltas.resize(3);
            
            // the auto pixel at relative coordinates (0,0)
            sofiPixelCalculation.inputRowDeltas[0] = -1;
            sofiPixelCalculation.inputColDeltas[0] = 0;
            sofiPixelCalculation.inputRowDeltas[1] = 0;
            sofiPixelCalculation.inputColDeltas[1] = 0;
            sofiPixelCalculation.inputRowDeltas[2] = +1;
            sofiPixelCalculation.inputColDeltas[2] = 0;
            sofiPixelCalculation.outputRowDelta = 0;
            sofiPixelCalculation.outputColDelta = 0;
            kernel[0].push_back(sofiPixelCalculation);
            sofiPixelCalculation.inputRowDeltas[0] = 0;
            sofiPixelCalculation.inputColDeltas[0] = -1;
            sofiPixelCalculation.inputRowDeltas[1] = 0;
            sofiPixelCalculation.inputColDeltas[1] = 0;
            sofiPixelCalculation.inputRowDeltas[2] = 0;
            sofiPixelCalculation.inputColDeltas[2] = +1;
            sofiPixelCalculation.outputRowDelta = 0;
            sofiPixelCalculation.outputColDelta = 0;
            kernel[0].push_back(sofiPixelCalculation);
            
            // AAB
            sofiPixelCalculation.inputRowDeltas[0] = 0;
            sofiPixelCalculation.inputColDeltas[0] = 0;
            sofiPixelCalculation.inputRowDeltas[1] = 0;
            sofiPixelCalculation.inputColDeltas[1] = 0;
            sofiPixelCalculation.inputRowDeltas[2] = 0;
            sofiPixelCalculation.inputColDeltas[2] = +1;
            sofiPixelCalculation.outputRowDelta = 0;
            sofiPixelCalculation.outputColDelta = +1;
            kernel[1].push_back(sofiPixelCalculation);
            
            // ABB
            sofiPixelCalculation.inputRowDeltas[0] = 0;
            sofiPixelCalculation.inputColDeltas[0] = 0;
            sofiPixelCalculation.inputRowDeltas[1] = 0;
            sofiPixelCalculation.inputColDeltas[1] = +1;
            sofiPixelCalculation.inputRowDeltas[2] = 0;
            sofiPixelCalculation.inputColDeltas[2] = +1;
            sofiPixelCalculation.outputRowDelta = 0;
            sofiPixelCalculation.outputColDelta = +2;
            kernel[2].push_back(sofiPixelCalculation);
            
            // AAC
            sofiPixelCalculation.inputRowDeltas[0] = 0;
            sofiPixelCalculation.inputColDeltas[0] = 0;
            sofiPixelCalculation.inputRowDeltas[1] = 0;
            sofiPixelCalculation.inputColDeltas[1] = 0;
            sofiPixelCalculation.inputRowDeltas[2] = +1;
            sofiPixelCalculation.inputColDeltas[2] = 0;
            sofiPixelCalculation.outputRowDelta = +1;
            sofiPixelCalculation.outputColDelta = 0;
            kernel[3].push_back(sofiPixelCalculation);
            
            // ACC
            sofiPixelCalculation.inputRowDeltas[0] = 0;
            sofiPixelCalculation.inputColDeltas[0] = 0;
            sofiPixelCalculation.inputRowDeltas[1] = +1;
            sofiPixelCalculation.inputColDeltas[1] = 0;
            sofiPixelCalculation.inputRowDeltas[2] = +1;
            sofiPixelCalculation.inputColDeltas[2] = 0;
            sofiPixelCalculation.outputRowDelta = +2;
            sofiPixelCalculation.outputColDelta = 0;
            kernel[4].push_back(sofiPixelCalculation);
            
            // AAD
            sofiPixelCalculation.inputRowDeltas[0] = 0;
            sofiPixelCalculation.inputColDeltas[0] = 0;
            sofiPixelCalculation.inputRowDeltas[1] = 0;
            sofiPixelCalculation.inputColDeltas[1] = 0;
            sofiPixelCalculation.inputRowDeltas[2] = +1;
            sofiPixelCalculation.inputColDeltas[2] = +1;
            sofiPixelCalculation.outputRowDelta = +1;
            sofiPixelCalculation.outputColDelta = +1;
            kernel[5].push_back(sofiPixelCalculation);
            
            // ADD
            sofiPixelCalculation.inputRowDeltas[0] = 0;
            sofiPixelCalculation.inputColDeltas[0] = 0;
            sofiPixelCalculation.inputRowDeltas[1] = +1;
            sofiPixelCalculation.inputColDeltas[1] = +1;
            sofiPixelCalculation.inputRowDeltas[2] = +1;
            sofiPixelCalculation.inputColDeltas[2] = +1;
            sofiPixelCalculation.outputRowDelta = +2;
            sofiPixelCalculation.outputColDelta = +2;
            kernel[6].push_back(sofiPixelCalculation);
            
            // BCC
            sofiPixelCalculation.inputRowDeltas[0] = 0;
            sofiPixelCalculation.inputColDeltas[0] = +1;
            sofiPixelCalculation.inputRowDeltas[1] = +1;
            sofiPixelCalculation.inputColDeltas[1] = 0;
            sofiPixelCalculation.inputRowDeltas[2] = +1;
            sofiPixelCalculation.inputColDeltas[2] = 0;
            sofiPixelCalculation.outputRowDelta = +2;
            sofiPixelCalculation.outputColDelta = +1;
            kernel[7].push_back(sofiPixelCalculation);
            
            // BBC
            sofiPixelCalculation.inputRowDeltas[0] = 0;
            sofiPixelCalculation.inputColDeltas[0] = +1;
            sofiPixelCalculation.inputRowDeltas[1] = 0;
            sofiPixelCalculation.inputColDeltas[1] = +1;
            sofiPixelCalculation.inputRowDeltas[2] = +1;
            sofiPixelCalculation.inputColDeltas[2] = 0;
            sofiPixelCalculation.outputRowDelta = +1;
            sofiPixelCalculation.outputColDelta = +2;
            kernel[8].push_back(sofiPixelCalculation);
            
            break;
        }
        default:
            throw std::runtime_error("Invalid case in getKernel");
            break;
    }
    return kernel;
}

ImagePtr SOFICorrector_Order2::doImageCorrection(ImagePtr imageToCorrect) {
	size_t nRows = imageToCorrect->rows();
	size_t nCols = imageToCorrect->cols();
	
	double avgOfAuto = 0.0, avgOfHoriz = 0.0, avgOfDiag = 0.0;
	double nAuto = 0.0, nHoriz = 0.0, nDiag = 0.0;
	
	ImagePtr correctedImage(new Image(*imageToCorrect));
	
	for (size_t j = 0; j < nCols; j+=1) {
		for (size_t i = 0; i < nRows; i+=1) {
			if ((i % 2 == 0) && (j % 2 == 0)) {
				avgOfAuto += (*imageToCorrect)(i, j);
				nAuto += 1.0;
				continue;
			}
			if ((i % 2 == 1) && (j % 2 == 1)) {
				avgOfDiag += (*imageToCorrect)(i, j);
				nDiag += 1.0;
				continue;
			}
			
			avgOfHoriz += (*imageToCorrect)(i, j);
			nHoriz += 1.0;
		}
	}
	
	avgOfAuto /= nAuto;
	avgOfDiag /= nDiag;
	avgOfHoriz /= nHoriz;
	
	double horizFactor = avgOfHoriz / avgOfAuto;
	double diagFactor = avgOfDiag / avgOfAuto;
	
	for (size_t j = 0; j < nCols; j+=1) {
		for (size_t i = 0; i < nRows; i+=1) {
			if ((i % 2 == 0) && (j % 2 == 0)) {
				continue;
			}
			if ((i % 2 == 1) && (j % 2 == 1)) {
				(*correctedImage)(i, j) /= diagFactor;
				continue;
			}
			
			(*correctedImage)(i, j) /= horizFactor;
		}
	}
	
	return correctedImage;
}

SOFIFrameVerifier_NoSaturation::SOFIFrameVerifier_NoSaturation(int storageType_rhs) {
	
	validLimitsSet = true;
	switch (storageType_rhs) {
            // some systems (e.g. the Nikon TIRF) never report absolutely saturated pixel
            // values, but instead report a different value close to it
            // so add in some arbitrary margin for the saturation detection
		case STORAGE_TYPE_INT8:
			this->saturationValue = 127.0 - 1.0;
			break;
		case STORAGE_TYPE_UINT8:
			this->saturationValue = 255.0 - 2.0;
			break;
		case STORAGE_TYPE_INT16:
			this->saturationValue = 32767.0 - 200.0;
			break;
		case STORAGE_TYPE_UINT16:
			this->saturationValue = 65535.0 - 200.0;
			break;
		case STORAGE_TYPE_INT32:
			this->saturationValue = 2147483647.0 - 500.0;
			break;
		case STORAGE_TYPE_UINT32:
			this->saturationValue = 4294967295.0 - 500.0;
			break;
		case STORAGE_TYPE_INT64:
			this->saturationValue = 9223372036854775807.0 - 1000.0;
			break;
		case STORAGE_TYPE_UINT64:
			this->saturationValue = 18446744073709551615.0 - 1000.0;
			break;
		default:
			validLimitsSet = false;
			break;
	}
}

int SOFIFrameVerifier_NoSaturation::isValidFrame(ImagePtr frame) {
	if (validLimitsSet == false)
		return 1;
	
	size_t nRows = frame->rows();
	size_t nCols = frame->cols();
	
	size_t nPixels = nRows * nCols;
	double localSaturationValue = this->saturationValue;
	
	double *dataPtr = frame->data();
	for (size_t i = 0; i < nPixels; i+=1) {
		if (dataPtr[i] >= localSaturationValue) {
			return 0;
		}
	}
	
	return 1;
}

int SOFIFrameVerifier_MaxPixelValue::isValidFrame(ImagePtr image) {
	size_t nRows = image->rows();
	size_t nCols = image->cols();
	
	size_t nPixels = nRows * nCols;
	double localMaxPixelValue = this->maxPixelValue;
	
	double *dataPtr = image->data();
	for (size_t i = 0; i < nPixels; i+=1) {
		if (dataPtr[i] > localMaxPixelValue) {
			return 0;
		}
	}
	
	return 1;
}

