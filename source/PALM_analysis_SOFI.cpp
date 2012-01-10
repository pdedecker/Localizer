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

void DoSOFIAnalysis(boost::shared_ptr<ImageLoader> imageLoader, boost::shared_ptr<ImageOutputWriter> outputWriter,
                    boost::shared_ptr<ImageOutputWriter> averageImageOutputWriter,
					std::vector<boost::shared_ptr<SOFIFrameVerifier> > frameVerifiers, boost::shared_ptr<ProgressReporter> progressReporter,
					size_t nFramesToSkip, size_t nFramesToInclude, int lagTime, int order, int crossCorrelate, int nFramesToGroup) {
	size_t nImages = imageLoader->getNImages();
	size_t blockSize = 50;
    
    int doAverage = averageImageOutputWriter.get() != NULL;
    
    if (nFramesToInclude == (size_t)-1)
        nFramesToInclude = nImages - nFramesToSkip;
    
    if (nFramesToSkip + nFramesToInclude > nImages)
        throw std::runtime_error("Too many images requested");
	
	if (nImages <= lagTime)
		throw std::runtime_error("Not enough images for the requested lagtime");
	
	if (nFramesToSkip > nImages - lagTime)
		throw std::runtime_error("Too many frames are being skipped (/SKIP flag)");
	
	if ((nFramesToGroup <= lagTime) && (nFramesToGroup != 0))
		throw std::runtime_error("Cannot group less frames than the lag time");
    
    // adjust to the range requested by the user
    nImages = nFramesToSkip + nFramesToInclude;
	
	boost::shared_ptr<SOFICalculator> sofiCalculator;
	if (crossCorrelate == 0) {
		sofiCalculator = boost::shared_ptr<SOFICalculator>(new SOFICalculator_AutoCorrelation(order, lagTime, blockSize));
	} else {
		sofiCalculator = boost::shared_ptr<SOFICalculator>(new SOFICalculator_CrossCorrelation(order, lagTime, blockSize));
	}
	
	size_t nImagesToProcess = nImages - nFramesToSkip;
	int nGroups;
	if (nFramesToGroup != 0) {
		nGroups = ceil((double)(nImagesToProcess) / (double)(nFramesToGroup));
	} else {
		nGroups = 1;
	}
	
	ImagePtr currentImage;
	imageLoader->spoolTo(nFramesToSkip);
	
	size_t nFramesInThisGroup, nFramesIncludedInCalculation;
	size_t nFramesProcessedInThisRun;
	size_t nFramesProcessedTotal = 0;
	ImagePtr outputImage;
    ImagePtr averageOutputImage;
	ImagePtr sofiImage;
    ImagePtr averageImage;
	size_t nImagesInBlock;
	int status;
	int isValidFrame;
	
	double weightOfThisBlock, sumOfBlockWeights;
	
	progressReporter->CalculationStarted();
	
	for (size_t currentGroup = 0; currentGroup < nGroups; currentGroup += 1) {
		if (nFramesToGroup > 0) {
			nFramesInThisGroup = std::min(nFramesToGroup, (int)(nImagesToProcess - (currentGroup * nFramesToGroup)));
		} else {
			nFramesInThisGroup = nImagesToProcess;
		}
		
		nFramesProcessedInThisRun = 0;
		nImagesInBlock = 0;
		sumOfBlockWeights = 0.0;
		while (nFramesProcessedInThisRun < nFramesInThisGroup) {
			nImagesInBlock = std::min(nFramesInThisGroup - nFramesProcessedInThisRun, blockSize);
			if (nImagesInBlock < lagTime + 1)
				break;
			
            nFramesIncludedInCalculation = 0;
			for (size_t i = 0; i < nImagesInBlock; ++i) {
				currentImage = imageLoader->readNextImage();
				
				// verify if this is an acceptable image
				isValidFrame = 1;
				for (size_t j = 0; j < frameVerifiers.size(); j+=1) {
					if (frameVerifiers[j]->isValidFrame(currentImage) != 1) {
						isValidFrame = 0;
						break;
					}
				}
				if (isValidFrame == 1) {
					sofiCalculator->addNewImage(currentImage);
                    nFramesIncludedInCalculation += 1;
				}
				nFramesProcessedInThisRun += 1;
				nFramesProcessedTotal += 1;
			}
			
            try {
                sofiCalculator->getResult(outputImage, averageOutputImage);
            }
            catch (SOFINoImageInCalculation e) {
                // insufficient images were included in this calculation to get a result
                // this probably means that the verifier rejected (nearly) all of them
                // continue with the other images and hope that they make up for it
                continue;
            }
            
            status = progressReporter->UpdateCalculationProgress(nFramesProcessedTotal, nImagesToProcess);
			if (status != 0) {
				progressReporter->CalculationAborted();
				throw USER_ABORTED("Abort requested by user");
			}
            
			weightOfThisBlock = static_cast<double>(nFramesIncludedInCalculation) / static_cast<double>(blockSize);
			sumOfBlockWeights += weightOfThisBlock;
				
			if (sofiImage.get() == NULL) {
				sofiImage = outputImage;
                *sofiImage *= weightOfThisBlock;
			} else {
				*sofiImage += *outputImage * weightOfThisBlock;
			}
            if (averageImage.get() == NULL) {
                averageImage = averageOutputImage;
                *averageImage *= weightOfThisBlock;
            } else {
                *averageImage += *averageOutputImage * weightOfThisBlock;
            }
		}
		
		if (sofiImage.get() != NULL) {
			*sofiImage /= sumOfBlockWeights;
			
			outputWriter->write_image(sofiImage);
            
            if (doAverage) {
                *averageImage /= sumOfBlockWeights;
                averageImageOutputWriter->write_image(averageImage);
            }
		}
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
    
    if (imageVector.size() == batchSize) {
        // it's time to calculate and store some output
        double thisBatchWeight = static_cast<double>(imageVector.size()) / batchSize;
        ImagePtr thisBatchAverage, thisBatchSOFI;
        this->performCalculation(thisBatchSOFI, thisBatchAverage);
        
        // set up the output images if required
        if (this->sofiImage.get() == NULL) {
            this->sofiImage = thisBatchSOFI;
        } else {
            *(this->sofiImage) += *thisBatchSOFI * thisBatchWeight;
        }
        
        if (this->averageImage.get() == NULL) {
            this->averageImage = thisBatchAverage;
        } else {
            *(this->averageImage) += *thisBatchAverage * thisBatchWeight;
        }
        
        sumOfWeightsIncluded += thisBatchWeight;
    }
}

void SOFICalculator::getResult(ImagePtr &calculatedSOFIImage, ImagePtr &calculatedAverageImage) {
    // first check if the output exists
    if ((this->averageImage.get() == NULL) || (this->sofiImage.get() == NULL))
        throw std::runtime_error("Requested a SOFI image even though there are no evaluations");
    
    // calculate the averaged SOFI and average image
    calculatedSOFIImage = ImagePtr(new Image(*sofiImage));
    calculatedAverageImage = ImagePtr(new Image(*averageImage));
    
    *calculatedSOFIImage /= sumOfWeightsIncluded;
    *calculatedAverageImage /= sumOfWeightsIncluded;
}

SOFICalculator_AutoCorrelation::SOFICalculator_AutoCorrelation(int order_rhs, int lagTime_rhs, size_t batchSize) :
    SOFICalculator(batchSize),
    order(order_rhs)
{
}

void SOFICalculator_AutoCorrelation::performCalculation(ImagePtr &calculatedSOFIImage, ImagePtr &calculatedAverageImage) {
	if (this->imageVector.size() == 0) {
		throw SOFINoImageInCalculation("Requested a SOFI image even though there are no evaluations");
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
	
	// since the lagtime is assumed to be zero for now, simply multiply the images together
    ImagePtr outputImage(new Image(firstImage->rows(), firstImage->cols()));
    ImagePtr subImage(new Image(firstImage->rows(), firstImage->cols()));
	
    outputImage->setConstant(0.0);
    for (int n = 0; n < this->imageVector.size(); ++n) {
        subImage->setConstant(1.0);
        for (int i = 0; i < order; ++i) {
            (*subImage).array() *= (*imageVector[n]).array();
        }
        *outputImage += *subImage;
    }
    
    *outputImage /= static_cast<double>(this->imageVector.size());
    
    calculatedSOFIImage = ImagePtr(new Image(*outputImage));
    calculatedAverageImage = ImagePtr(new Image(*thisBatchAverageImage));
	
	// now reset everything for the next calculation
	// before returning
	this->imageVector.clear();
}

SOFICalculator_CrossCorrelation::SOFICalculator_CrossCorrelation(int order_rhs, int lagTime_rhs, size_t batchSize) :
    SOFICalculator(batchSize),
    order(order_rhs)
{
}

void SOFICalculator_CrossCorrelation::performCalculation(ImagePtr &calculatedSOFIImage, ImagePtr &calculatedAverageImage) {
    // verify that there are enough images in the input buffer
    if (this->imageVector.size() == 0)
        throw SOFINoImageInCalculation("no images for SOFICalculator_CrossCorrelation::getResult()");
    
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
    
    // set everything up
    XCSOFIKernelProvider kernelProvider;
    int firstRow, firstCol, lastRow, lastCol;
    int nRowsOutput, nColsOutput;
    int outputRow, outputCol;
    int nKernelRows, nKernelCols;
    
    kernelProvider.getSizeOfOutputImage(order, firstImage->rows(), firstImage->cols(), nRowsOutput, nColsOutput);
    kernelProvider.getCoordinatesOfFirstUsableInputPixel(order, firstImage->rows(), firstImage->cols(), firstRow, firstCol);
    kernelProvider.getCoordinatesOfLastUsableInputPixel(order, firstImage->rows(), firstImage->cols(), lastRow, lastCol);
    boost::shared_array<std::vector<SOFIPixelCalculation> > kernel = kernelProvider.getKernel(order, 0, nKernelRows, nKernelCols);
    int nPixelsInKernel = nKernelRows * nKernelCols;
    
    // allocate the output image
    ImagePtr outputImage(new Image(nRowsOutput, nColsOutput));
    
    double currentVal, summedVal;
    ImagePtr currentImage;
    // main calculation
    // loop over all images
    for (int n = 0; n < this->imageVector.size(); ++n) {
        currentImage = imageVector[n];
        // loop over all pixels in the input
        for (int i = firstRow; i <= lastRow; ++i) {
            for (int j = firstCol; j <= lastCol; ++j) {
                // loop over all the kernel
                for (int k = 0; k < nPixelsInKernel; ++k) {
                    // loop over all calculations within this kernel pixel
                    std::vector<SOFIPixelCalculation> *calculationsForThisPixel = &(kernel[k]);
                    summedVal = 0;
                    for (int calculationIndex = 0; calculationIndex < calculationsForThisPixel->size(); ++calculationIndex) {
                        // and finally, loop over the different operations that each input requires
                        SOFIPixelCalculation *thisCalculation = &((*calculationsForThisPixel)[calculationIndex]);
                        currentVal = 1;
                        for (int multiplicationIndex = 0; multiplicationIndex < thisCalculation->inputRowDeltas.size(); ++multiplicationIndex) {
                            currentVal *= (*currentImage)(i + thisCalculation->inputRowDeltas[multiplicationIndex], j + thisCalculation->inputColDeltas[multiplicationIndex]);
                        }
                        summedVal += currentVal;
                    }
                    summedVal /= static_cast<double>(calculationsForThisPixel->size());
                    (*calculationsForThisPixel)[0].getOutputPixelCoordinates(order, i, j, outputRow, outputCol);
                    (*outputImage)(outputRow + (*calculationsForThisPixel)[0].outputRowDelta, outputCol + (*calculationsForThisPixel)[0].outputColDelta) += summedVal;
                }
            }
        }
    }
    
    // take the average
    (*outputImage) /= static_cast<double>(this->imageVector.size());
    if ((*outputImage).maxCoeff() > 1e50)
        XOPNotice("Too large coeff in outputImage\r");
    
    calculatedSOFIImage = ImagePtr(new Image(*outputImage));
    calculatedAverageImage = ImagePtr(new Image(*thisBatchAverageImage));
    
    this->imageVector.clear();
    return;
                                       
	// now normalize the pixel by requiring that the mean of every kind of pixel is the same
    /*int nKindsOfPixels = nPixelsInKernel;
    int kindOfRow, kindOfCol;
    Eigen::MatrixXd pixelAverages(nKernelRows, nKernelCols);
    pixelAverages.setConstant(0.0);
    int nPixelsOfEachKind = nRowsOutput * nColsOutput / nKindsOfPixels;
    for (int i = 0; i < nRowsOutput; ++i) {
        for (int j = 0; j < nColsOutput; ++j) {
            kindOfRow = i % nKernelRows;
            kindOfCol = j % nKernelCols;
            pixelAverages(kindOfRow, kindOfCol) += (*outputImage)(i, j);
        }
    }
    
    pixelAverages /= static_cast<double>(nPixelsOfEachKind);
    
    // now determine correction factors, arbitrary to the first element
    for (int i = 0; i < nKernelRows; ++i) {
        for (int j = 0; j < nKernelCols; ++j) {
            if (i == 0 && j == 0)
                continue;
            pixelAverages(i, j) = pixelAverages(i, j) / pixelAverages(0, 0);
        }
    }
    pixelAverages(0, 0) = 1;
    
    // allocate the corrected image
    ImagePtr correctedImage(new Image(*outputImage));
    // and correct it
    for (int i = 0; i < nRowsOutput; ++i) {
        for (int j = 0; j < nColsOutput; ++j) {
            kindOfRow = i % nKernelRows;
            kindOfCol = j % nKernelCols;
            (*correctedImage)(i, j) /= pixelAverages(kindOfRow, kindOfCol);
        }
    }
    
    this->imageVector.clear();
    this->averageImage->setConstant(0.0);
	return correctedImage;*/
}

void SOFIPixelCalculation::getOutputPixelCoordinates(int order, int inputRow, int inputCol, int &outputRow, int &outputCol) {
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
            // todo: error
            break;
    }
}

void XCSOFIKernelProvider::getCoordinatesOfFirstUsableInputPixel(int order, int nRowsInput, int nColsInput, int &firstRow, int &firstCol) {
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
            // TODO: error
            break;
    }
}

void XCSOFIKernelProvider::getCoordinatesOfLastUsableInputPixel(int order, int nRowsInput, int nColsInput, int &lastRow, int &lastCol) {
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
            // TODO: error
            break;
    }
}

void XCSOFIKernelProvider::getSizeOfOutputImage(int order, int nRowsInput, int nColsInput, int &nRows, int &nCols) {
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
            // todo: ERROR
            break;
    }
}

boost::shared_array<std::vector<SOFIPixelCalculation> > XCSOFIKernelProvider::getKernel(int order, int lagTime, int &nRows, int &nCols) {
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
            sofiPixelCalculation.imageIndices.resize(2);
            sofiPixelCalculation.imageIndices[0] = 0;
            sofiPixelCalculation.imageIndices[1] = 0;
            
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
            sofiPixelCalculation.imageIndices.resize(3);
            sofiPixelCalculation.imageIndices[0] = 0;
            sofiPixelCalculation.imageIndices[1] = 0;
            sofiPixelCalculation.imageIndices[2] = 0;
            
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
            // TODO error
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

/*ImagePtr SOFICorrector_Order2::doImageCorrection(ImagePtr imageToCorrect) {
	// estimate the PSF standard deviation
	double optimalPSFStdDev = determinePSFStdDev(imageToCorrect);
	
	ImagePtr correctedImage = performPSFCorrection(imageToCorrect.get(), optimalPSFStdDev);
	return correctedImage;
}

double SOFICorrector_Order2::determinePSFStdDev(ImagePtr imageToCorrect) {
	const gsl_min_fminimizer_type * minizerType = gsl_min_fminimizer_brent;
	gsl_min_fminimizer *minimizer = gsl_min_fminimizer_alloc(minizerType);
	if (minimizer == NULL) {
		throw std::bad_alloc();
	}
	
	gsl_function func;
	func.function = SOFICorrector_Order2::functionToMinimize;
	func.params = (void *)(&(*imageToCorrect));
	
	int status;
	
	// try to initialize the minimizer
	// somewhat annoyingly, the solver refuses to do anything
	// unless the boundaries and starting value indicate a minimum
	// so if we don't get lucky right away then search some more
	// for a starting guess
	double initialPSFWidth = 2.0;	// TODO: provide this initial value in a more sensible way
	status = gsl_min_fminimizer_set(minimizer, &func, initialPSFWidth, 0.1, 10.0);
	if (status == GSL_EINVAL) {
		double updatedInitialGuess = 0.1 + 0.5;
		do {
			status = gsl_min_fminimizer_set(minimizer, &func, updatedInitialGuess, 0.1, 10);
			updatedInitialGuess += 0.5;
		} while ((status == GSL_EINVAL) && (updatedInitialGuess < 10.0));
	}
	
	if (status == GSL_EINVAL) {
		// we didn't find any acceptable interval
		throw std::runtime_error("Unable to find acceptable starting guess for the PSF correction in 2nd order XC");
	}
	
	int maxIterations = 100, iterations = 0;
	for (;;) {
		status = gsl_min_fminimizer_iterate(minimizer);
		if (status != GSL_SUCCESS) {
			gsl_min_fminimizer_free(minimizer);
			throw std::runtime_error("gsl_min_fminimizer_iterate reported an error");
		}
		
		status = gsl_min_test_interval(gsl_min_fminimizer_x_lower(minimizer), gsl_min_fminimizer_x_upper(minimizer), 0.001, 0.0);
		if (status == GSL_SUCCESS)
			break;
		
		iterations += 1;
		if (iterations > maxIterations)
			break;
	}
	
	double minimum = gsl_min_fminimizer_minimum(minimizer);
	gsl_min_fminimizer_free(minimizer);
	
	return minimum;
}

ImagePtr SOFICorrector_Order2::performPSFCorrection(Image *image, double psfStdDev) {
	size_t nRows = image->rows();
	size_t nCols = image->cols();
	
	double autoPixelFactor = exp(- sqrt(2.0) / (2.0 * psfStdDev * psfStdDev));
	double horizontalFactor = exp(- (1.0 / 2.0) / (2.0 * psfStdDev * psfStdDev));
	double diagonalFactor = exp(- 1.0 / (2.0 * psfStdDev * psfStdDev));	// the 1.0 comes from sqrt(2.0) / sqrt(2.0)
	
	ImagePtr correctedImage(new Image(*image));
	
	// only loop over the crosscorrelation pixels
	for (size_t j = 0; j < nCols; ++j) {
		for (size_t i = 0; i < nRows; ++i) {
			if ((i % 2 == 0) && (j % 2 == 0)) {
				// autocorrelation pixel
				(*correctedImage)(i, j) /= autoPixelFactor;
				continue;
			}
			if ((i % 2 == 1) && (j % 2 == 1)) {
				// this is a diagonal crosscorrelation pixel
				(*correctedImage)(i, j) /= diagonalFactor;
			} else {
				(*correctedImage)(i, j) /= horizontalFactor;
			}
		}
	}
	
	return correctedImage;
}

double SOFICorrector_Order2::functionToMinimize(double psfStdDev, void *params) {
	Image* image = (Image *)params;
	
	size_t nRows = image->rows();
	size_t nCols = image->cols();
	
	ImagePtr correctedImage = performPSFCorrection(image, psfStdDev);
	
	// calculate the relative variance of the mean
	double avg, avgOfSquares, variance;
	
	avg = (*image).sum() / (double)(nRows * nCols);
	avgOfSquares = (*image).array().square().sum() / (double)(nRows * nCols);
	variance = avgOfSquares - avg * avg;
	
	return variance / (avg * avg);
}*/

SOFIFrameVerifier_NoSaturation::SOFIFrameVerifier_NoSaturation(int storageType_rhs) {
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
		case STORAGE_TYPE_INT4:
		case STORAGE_TYPE_UINT4:
		case STORAGE_TYPE_FP32:
		case STORAGE_TYPE_FP64:
			throw std::runtime_error("invalid storage type for saturation verification");
			break;
		default:
			throw std::runtime_error("unknown storage type for saturation verification");
			break;
	}
}

int SOFIFrameVerifier_NoSaturation::isValidFrame(ImagePtr frame) {
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

