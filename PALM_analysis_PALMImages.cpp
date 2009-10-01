/*
 *  PALM_analysis_PALMImages.cpp
 *  PALM analysis XOP
 *
 *  Created by Peter Dedecker on 01/10/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "PALM_analysis_PALMImages.h"


boost::shared_ptr<PALMVolume <unsigned short> > calculate_PALM_bitmap_image(boost::shared_ptr<PALMMatrix<double> > positions, boost::shared_ptr<PALMMatrix<double> > colors, boost::shared_ptr<PALMBitmapImageDeviationCalculator> deviationCalculator,
																			size_t xSize, size_t ySize, size_t imageWidth, size_t imageHeight, int normalizeColors) {
	size_t nPositions = positions->getXSize();
	size_t nColors = colors->getXSize();
	size_t nFrames = (size_t)(positions->get(nPositions - 1, 0)) + 1;
	size_t currentFrame;
	
	double imageScaleFactor = (double)(imageWidth - 1) / (double)xSize;
	double maxAmplitude = 0;
	
	double currentX, currentY, currentAmplitude;
	size_t centerX, centerY;
	long startX, endX, startY, endY;
	double deviation, currentIntensity;
	double distanceXSquared, distanceYSquared;
	size_t colorIndex;
	double currentColors[3];
	double summedIntensity;
	
	boost::shared_ptr<PALMVolume <unsigned short> > outputImage(new PALMVolume <unsigned short>(imageWidth, imageHeight, 3));	// 3 layers because it will be a direct color image
	boost::shared_ptr<PALMMatrix<double> > totalIntensities(new PALMMatrix<double>(imageWidth, imageHeight));	// keep track of the total intensities in each pixel
	
	outputImage->set_all(0);
	totalIntensities->set_all(0);
	
	if (normalizeColors != 0) {
		// get the position with the maximum amplitude
		// those positions will have the 'full' colors, the colors of the other positions will be scaled relative to it
		for (size_t n = 0; n < nPositions; ++n) {
			if ((*positions)(n, 1) > maxAmplitude) {
				maxAmplitude = positions->get(n,1);
			}
		}
	}
	
	for (size_t n = 0; n < nPositions; ++n) {
		currentFrame = (size_t)((*positions)(n, 0) + 0.5);
		currentAmplitude = positions->get(n, 1);
		currentX = positions->get(n, 3);
		currentY = positions->get(n, 4);
		
		if ((currentAmplitude < 0) || (currentX < 0) || (currentX >= xSize) || (currentY < 0) || (currentY >= ySize)) {
			continue;
		}
		
		if (normalizeColors == 0) {
			currentAmplitude = 1.0;	// every position is equally important when we don't do scaling
		}
		
		centerX = (size_t)(currentX * imageScaleFactor + 0.5);
		centerY = (size_t)(currentY * imageScaleFactor + 0.5);
		deviation = (deviationCalculator->getDeviation(positions, n) * imageScaleFactor);
		
		startX = floor((double)centerX - 4.0 * deviation);
		startY = floor((double)centerY - 4.0 * deviation);
		endX = ceil((double)centerX + 4.0 * deviation);
		endY = ceil((double)centerY + 4.0 * deviation);
		
		if (startX < 0)
			startX = 0;
		if (endX >= imageWidth)
			endX = imageWidth - 1;
		if (startY < 0)
			startY = 0;
		if (endY >= imageHeight)
			endY = imageHeight - 1;
		
		colorIndex = (size_t)((double)currentFrame / (double)nFrames * (double)(nColors - 1) + 0.5);
		
		for (size_t i = startX; i <= endX; ++i) {
			for (size_t j = startY; j < endY; ++j) {
				distanceXSquared = ((double)i - (double)centerX) * ((double)i - (double)centerX);
				distanceYSquared = ((double)j - (double)centerY) * ((double)j - (double)centerY);
				currentIntensity = currentAmplitude * exp(- (distanceXSquared + distanceYSquared) / (2 * deviation * deviation));
				
				if (normalizeColors != 0) {
					currentColors[0] = (*colors)(colorIndex, 0) * currentIntensity / maxAmplitude;	// Simplification of colors->get(colorIndex, 0) * currentIntensity / currentAmplitude * currentAmplitude / maxAmplitude
					currentColors[1] = (*colors)(colorIndex, 1) * currentIntensity / maxAmplitude;
					currentColors[2] = (*colors)(colorIndex, 2) * currentIntensity / maxAmplitude;
				} else {
					currentColors[0] = (*colors)(colorIndex, 0) * currentIntensity / currentAmplitude;
					currentColors[1] = (*colors)(colorIndex, 1) * currentIntensity / currentAmplitude;
					currentColors[2] = (*colors)(colorIndex, 2) * currentIntensity / currentAmplitude;
				}
				
				summedIntensity = currentIntensity + (*totalIntensities)(i, j);
				
				(*outputImage)(i, j, 0) = (unsigned short)(currentIntensity / summedIntensity * currentColors[0] + totalIntensities->get(i, j) / summedIntensity * outputImage->get(i, j, 0));
				(*outputImage)(i, j, 1) = (unsigned short)(currentIntensity / summedIntensity * currentColors[1] + totalIntensities->get(i, j) / summedIntensity * outputImage->get(i, j, 1));
				(*outputImage)(i, j, 2) = (unsigned short)(currentIntensity / summedIntensity * currentColors[2] + totalIntensities->get(i, j) / summedIntensity * outputImage->get(i, j, 2));
				
				(*totalIntensities)(i,j) = totalIntensities->get(i, j) + currentIntensity;
			}
		}
	}
	
	return outputImage;
}

boost::shared_ptr<PALMVolume <unsigned short> > calculate_PALM_bitmap_image_parallel(boost::shared_ptr<PALMMatrix<double> > positions, boost::shared_ptr<PALMMatrix<double> > colors, boost::shared_ptr<PALMBitmapImageDeviationCalculator> deviationCalculator,
																					 size_t xSize, size_t ySize, size_t imageWidth, size_t imageHeight, int normalizeColors) {
	vector<boost::shared_ptr<boost::thread> > threads;
	boost::shared_ptr<boost::thread> singleThreadPtr;
	vector<size_t> startPositions;
	vector<size_t> endPositions;
	vector<boost::shared_ptr<calculate_PALM_bitmap_image_ThreadStartParameters> > threadData;
	boost::shared_ptr<PALMVolume <unsigned short> > outputImage;
	boost::shared_ptr<PALMMatrix<double> > totalIntensities;
	
	size_t nPositions = positions->getXSize();
	size_t numberOfProcessors = boost::thread::hardware_concurrency();
	if (numberOfProcessors == 0) {
		numberOfProcessors = 1;
	}
	size_t nThreads = numberOfProcessors;
	size_t nPositionsPerThread;
	size_t currentThreadStart;
	size_t currentThreadEnd;
	size_t nFrames = (size_t)(positions->get(nPositions - 1, 0)) + 1;
	
	double maxAmplitude = 0;
	double imageScaleFactor = (double)(imageWidth - 1) / (double)xSize;
	double summedIntensity;
	
	nThreads = nPositions / 5;
	if (nThreads > numberOfProcessors) {
		nThreads = numberOfProcessors;
	}
	
	if (nThreads == 0) {
		nThreads = 1;
	}
	
	nPositionsPerThread = nPositions / nThreads;
	
	if (normalizeColors != 0) {
		// get the position with the maximum amplitude
		// those positions will have the 'full' colors, the colors of the other positions will be scaled relative to it
		for (size_t n = 0; n < nPositions; ++n) {
			if ((*positions)(n,1) > maxAmplitude) {
				maxAmplitude = (*positions)(n,1);
			}
		}
	}
	
	// set up the vector containing the data for the threads
	threadData.clear();
	currentThreadEnd = (size_t)-1;
	for (size_t i = 0; i < nThreads; ++i) {
		if (i == nThreads - 1) {	// the last thread is special
			currentThreadStart = currentThreadEnd + 1;
			currentThreadEnd = nPositions - 1;
		} else {
			currentThreadStart = currentThreadEnd + 1;
			currentThreadEnd = currentThreadStart + nPositionsPerThread;
		}
		
		boost::shared_ptr<calculate_PALM_bitmap_image_ThreadStartParameters> singleThreadStartParameter(new calculate_PALM_bitmap_image_ThreadStartParameters);
		singleThreadStartParameter->positions = positions;
		singleThreadStartParameter->image = boost::shared_ptr<PALMVolume <unsigned short> > (new PALMVolume <unsigned short>(imageWidth, imageHeight, 3));
		singleThreadStartParameter->totalIntensities = boost::shared_ptr<PALMMatrix<double> > (new PALMMatrix<double>(imageWidth, imageHeight));
		singleThreadStartParameter->colors = colors;
		singleThreadStartParameter->deviationCalculator = deviationCalculator;
		singleThreadStartParameter->normalizeColors = normalizeColors;
		singleThreadStartParameter->nFrames = nFrames;
		singleThreadStartParameter->startIndex = currentThreadStart;
		singleThreadStartParameter->endIndex = currentThreadEnd;
		singleThreadStartParameter->imageWidth = imageWidth;
		singleThreadStartParameter->imageHeight = imageHeight;
		singleThreadStartParameter->xSize = xSize;
		singleThreadStartParameter->ySize = ySize;
		singleThreadStartParameter->maxAmplitude = maxAmplitude;
		singleThreadStartParameter->scaleFactor = imageScaleFactor;
		
		threadData.push_back(singleThreadStartParameter);
	}
	
	// now start the threads
	threads.clear();
	for (size_t j = 0; j < nThreads; ++j) {
		singleThreadPtr = boost::shared_ptr<boost::thread>(new boost::thread(&calculate_PALM_bitmap_image_ThreadStart, threadData.at(j)));
		threads.push_back(singleThreadPtr);
	}
	
	// wait for the threads to finish
	for (size_t j = 0; j < nThreads; ++j) {
		threads.at(j)->join();
	}
	
	// combine the individual images into the final output image
	outputImage = threadData.at(0)->image;
	totalIntensities = threadData.at(0)->totalIntensities;
	
	for (size_t n = 1; n < nThreads; ++n) {
		for (size_t i = 0; i < imageWidth; ++i) {
			for (size_t j = 0; j < imageHeight; ++j) {
				summedIntensity = totalIntensities->get(i,j) + threadData[n]->totalIntensities->get(i,j);
				(*outputImage)(i, j, 0) = (unsigned short)((*totalIntensities)(i,j) / summedIntensity * (*outputImage)(i, j, 0) + threadData[n]->totalIntensities->get(i,j) / summedIntensity * threadData[n]->image->get(i, j, 0));
				(*outputImage)(i, j, 1) = (unsigned short)((*totalIntensities)(i,j) / summedIntensity * (*outputImage)(i, j, 1) + threadData[n]->totalIntensities->get(i,j) / summedIntensity * threadData[n]->image->get(i, j, 0));
				(*outputImage)(i, j, 2) = (unsigned short)((*totalIntensities)(i,j) / summedIntensity * (*outputImage)(i, j, 2) + threadData[n]->totalIntensities->get(i,j) / summedIntensity * threadData[n]->image->get(i, j, 0));
				(*totalIntensities)(i,j) = summedIntensity;
			}
		}
	}
	
	return outputImage;
}


void calculate_PALM_bitmap_image_ThreadStart(boost::shared_ptr<calculate_PALM_bitmap_image_ThreadStartParameters> startParameters) {
	boost::shared_ptr<PALMMatrix<double> > positions = startParameters->positions;
	boost::shared_ptr<PALMVolume <unsigned short> > image = startParameters->image;
	boost::shared_ptr<PALMMatrix<double> > totalIntensities = startParameters->totalIntensities;
	boost::shared_ptr<PALMMatrix<double> > colors = startParameters->colors;
	
	boost::shared_ptr<PALMBitmapImageDeviationCalculator> deviationCalculator = startParameters->deviationCalculator;
	
	size_t nColors = colors->getXSize();
	size_t startIndex = startParameters->startIndex;
	size_t endIndex = startParameters->endIndex;
	size_t nFrames = startParameters->nFrames;
	size_t currentFrame;
	int normalizeColors = startParameters->normalizeColors;
	
	size_t imageWidth = startParameters->imageWidth;
	size_t imageHeight = startParameters->imageHeight;
	size_t xSize = startParameters->xSize;
	size_t ySize = startParameters->ySize;
	double imageScaleFactor = startParameters->scaleFactor;
	double maxAmplitude = startParameters->maxAmplitude;
	
	double currentX, currentY, currentAmplitude;
	size_t centerX, centerY;
	long startX, endX, startY, endY;
	double deviation, currentIntensity;
	double distanceXSquared, distanceYSquared;
	size_t colorIndex;
	double currentColors[3];
	double summedIntensity;
	
	image->set_all(0);
	totalIntensities->set_all(0);
	
	for (size_t n = startIndex; n <= endIndex; ++n) {
		currentFrame = (size_t)((*positions)(n, 0) + 0.5);
		currentAmplitude = (*positions)(n, 1);
		currentX = (*positions)(n, 3);
		currentY = (*positions)(n, 4);
		
		if ((currentAmplitude < 0) || (currentX < 0) || (currentX >= xSize) || (currentY < 0) || (currentY >= ySize)) {
			continue;
		}
		
		if (normalizeColors == 0) {
			currentAmplitude = 1.0;	// every position is equally important when we don't do scaling
		}
		
		centerX = (size_t)(currentX * imageScaleFactor + 0.5);
		centerY = (size_t)(currentY * imageScaleFactor + 0.5);
		deviation = (deviationCalculator->getDeviation(positions, n) * imageScaleFactor);
		
		startX = floor((double)centerX - 4.0 * deviation);
		startY = floor((double)centerY - 4.0 * deviation);
		endX = ceil((double)centerX + 4.0 * deviation);
		endY = ceil((double)centerY + 4.0 * deviation);
		
		if (startX < 0)
			startX = 0;
		if (endX >= imageWidth)
			endX = imageWidth - 1;
		if (startY < 0)
			startY = 0;
		if (endY >= imageHeight)
			endY = imageHeight - 1;
		
		colorIndex = (size_t)((double)currentFrame / (double)nFrames * (double)(nColors - 1) + 0.5);
		
		for (size_t i = startX; i <= endX; ++i) {
			for (size_t j = startY; j < endY; ++j) {
				distanceXSquared = ((double)i - (double)centerX) * ((double)i - (double)centerX);
				distanceYSquared = ((double)j - (double)centerY) * ((double)j - (double)centerY);
				currentIntensity = currentAmplitude * exp(- (distanceXSquared + distanceYSquared) / (2 * deviation * deviation));
				
				if (normalizeColors != 0) {
					currentColors[0] = colors->get(colorIndex, 0) * currentIntensity / maxAmplitude;	// Simplification of colors->get(colorIndex, 0) * currentIntensity / currentAmplitude * currentAmplitude / maxAmplitude
					currentColors[1] = colors->get(colorIndex, 1) * currentIntensity / maxAmplitude;
					currentColors[2] = colors->get(colorIndex, 2) * currentIntensity / maxAmplitude;
				} else {
					currentColors[0] = colors->get(colorIndex, 0) * currentIntensity / currentAmplitude;
					currentColors[1] = colors->get(colorIndex, 1) * currentIntensity / currentAmplitude;
					currentColors[2] = colors->get(colorIndex, 2) * currentIntensity / currentAmplitude;
				}
				
				summedIntensity = currentIntensity + totalIntensities->get(i, j);
				
				image->set(i, j, 0, (unsigned short)(currentIntensity / summedIntensity * currentColors[0] + totalIntensities->get(i, j) / summedIntensity * image->get(i, j, 0)));
				image->set(i, j, 1, (unsigned short)(currentIntensity / summedIntensity * currentColors[1] + totalIntensities->get(i, j) / summedIntensity * image->get(i, j, 1)));
				image->set(i, j, 2, (unsigned short)(currentIntensity / summedIntensity * currentColors[2] + totalIntensities->get(i, j) / summedIntensity * image->get(i, j, 2)));
				
				totalIntensities->set(i,j, (totalIntensities->get(i, j) + currentIntensity));
			}
		}
	}
}
