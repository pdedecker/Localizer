#ifndef NEWSOFI_H
#define NEWSOFI_H

#include <vector>
#include <memory>

#include "PALM_analysis_defines.h"
#include "PALM_analysis_FileIO.h"
#include "PALM_analysis_ProgressReporting.h"
#include "NewSOFIKernels.h"

void DoNewSOFI(std::shared_ptr<ImageLoader> imageLoader, std::shared_ptr<ProgressReporter> progressReporter, const int order, std::vector<ImagePtr>& sofiOutputImages);

void DoNewSOFI2(std::shared_ptr<ImageLoader> imageLoader, std::shared_ptr<ProgressReporter> progressReporter, const int order, std::vector<ImagePtr>& sofiOutputImages);

void RawSOFIWorker(std::shared_ptr<ImageLoader> imageLoader, const int firstImageToProcess, const int lastImageToProcess, int &imagesProcessedSoFar, const int& totalNumberOfImagesToProcess, std::shared_ptr<ProgressReporter> progressReporter, const std::vector<std::pair<int, std::vector<SOFIKernel> > >& orders, std::map<PixelCombination,ImagePtr,ComparePixelCombinations>& pixelMap, std::vector<ImagePtr>& sofiImages);

#endif
