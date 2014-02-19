#ifndef NEWSOFI_H
#define NEWSOFI_H

#include <vector>
#include <memory>

#include "PALM_analysis_defines.h"
#include "PALM_analysis_FileIO.h"
#include "PALM_analysis_ProgressReporting.h"
#include "SOFIFrameVerifiers.h"

void DoNewSOFI(std::shared_ptr<ImageLoader> imageLoader, const std::vector<std::shared_ptr<SOFIFrameVerifier> >& frameVerifiers, std::shared_ptr<ProgressReporter> progressReporter, const int order, std::vector<ImagePtr>& sofiOutputImages, bool wantAverageImage, ImagePtr& averageImage);

#endif
