/*
 Copyright 2008-2016 Peter Dedecker.
 
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

#include "SOFICommandLine.h"

#include <iostream>

#include "FileIO.h"
#include "NewSOFI.h"
#include "ProgressReporting.h"
#include "Defines.h"


int main(int argc, char** argv) {
    // need to have 1 argument: path to the data file
    if (argc != 2) {
        std::cerr << "need a file path as argument";
        return 1;
    }
    
    std::string filePath(argv[1]);
    
    try {
        std::shared_ptr<ImageLoader> imageLoader = GetImageLoader(filePath);
        std::shared_ptr<ProgressReporter> progressReporter(new ProgressReporter_stdout());
        SOFIOptions options;
        std::vector<ImagePtr> sofiImages;
        std::string outputFilePath = filePath.substr(0, filePath.length() - 4) + "_sofi.tif";   // simple-minded and dangerous
        
        DoNewSOFI(imageLoader, options, progressReporter, sofiImages);
        LocalizerTIFFImageOutputWriter outputWriter(outputFilePath, 1, false, kFP32);
        outputWriter.write_image(sofiImages.at(0));
    }
    catch (std::runtime_error e) {
        std::cerr << "caught runtime error: " << e.what() << std::endl;
        return 1;
    }
    catch (...) {
        std::cerr << "caught unknown exception" << std::endl;
        return 1;
    }
    
    return 0;
}
