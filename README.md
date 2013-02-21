README
======

Welcome to the source code for Localizer, which implements various kinds of analyses related to superresolution microscopy. The routines were lovingly hand-crafted by Peter Dedecker (<peter.dedecker@chem.kuleuven.be>), who currently resides at the [Department of Chemistry](http://www.chem.kuleuven.be/department/department_en.html) at the [University of Leuven](http://www.kuleuven.be). This repository is available at at [https://bitbucket.org/pdedecker/localizer](https://bitbucket.org/pdedecker/localizer).

A paper on Localizer is now available from the [Journal of Biomedical Optics](http://spie.org/x866.xml). [Download the paper here](http://sushi.chem.kuleuven.be/LocalizerJBO.pdf) after checking that you have access to the journal. Please cite this paper if Localizer turns out to be useful for your research.

*The fastest way to get started with Localizer is to [download the ready-made binaries appropriate for your platform](http://sushi.chem.kuleuven.be/svn/Localizer). This is all that you need to use the software*. If you have trouble accessing that link then try [this one](http://134.58.38.13/svn/Localizer).

At present there are binaries for [Igor Pro](http://www.wavemetrics.com) and for [Matlab](http://www.themathworks.com). There is also a command-line interface that allows Localizer to function in a commandline environment, which is mostly aimed at experienced users and will require compilation (see below).

If you are unfamiliar with programming the Igor Pro plugin will be the best option, since it comes with a complete graphical interface to use for analyzing your PALM/STORM or SOFI data, including visualization, clustering analysis, drift correction, etc. Igor Pro is commercial software, but [a free 30-day trial version is available for direct download](http://www.wavemetrics.com/support/demos.htm). Localizer will function without restrictions during the trial period. After the trial period Localizer will continue to function normally, but you will be unable to save the data or figures it generates.

Instructions are provided on the [same website](http://sushi.chem.kuleuven.be/svn/Localizer). For the Igor plugin, take a look at the "Getting Started" document, while for the Matlab plugin you can examine the example ".m" files.

Example data to play around with can be found [here for SOFI](http://sushi.chem.kuleuven.be/SOFI.tif.zip) and [here for localization microscopy](http://sushi.chem.kuleuven.be/PALM.pde.zip). The SOFI data shows a Hek cell transfected using Lyn-Dronpa, while the localization data shows a fixed HeLa cell that expressed Lyn-Dendra2.

At present I do not have a Linux or Macintosh version of the Matlab plugin. That is simply because I only have a Windows version of Matlab. Given Linux and/or Mac Matlab installations, it should be a rather easy to create suitable plugins. If this interests you and you have these installations, let me know!

Downloading and modifying the source code
=========================================
The Localizer code is licensed under the GPL license, which means that you are free to modify and distribute the code, provided that you provide proper attribution, and that any such code is governed by a compatible license (which means that your code must be similarly free). If this license is restrictive for your application, contact me and we may be able to work it out.

You can download the code using the links above. I use three different systems to compile the code:

*   [Xcode](https://developer.apple.com/xcode/) on Macintosh.
*   [Visual C++ Express](https://www.microsoft.com/visualstudio/en-us/products/2010-editions/visual-cpp-express) on Windows. For 64-bit development you will also want to install the [Windows SDK 7.1](https://www.microsoft.com/en-us/download/details.aspx?id=8279).
*   A simple shell script on Linux.

All of these tools are freely available, and can be downloaded from the provided links. The shell script is included with the source code. To compile the Igor Pro and Matlab plugins you will also need a copy of the [Igor XOP Toolkit](http://www.wavemetrics.com/products/xoptoolkit/xoptoolkit.htm) and/or a working Matlab installation. This is important – without the XOP Toolkit or a Matlab installation you will not be able to compile the plugins!

If you modify or improve this code in any way, I would appreciate if you could communicate these changes back to me so that it can be included for the benefit of other users (with proper attribution). For this it might help if you are familiar with [git](http://git-scm.com/).

License
=======
You're free to distribute this code or the compiled products in any way, provided that you provide proper attribution in your code and/or the products you derive from them. You're also free to include whatever parts of this code in your own work, provided that any distribution of this work is governed by a license compatible with the GPLv3, and subject to the same attribution clause. "Proper attribution" is a clearly visible notice containing e.g. "this product is based on Localizer, created by Peter Dedecker at the University of Leuven, available at https://bitbucket.org/pdedecker/localizer".

Compiling the plugin for Igor Pro
=================================
1.  Make sure that you have a copy of [Xcode](https://developer.apple.com/xcode/) (Mac) or [Visual C++ Express](https://www.microsoft.com/visualstudio/en-us/products/2010-editions/visual-cpp-express) and the [Windows SDK 7.1](https://www.microsoft.com/en-us/download/details.aspx?id=8279) (Windows).
1.  Get a copy of the [Igor XOP Toolkit](http://www.wavemetrics.com/products/xoptoolkit/xoptoolkit.htm) if you don't already have it, or make sure that you have the latest version otherwise.
1.  Download the Localizer source code and place the folder in the XOP Toolkit 6/IgorXOPs/ folder (in the same location as the examples provided by WaveMetrics).
1.  On Windows, open Visual Studio/Localizer.sln. On Mac, open Xcode/Localizer.
1.  Build the 'Localizer' projects in Visual Studio or XCode. Everything should compile just fine.
1.  Make whatever changes you see fit!

Compiling the plugin for Matlab
===============================
1.  Make sure that you have a copy of [Visual C++ Express](https://www.microsoft.com/visualstudio/en-us/products/2010-editions/visual-cpp-express), the [Windows SDK 7.1](https://www.microsoft.com/en-us/download/details.aspx?id=8279) (Windows), and a working copy of Matlab.
1.  Download the Localizer source code and place it in some convenient location.
1.  Open Visual Studio/Localizer.sln.
1.  Edit the project so that it knows where to find the required header files and libraries. 
    1.  Right click the 'LocalizerMatlab' project, and choose properties.
    1.  Make sure that you choose the correct platform in the top-right combobox: "Win32" if you're using a 32-bit version of Matlab, "x64" otherwise.
    1.  In the left list, click "C/C++" and then "General".
    1.  Edit the "Additional Include Directories" setting to point to your Matlab include directory. For example, if it currently reads "C:\Program Files\MATLAB\R2010b\extern\include", you may have to change it to point to a "R2011a" folder. Do not remove the other paths!
    1. In the left list, click "Linker" and then "General". Modify the "Additional Library Directories" to point to the correct folder for your Matlab installation.
1.  Build the 'LocalizerMatlab' project.
1.  Make whatever changes you see fit!

Compiling the command-line tool
===============================
Instructions on how to do this are in the "LinuxCompile.sh" script included with the source code.
