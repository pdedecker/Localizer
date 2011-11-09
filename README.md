Localizer readme
================

Welcome to the source code for Localizer, which implements various kinds of analyses related to superresolution microscopy. The routines were lovingly hand-crafted by [Peter Dedecker](peter.dedecker@chem.kuleuven.be), who currently resides at the [Department of Chemistry](http://chem.kuleuven.be) at the [University of Leuven](http://www.kuleuven.be). This repository is available at https://github.com/pdedecker/Localizer.

Localizer is written in C++. This means that it runs very fast, can run on virtually any hardware, and does not require you to purchase any software to run. It was also designed to be modular. This means two things:

1. Localizer is fairly self contained and exposes its functionality through a small 'surface area'. This makes it especially suited for embedding in higher-level languages, using for example the MEX functionality in Matlab or the XOP Toolkit for Igor Pro. Right now only two wrappers are available: one that allows users to call Localizer from the Igor Pro language, and one that allows Localizer to be called from the command line (e.g. the terminal on Mac or Linux). I would be very interested in getting MEX bindings working, but haven't pursued this because I'm not all that familiar with Matlab. Do let me know if you would like to assist in this!
2. Localizer uses a 'building block' analogy for its algorithms. That means that you can introduce new algorithms for e.g. emitter localization or segmentation without a needing to change the source code. It also means that you could conceivably add an entirely new analysis strategy while making use of Localizer's functionality for e.g. reading CCD images.

This file accompanies the source code repository for Localizer. To compile your own version of the software, you will either need a copy of [XCode](https://developer.apple.com/xcode/) (on Macintosh) or [Visual Studio for C++](http://msdn.microsoft.com/en-us/express/future/bb421473). At present I am using XCode 4.1 and Visual C++ 2008 Express, both of which are freely available. You may be able to use different versions, but I have not tested this.

At present I have created front-ends for [Igor Pro](http://www.wavemetrics.com) and for a stand-alone command-line version. However, the command-line version can only be used on Macintosh or Linux right now, I have not looked into enabling this on Windows.

This also means that you will need to have a copy of the [Igor XOP Toolkit](http://www.wavemetrics.com/products/xoptoolkit/xoptoolkit.htm), unless you are using a Macintosh or Linux machine and are interested in the commandline version. You need version 6 of the Toolkit or later (free upgrade if you have ever purchased a license for the toolkit).

Visual Studio and XCode need to know where to find the file related to the XOP Toolkit. The easiest way to do this is to the place the project folder (the on this file is in) in the XOP Toolkit 6/IgorXOPs/ folder (next to the example projects provided by WaveMetrics). Once that is done you should be able to open the XCode or Visual Studio project, hit 'build', and obtain a working copy.

License
=======
I have made the Localizer software available under the GPL license. That means that you are free to modify and distribute the code, provided that you provide proper attribution, and that any such code is governed by a GPL-compatible license (which means that your code must be similarly free). If this license is restrictive for your application, contact me and we may be able to work it out.

If you modify or improve this code in any way, I would appreciate if you could communicate these changes back to me so that it can be included for the benefit of other users (with proper attribution). For this it might help if you are familiar with [git](http://git-scm.com/).

I'm not much for source code. Do you have something that I can start using right away?
======================================================================================

I have written a graphical user interface for Igor Pro that provides a fully-featured and capable program to use for analyzing your PALM/STORM or SOFI data, including visualization, clustering analysis, drift correction, etc. Igor Pro is commercial software, but [a free 30-day trial version is available for direct download](http://www.wavemetrics.com/support/demos.htm). Localizer will function without restrictions during the trial period. After the trial period Localizer will continue to function normally, but you will be unable to save the data or figures it generates.

You can download the necessary files [here](https://sushi.chem.kuleuven.be/SVN/Localizer).