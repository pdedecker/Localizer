#! /bin/bash
# Run this script to compile the commandline version of Localizer. You need development versions of the following:
# The GNU Scientific Library (GSL)
# Boost and libraries boost_thread, boost_date_time, boost_program_options and boost_system
# FFTW3
# Eigen3
# cblas
# atlas
# libtiff
# these should all be available in the package repository associated with your linux distribution.
# after compilation you will find a binary called LocalizationAnalysis, that you can then copy to e.g. usr/local/lib.
# for usage instructions, execute "LocalizationAnalysis --help"
echo "Compiling LocalizationAnalysis"
g++ -o LocalizationAnalysis -O2 -ffast-math -D NDEBUG -mtune=core2 -lgsl -lcblas -latlas -lboost_date_time -lboost_filesystem -lboost_thread -lfftw3 -ltiff -lboost_program_options ./source/PALM_analysis_CommandLine.cpp ./source/PALM_analysis.cpp ./source/PALM_analysis_Convolver.cpp ./source/PALM_analysis_FileIO.cpp ./source/PALM_analysis_Localization.cpp ./source/PALM_analysis_PALMImages.cpp ./source/PALM_analysis_ParticleFinding.cpp ./source/PALM_analysis_Processing.cpp ./source/PALM_analysis_segmentation.cpp ./source/PALM_analysis_storage.cpp ./source/PALM_analysis_MatrixRecycler.cpp
