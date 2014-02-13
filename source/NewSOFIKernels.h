#ifndef NEWSOFIKERNELS_H
#define NEWSOFIKERNELS_H

#include <vector>
#include <utility>

//typedef std::vector<std::pair<int, int> > PixelCombination;
//typedef std::vector<PixelCombination> Subset;
//typedef std::vector<Subset> Partition;
//typedef std::vector<Partition> GroupOfPartitions;

typedef std::vector<std::pair<int, int> > PixelCombination;
typedef std::vector<PixelCombination> Partition;
typedef std::vector<Partition> GroupOfPartitions;

class SOFIPixelCombination {
public:
    int outputDeltaX;
    int outputDeltaY;
    std::vector<PixelCombination> combinations;
};

class SOFIKernel {
public:
    int outputDeltaX;
    int outputDeltaY;
    std::vector<GroupOfPartitions> combinations;
};

std::vector<SOFIKernel> PixelCombinationsToKernels(const std::vector<SOFIPixelCombination>& sofiPixelCombination);

class ComparePixelCombinations {
public:
    bool operator() (const PixelCombination& comb1, const PixelCombination& comb2) const;
};

GroupOfPartitions AllPartitions(PixelCombination pixelCombination);
std::string PrintPartition(const Partition& partition);
std::string ParititionTest();

std::vector<SOFIPixelCombination> sofiPixelCombinations2();
std::vector<SOFIPixelCombination> sofiPixelCombinations3();
std::vector<SOFIPixelCombination> sofiPixelCombinations4();
std::vector<SOFIPixelCombination> sofiPixelCombinations5();
std::vector<SOFIPixelCombination> sofiPixelCombinations6();

#endif