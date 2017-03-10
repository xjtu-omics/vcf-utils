#ifndef SHARED_FUNCTIONS_H
#define SHARED_FUNCTIONS_H

#include <string>

// changes a chromosome name into an integer, to sort chromosomes, even those
// with "weird" names like "chrX" and "chrY", into the desired sequence
int chromosomeNameToIndex(const std::string& chromosomeName);
bool comesBefore(const std::string& firstLine, const std::string& secondLine);
std::string intToString(int i);
void Require(bool requirementMet, std::string errorMessage);
bool StringStartsWith(const std::string& stringToBeAssessed, const std::string& putativeStart);

#endif // SHARED_FUNCTIONS_H
