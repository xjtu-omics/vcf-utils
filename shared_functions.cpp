
#include <cstdlib>
#include <iostream> // debugging
#include <sstream>

#include "shared_functions.h"

#include "event.h"


int chromosomeNameToIndex(const std::string& chromosomeName) {
  //std::cout << chromosomeName << std::endl;
  std::string chromIdPart = chromosomeName.substr(3); // eliminate 'chr'
  if (chromIdPart == "X") {
    return 100;
  } else if (chromIdPart == "Y") {
    return 101;
  } else if (chromIdPart == "M") {
    return 102;
  } else {
    return atoi(chromIdPart.c_str());
  }
}


/* Returns whether a string starts with a certain other string, so if
   'stringToBeAssessed' is 'albert' and 'putativeStart' is 'al', this function
   returns true. */
bool StringStartsWith(const std::string& stringToBeAssessed,
    const std::string& putativeStart) {
  int subStringLength = putativeStart.length();
  std::string startOfAssessableString =
      stringToBeAssessed.substr(0,subStringLength);
  return (startOfAssessableString.compare(putativeStart) == 0 );
}

bool comesBefore(const std::string& firstLine, const std::string& secondLine) {
  Event firstEvent(firstLine);
  Event secondEvent(secondLine);
  return (firstEvent < secondEvent);
}

std::string intToString(int i) {
  std::stringstream ss;
  ss << i;
  std::string output;
  ss >> output;
  return output;
}

/** Utility function that halts/crashes the program, helps to catch bugs early.
**/
void Require(bool requirementMet, std::string errorMessage) {
  if (!requirementMet) {
    std::cerr << errorMessage << std::endl;
    exit(-1);
  }
}
