#include <sstream>

#include "event.h"

#include "shared_functions.h"


Event::Event(const std::string& line) {
  //std::cout << "Constructing event from line: " << line << std::endl;
  std::stringstream ss;
  ss << line;
  std::string chromosomeName;
  ss >> chromosomeName;
  chromosomeIndex_ = chromosomeNameToIndex(chromosomeName);
  ss >> position_;
  std::string dummy;
  ss >> dummy;
  ss >> refAllele_;
  ss >> altAllele_;
}

bool operator<(const Event& leftEvent, const Event& rightEvent) {
  if (leftEvent.chromosomeIndex_ != rightEvent.chromosomeIndex_) {
    return (leftEvent.chromosomeIndex_ < rightEvent.chromosomeIndex_);
  } else if (leftEvent.position_ != rightEvent.position_) {
    return (leftEvent.position_ < rightEvent.position_);
  } else if (leftEvent.refAllele_ != rightEvent.refAllele_) {
    return (leftEvent.refAllele_ < rightEvent.refAllele_);
  } else {
    return (leftEvent.altAllele_ < rightEvent.altAllele_);
  }
}
