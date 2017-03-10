/**
  vcf_compare.cpp

  Purpose: finds all events in the first file that are represented by a similar event 
  (similar location, same or similar SV-length) in the second file. Additional arguments
  indicate what difference in location is considered similar enough (or rather: which is
  the minimum distance at which events are considered dissimilar), and whether one should
  ignore SV lengths in the comparison ('same_len' or 'ignore_len') 

  Usage: ./compare first_vcf second_vcf wiggle_room_bp whether_compare_lengths merged_vcf
  Example: ./compare pacbio_deletions.vcf freebayes_deletions.vcf 10 same_len pacbio_del_found_by_freebayes.vcf

  Contact data: Eric-Wubbo Lameijer, Xi'an Jiaotong University, eric_wubbo@hotmail.com
**/

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "shared_functions.h"

enum EventType { INS, DEL, SNP, RPL };

bool isInsertion(const std::string& ref, const std::string& alt) {
  return ((ref.length() == 1) && (alt.length() > 1 ));
}

bool isDeletion(const std::string& ref, const std::string& alt) {
  return ((ref.length() > 1) && (alt.length() == 1 ));
}

bool isPureInsertion(const std::string& ref, const std::string& alt) {
  return (isInsertion(ref,alt) && ref[0] == alt[0]);
}

bool isPureDeletion(const std::string& ref, const std::string& alt) {
  return (isDeletion(ref,alt) && ref[0] == alt[0]);
}

class Coordinate {

friend bool operator<(const Coordinate& leftCoordinate, const Coordinate& rightCoordinate);

public:
  Coordinate();
  Coordinate(const std::string& chromosomeName, int position);
  Coordinate(int chromosomeIndex, int position);
  Coordinate(const Coordinate& otherCoordinate);

  bool withinDistance(const Coordinate& otherCoordinate, int distance) const;
  int getDistanceBetween(const Coordinate& otherCoordinate) const; // warning: needs to be on same chromosome!
  Coordinate getDecreasedCoordinate(int distance) const;
  std::string getChromosomeName() const;
  int getPosition() const;

private:
  int m_chromosomeIndex;
  int m_position;
};

// dummy constructor, to avoid hassles with pointers elsewhere
Coordinate::Coordinate() :
  m_chromosomeIndex(0), m_position(0) {
}

Coordinate::Coordinate(const std::string& chromosomeName, int position) : 
  m_chromosomeIndex(chromosomeNameToIndex(chromosomeName)), m_position(position) {
}

Coordinate::Coordinate(int chromosomeIndex, int position) : 
  m_chromosomeIndex(chromosomeIndex), m_position(position) {
}

int Coordinate::getDistanceBetween(const Coordinate& otherCoordinate) const {
  Require(m_chromosomeIndex == otherCoordinate.m_chromosomeIndex, 
    "Coordinate::getDistanceBetween() error: can't compare coordinates on different chromosomes!");
  return abs(m_position - otherCoordinate.m_position);
}

bool operator<(const Coordinate& leftCoordinate, const Coordinate& rightCoordinate) {
  if (leftCoordinate.m_chromosomeIndex != rightCoordinate.m_chromosomeIndex) {
    return (leftCoordinate.m_chromosomeIndex < rightCoordinate.m_chromosomeIndex);
  } else {
    return (leftCoordinate.m_position < rightCoordinate.m_position);
  }
}

Coordinate::Coordinate(const Coordinate& otherCoordinate) :
  m_chromosomeIndex(otherCoordinate.m_chromosomeIndex), 
  m_position(otherCoordinate.m_position) {
}

std::string Coordinate::getChromosomeName() const {
  std::string id = "";
  if (m_chromosomeIndex == 100) {
    id = "X";
  } else if (m_chromosomeIndex == 101) {
    id = "Y";
  } else if (m_chromosomeIndex == 102) {
    id = "M";
  } else {
    id = intToString(m_chromosomeIndex);
  }
  return "chr"+id;
}

int Coordinate::getPosition() const {
  return m_position;
}

/** returns a clone of this coordinate, but with decreased start position **/
Coordinate Coordinate::getDecreasedCoordinate(int distance) const {
  int newPosition = m_position - distance;
  if (newPosition < 1) {
    newPosition = 1;
  }
  Coordinate newCoordinate(m_chromosomeIndex, newPosition);
  return newCoordinate;
}

/** Is this coordinate within distance "distance" of the other coordinate?
    For example 14 and 16 are within distance 3 of each other, but not within distance 2. **/
bool Coordinate::withinDistance(const Coordinate& otherCoordinate, int distance) const {
  if (m_chromosomeIndex != otherCoordinate.m_chromosomeIndex) {
    return false;
  } else {
    return abs(m_position - otherCoordinate.m_position) < distance;
  }
}

class Event {
  friend bool operator<(const Event& leftEvent, const Event& rightEvent);
  friend std::ostream& operator<<(std::ostream& os, const Event& event);

public:
  Event(const std::string& line);
  Event(const Coordinate& coordinate, const std::string& ref, const std::string& alt);
  Coordinate getCoordinate() const;
  EventType getType() const;
  int getSize() const;

private:
  Coordinate m_coordinate;
  std::string m_ref;
  std::string m_alt;
};

Event::Event(const std::string& line) {
  std::stringstream ss;
  ss << line;
  std::string chromosomeName;
  ss >> chromosomeName;
  int position;
  ss >> position;
  Coordinate tempCoordinate(chromosomeName, position);
  m_coordinate = tempCoordinate;

  std::string dummy;
  ss >> dummy;

  ss >> m_ref;
  ss >> m_alt;
}

Event::Event(const Coordinate& coordinate, const std::string& ref, const std::string& alt) {
  m_coordinate = coordinate;
  m_ref = ref;
  m_alt = alt;
}

Coordinate Event::getCoordinate() const {
  return m_coordinate;
}

bool operator<(const Event& leftEvent, const Event& rightEvent) {
  return (leftEvent.getCoordinate() < rightEvent.getCoordinate());
}

EventType Event::getType() const {
  if (m_ref.length() == 1 ) {
    if (m_alt.length() == 1) {
      return SNP;
    } else {
      if (isPureInsertion(m_ref, m_alt)) {
        return INS;
      }
    }
  } else {
    // reference length > 0
    if (isPureDeletion(m_ref, m_alt)) {
      return DEL;
    }
  }
  return RPL;
}

int Event::getSize() const {
  return abs(m_ref.length() - m_alt.length());
} 

std::ostream& operator<<(std::ostream& os, const Event& event) {
  os << event.getCoordinate().getChromosomeName() << "\t" << 
        event.getCoordinate().getPosition() << "\t" <<
        event.m_ref << "\t" <<
        event.m_alt << std::endl;
  return os;
}

bool sufficientlySimilar(const Event& currentEvent, const Event& soughtEvent, 
    int differenceDefiningDistance, bool requireSameSize) {
  if (currentEvent.getType() != soughtEvent.getType()) {
    return false;
  }
  if (requireSameSize) {
    if (currentEvent.getSize() != soughtEvent.getSize()) {
      return false;
    }
  }
  return (currentEvent.getCoordinate().getDistanceBetween(soughtEvent.getCoordinate()) < differenceDefiningDistance);
}

void transformFile(const std::string& nameOfComparedFile, const std::string& nameOfComparisonFile, 
    int wiggleRoom, bool requireIdenticalLengths, const std::string& nameOfOutputFile) {
  std::ifstream comparedFile(nameOfComparedFile.c_str());
  std::ifstream comparisonFile(nameOfComparisonFile.c_str());
  std::ofstream outputFile(nameOfOutputFile.c_str());

  std::string oldChrom = "";
  std::string oldPos = "";

  std::vector<Event> events;

  while (!comparisonFile.eof()) {
    std::string line;
    std::stringstream buffer_ss;
    getline(comparisonFile, line );
    if (line.length() == 0) {
      break;
    }

    // skip lines beginning with '#'
    const char START_OF_COMMENT_CHAR = '#';
    if (StringStartsWith(line,"#" )) {
      continue;
    } else {
      events.push_back(line);
    }
  }

  while (!comparedFile.eof()) {
    std::string line;
    std::stringstream buffer_ss;
    getline(comparedFile, line );
    if (line.length() == 0) {
      break;
    }

    // skip lines beginning with '#'
    const char START_OF_COMMENT_CHAR = '#';
    if (StringStartsWith(line,"#" )) {
      outputFile << line << "\n";
      continue;
    } 
    
    Event currentEvent(line);
    Coordinate lowerSearchBound = currentEvent.getCoordinate().getDecreasedCoordinate(wiggleRoom);
    Event lowerSearchDummyEvent(lowerSearchBound, "", "");

    // Now check for nearby events in the comparison-file
    std::vector<Event>::iterator firstEventToSearch;
    firstEventToSearch = std::lower_bound(events.begin(), events.end(), lowerSearchDummyEvent);
    std::vector<Event>::iterator eventToSearch = firstEventToSearch;
    std::vector<Event>::iterator bestEventIt = firstEventToSearch;
    int minDistance = wiggleRoom;
    while (eventToSearch->getCoordinate().withinDistance(currentEvent.getCoordinate(), wiggleRoom)) {
      if (sufficientlySimilar(*eventToSearch, currentEvent, wiggleRoom, requireIdenticalLengths)) {
         bestEventIt = eventToSearch;
         minDistance = eventToSearch->getCoordinate().getDistanceBetween(currentEvent.getCoordinate());
      } // if this event is a better match
      ++eventToSearch;
    } // while the event is still in the valid space
    if (minDistance < wiggleRoom) {
      outputFile << line << "\n";
    }
  }
  
  comparedFile.close();
  comparisonFile.close();
  outputFile.close();
}

bool argumentsCorrect(int argc, char** argv) {
  if (argc != 6) {
    return false;
  }
  if (atoi(argv[3]) == 0) {
    // wiggle room cannot be 0 or a string
    return false;
  }
  std::string lengthConsideration = argv[4];
  if (lengthConsideration != "same_len" && lengthConsideration != "ignore_len") {
    return false;
  }
  return true;
}


int main(int argc, char** argv) {
  if (!argumentsCorrect(argc,argv)) {
    std::cout <<
      "vcf_compare\n"
      "\n"
      "Purpose: finds all events in the first file that are represented by a similar event "
      "(similar location, same or similar SV-length) in the second file. Additional arguments "
      "indicate what difference in location is considered similar enough (or rather: which is "
      "the minimum distance at which events are considered dissimilar), and whether one should "
      "ignore SV lengths in the comparison ('same_len' or 'ignore_len').\n" 
      "\n"
      "Usage: ./compare first_vcf second_vcf wiggle_room_bp whether_compare_lengths merged_vcf\n"
      "Example: ./compare pacbio_deletions.vcf freebayes_deletions.vcf 10 same_len pacbio_del_found_by_freebayes.vcf\n"
      "\n"
      "Contact data: Eric-Wubbo Lameijer, Xi'an Jiaotong University, eric_wubbo@hotmail.com\n\n";
    return -1;
  } else {
    std::string nameOfFirstInputFile = argv[1];
    std::string nameOfSecondInputFile = argv[2];
    int wiggleRoom = atoi(argv[3]);
    std::string lengthConsideration = argv[4];
    bool requireIdenticalLengths = (lengthConsideration == "same_len");
    std::string nameOfOutputFile = argv[5];
    transformFile(nameOfFirstInputFile, nameOfSecondInputFile, wiggleRoom, requireIdenticalLengths, nameOfOutputFile);
    return 0;
  }
}
