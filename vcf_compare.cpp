/**
  vcf_compare.cpp

  Purpose: finds all events in the first file that are represented by a similar event 
  (similar location, same SV-length) in the second file.

  Usage: ./compare first_vcf second_vcf merged_vcf 
  Example: ./compare pacbio_deletions.vcf freebayes_deletions.vcf pacbio_del_found_by_freebayes.vcf

  Contact data: Eric-Wubbo Lameijer, Xi'an Jiaotong University, eric_wubbo@hotmail.com
**/

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

#include <string>
#include <vector>

enum EventType { INS, DEL, SNP, RPL };

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

bool isInsertion(const std::string& ref, const std::string& alt) {
  return ((ref.length() == 1) && (alt.length() > 1 ));
}

bool isDeletion(const std::string& ref, const std::string& alt) {
  return ((ref.length() > 1) && (alt.length() == 1 ));
}

bool isPure(const std::string& ref, const std::string& alt) {
  return (ref[0] == alt[0]);
}

bool isPureInsertion(const std::string& ref, const std::string& alt) {
  return (isInsertion(ref,alt) && isPure(ref,alt));
}

bool isPureDeletion(const std::string& ref, const std::string& alt) {
  return (isDeletion(ref,alt) && isPure(ref,alt));
}

int extractChromNo(const std::string& line) {
  std::stringstream ss;
  ss << line;
  std::string chromAsString;
  ss >> chromAsString;
  std::string chromIdPart = chromAsString.substr(3); // eliminate 'chr'
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

int chromNameToNumber(const std::string& chromNameAsString) {
  std::string chromIdPart = chromNameAsString.substr(3); // eliminate 'chr'
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

int extractPos(const std::string& line) {
  std::stringstream ss;
  ss << line;
  std::string chromAsString;
  ss >> chromAsString;
  int position;
  ss >> position;
  return position;
}

bool comesBefore(const std::string& firstLine, const std::string& secondLine) {
  int chromNoFirstEvent = extractChromNo(firstLine);
  int chromNoSecondEvent = extractChromNo(secondLine);
  if (chromNoFirstEvent != chromNoSecondEvent) {
    return (chromNoFirstEvent < chromNoSecondEvent);
  }
  int posFirstEvent = extractPos(firstLine);
  int posSecondEvent = extractPos(secondLine);
  return (posFirstEvent < posSecondEvent);
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
  m_chromosomeIndex(chromNameToNumber(chromosomeName)), m_position(position) {
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

void transformFile(const std::string& nameOfComparedFile, const std::string& nameOfComparisonFile, int wiggleRoom,
  const std::string& nameOfOutputFile) {
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
      if (eventToSearch->getType() == currentEvent.getType() && 
          //eventToSearch->getSize() == currentEvent.getSize() &&
          eventToSearch->getCoordinate().getDistanceBetween(currentEvent.getCoordinate()) < minDistance
          ) {
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


int main(int argc, char** argv) {
  if (argc != 5) {
    std::cout <<
      "vcf_compare\n"
      "\n"
      "Purpose: finds all events in the first file that are represented by a similar event "
      "(similar location, same SV-length) in the second file.\n"
      "\n"
      "Usage: ./compare first_vcf second_vcf wiggle_room_bp merged_vcf\n"
      "Example: ./compare pacbio_deletions.vcf freebayes_deletions.vcf 10 pacbio_del_found_by_freebayes.vcf\n"
      "\n"
      "Contact data: Eric-Wubbo Lameijer, Xi'an Jiaotong University, eric_wubbo@hotmail.com\n\n";
    return -1;
  } else {
    std::string nameOfFirstInputFile = argv[1];
    std::string nameOfSecondInputFile = argv[2];
    int wiggleRoom = atoi(argv[3]);
    std::string nameOfOutputFile = argv[4];
    transformFile(nameOfFirstInputFile, nameOfSecondInputFile, wiggleRoom, nameOfOutputFile);
    return 0;
  }
}
