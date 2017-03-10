/** 
  vcf_filter_eventtypes.cpp

  Purpose: does some basic filtering, like only keeping events with a certain
  min length, max length, or event type (INS/DEL/SNP/ALL)
 
  Usage: ./filter_eventtypes input_vcf min_size max_size event_type output_vcf
  Example: ./filter_eventtypes pacbio_hanchild.vcf 1 1000 ALL pacbio_hanchild_maxsize1000.vcf

  Contact data: Eric-Wubbo Lameijer, Xi'an Jiaotong University, eric_wubbo@hotmail.com
**/

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

enum EventType {SNP, INS, DEL, ALL};

EventType stringToEventType(const std::string& eventTypeAsString) {
  if (eventTypeAsString == "SNP") {
    return SNP;
  } else if (eventTypeAsString == "INS") {
    return INS;
  } else if (eventTypeAsString == "DEL") {
    return DEL;
  } else if (eventTypeAsString == "ALL") {
    return ALL;
  } else {
    std::cerr << eventTypeAsString << " is not recognized as an eventtype." << std::endl;
    exit(-1);
  }
}

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
  return (isDeletion(ref,alt) && alt[0] == ref[0]);
}

int changeInSize(const std::string& ref, const std::string& alt) {
  return abs(ref.length() - alt.length());
}

bool isHomopolymer(const std::string& ref, const std::string& alt) {
  if (isPureInsertion(ref,alt) || isPureDeletion(ref,alt)) {
    std::string alleleToBeInvestigated = (ref.length() > alt.length()) ? ref : alt;
    std::string sequenceToBeInvestigated = alleleToBeInvestigated.substr(1);
    char homopolyCandidateChar = sequenceToBeInvestigated[0];
    for (int position = 1; position < sequenceToBeInvestigated.length(); ++position) {
      if (sequenceToBeInvestigated[position] != homopolyCandidateChar) {
        return false;
      }
    }
    return true;
  } else {
    // not a neat insertion or deletion? Don't judge it a homopolymer
    return false;
  }
}

bool isEventType(const std::string& ref, const std::string& alt, EventType eventType) {
  if (eventType == ALL) {
    return true;
  } else if (eventType == SNP) {
    return (ref.length() == 1 && alt.length() == 1);
  } else if (eventType == DEL) {
    return isPureDeletion(ref, alt);
  } else if (eventType == INS) {
    return isPureInsertion(ref, alt);
  } else {
    std::cerr << "unknown eventtype!" << std::endl;
    return false;
  }
}

void transformFile(const std::string& nameOfInputFile, int minSize, int maxSize, EventType eventType, const std::string& nameOfOutputFile) {
  std::ifstream inputFile(nameOfInputFile.c_str());
  std::ofstream outputFile(nameOfOutputFile.c_str());

  while (!inputFile.eof()) {
    std::string line;
    std::stringstream buffer_ss;
    getline(inputFile, line );
    if (line.length() == 0) {
      break;
    }

    // skip lines beginning with '#'
    const char START_OF_COMMENT_CHAR = '#';
    if (line[0] == START_OF_COMMENT_CHAR ) {
      outputFile << line << "\n";
      continue;
    }

    buffer_ss << line;
    std::string dummy;
    
    for (int i = 0; i < 3; i++ ) {
       buffer_ss >> dummy;
    }
    std::string ref;
    buffer_ss >> ref;
 
    std::string alt;
    buffer_ss >> alt;

    int sizeChange = changeInSize(ref,alt);

    /*if (sizeChange == 1 || 
        (sizeChange < 4) && isHomopolymer(ref,alt)) {*/
    if (sizeChange < minSize || sizeChange > maxSize || !isEventType(ref, alt, eventType)) {
      std::cout << "Filtered out: " << ref << ", " << alt << std::endl;
    } else {
      outputFile << line << std::endl;
    }
  }     
  inputFile.close();
  outputFile.close();
}

bool commandLineArgumentsValid(int argc, char** argv) {
  if (argc != 6) {
    return false;
  }
  std::string eventType = argv[4];
  if (eventType != "SNP" && eventType != "INS" && eventType != "DEL" && eventType != "ALL") {
    return false;
  }
  return true;
}


int main(int argc, char** argv) {
  if (!commandLineArgumentsValid(argc, argv)) {
    std::cout <<
      "filter_eventtypes\n"
      "\n"
      "Purpose: does some basic filtering, like only keeping events with a certain "
      "min length, max length, or event type (INS/DEL/SNP/ALL)\n"
      "\n"
      "Usage: ./filter_eventtypes input_vcf min_size max_size event_type output_vcf\n"
      "Example: ./filter_eventtypes pacbio_hanchild.vcf 1 1000 ALL pacbio_hanchild_maxsize1000.vcf\n"
      "\n"
      "Contact data: Eric-Wubbo Lameijer, Xi'an Jiaotong University, eric_wubbo@hotmail.com\n\n";
    return -1;
  } else {
    std::string nameOfInputFile = argv[1];
    int minLength = atoi(argv[2]);
    int maxLength = atoi(argv[3]);
    std::string eventTypeString = argv[4];
    EventType eventType = stringToEventType(eventTypeString);
    std::string nameOfOutputFile = argv[5];
    transformFile(nameOfInputFile, minLength, maxLength, eventType, nameOfOutputFile);   	
    return 0;
  }
}
