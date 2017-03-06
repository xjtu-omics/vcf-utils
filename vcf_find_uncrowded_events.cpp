/**
  vcf_find_uncrowded_events.cpp

  Purpose: finds all events that have no events X basepairs before or behind them. So
  if X=1 will from events chr1 100.. chr1 102 ... chr1 103 ... chr1 105 eliminate 
  102 and 103, as their difference in position is less or equal to 1.

  usage: ./find_uncrowded input_vcf free_space
  example: ./find_uncrowded pindel_hanchild.vcf 100 > pindel_hanchild_uncrowded100.txt

  contact data: Eric-Wubbo Lameijer, Xi'an Jiaotong University, eric_wubbo@hotmail.com
**/

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

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

void transformFile(const std::string& nameOfInputFile, int windowSize) {
  std::ifstream inputFile(nameOfInputFile.c_str());

  std::string oldChrom = "";
  int oldPos = 0;
  bool predecessorPreGapOk = false;

  // Basically: if in the original file, an event has no event within X basepairs before or behind them,
  // select the event

  // Note that to know that, you must know both the previous and next line of the event
  // As reading the next line is a bit tricky, let's keep the last coordinates in memory (that's what we need to remember anyway)
  // which is automatically done with oldChrom and oldPos

  // as well as if they themselves were sufficiently distant from their predecessor

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
      continue;
    }

    buffer_ss << line;

    std::string chrom;
    buffer_ss >> chrom;

    int pos;
    buffer_ss >> pos;

    if (chrom != oldChrom || pos > oldPos + windowSize) {
       // well, the new gap is big enough
       if (predecessorPreGapOk) {
         std::cout << oldChrom << ":" << oldPos << std::endl;
       }
       predecessorPreGapOk = true; // in any case make it true
    } else {
       predecessorPreGapOk = false;
    }

    oldChrom = chrom;
    oldPos = pos;
  }    
  // save last event - if applicable
  if (predecessorPreGapOk) {
    std::cout << oldChrom << ":" << oldPos << std::endl;
  } 
  inputFile.close();
}


int main(int argc, char** argv) {
  //std::cout << "Converting the input VCF to output VCF";
  if (argc != 3) { 
    std::cout <<
      "find_uncrowded\n"
      "\n"
      "Purpose: finds all events that have no events X basepairs before or behind them. So "
      "if X=1 will from events chr1 100.. chr1 102 ... chr1 103 ... chr1 105 eliminate "
      "102 and 103, as their difference in position is less or equal to 1.\n"
      "\n"
      "usage: ./find_uncrowded input_vcf free_space\n"
      "example: ./find_uncrowded pindel_hanchild.vcf 100 > pindel_hanchild_uncrowded100.tx\n"
      "\n"
      "contact data: Eric-Wubbo Lameijer, Xi'an Jiaotong University, eric_wubbo@hotmail.com\n\n";
  } else {
    std::string nameOfInputFile = argv[1];
    int windowSize = atoi(argv[2]);
    transformFile(nameOfInputFile, windowSize);  	
  }
  return 0;
}
