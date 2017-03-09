/**
  vcf_fuse.cpp

  Purpose: takes two VCF files which have the same event(s) in the same sequence, and merges/fuses them,
  producing a VCF file that has the combined events sorted into the correct places.
  Note that if an event occurs in both VCF files, it will also occur twice (two lines, right above each other)
  in the resulting VCF.

  usage: ./fuse first_vcf second_vcf merged_vcf 
  example: ./fuse pindel_del.vcf freebayes_del.vcf pindel_freebayes_merged_del.vcf

  contact data: Eric-Wubbo Lameijer, Xi'an Jiaotong University, eric_wubbo@hotmail.com
**/

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "shared_functions.h"


void transformFile(const std::string& nameOfFirstInputFile, const std::string& nameOfSecondInputFile, const std::string& nameOfOutputFile) {
  std::ifstream firstInputFile(nameOfFirstInputFile.c_str());
  std::ifstream secondInputFile(nameOfSecondInputFile.c_str());
  std::ofstream outputFile(nameOfOutputFile.c_str());

  std::string oldChrom = "";
  std::string oldPos = "";

  std::vector<std::string> events;

  while (!firstInputFile.eof()) {
    std::string line;
    std::stringstream buffer_ss;
    getline(firstInputFile, line );
    if (line.length() == 0) {
      break;
    }

    // skip lines beginning with '#'
    const char START_OF_COMMENT_CHAR = '#';
    if (StringStartsWith(line,"##" )) {
      outputFile << line << "\n";
      continue;
    } else if (StringStartsWith(line,"#")) { // "#CHROM
      // use the #CHROM of the second file
      continue;
    } else {
      events.push_back(line);    
    }
  }

  while (!secondInputFile.eof()) {
    std::string line;
    std::stringstream buffer_ss;
    getline(secondInputFile, line );
    if (line.length() == 0) {
      break;
    }

    // skip lines beginning with '#'
    const char START_OF_COMMENT_CHAR = '#';
    if (StringStartsWith(line,"#" )) {
      outputFile << line << "\n";
      continue;
    } else {
      events.push_back(line);
    }
  }
  sort(events.begin(), events.end(), comesBefore);
  int numberOfEvents = events.size();
  for (int i = 0; i < events.size(); ++i) {
    outputFile << events[i] << std::endl;
  }
  
  firstInputFile.close();
  secondInputFile.close();
  outputFile.close();
}


int main(int argc, char** argv) {
  if (argc == 1) {
    std::cout <<
      "fuse\n"
      "\n"
      "Purpose: takes two VCF files which have the same event(s) in the same sequence, and merges/fuses them, "
      "producing a VCF file that has the combined events sorted into the correct places.\n"
      "Note that if an event occurs in both VCF files, it will also occur twice (two lines, right above each other) "
      "in the resulting VCF.\n"
      "\n"
      "usage: ./fuse first_vcf second_vcf merged_vcf\n"
      "example: ./fuse pindel_del.vcf freebayes_del.vcf pindel_freebayes_merged_del.vcf\n"
      "\n"
      "contact data: Eric-Wubbo Lameijer, Xi'an Jiaotong University, eric_wubbo@hotmail.com\n\n";
  } else if (argc < 4) {
    std::cout << "Invalid number of arguments. At least three arguments "
        "are needed, the name of the first input file, the name of the second input file, and the name of the "
        "output file.";
    return -1;
  } else {
    std::string nameOfFirstInputFile = argv[1];
    std::string nameOfSecondInputFile = argv[2];
    std::string nameOfOutputFile = argv[3];
    transformFile(nameOfFirstInputFile, nameOfSecondInputFile, nameOfOutputFile);
    return 0;
  }
}
