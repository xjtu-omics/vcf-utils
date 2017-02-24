/**
  vcf_sort.cpp

  Purpose: sorts a VCF file into the order chr1, chr2...chr22, chrX, chrY, chrM, as not
  all VCF files have this format (some have a format like chr1, chr11, chr12...chr19, chr2...)

  usage: ./sort original_vcf sorted_vcf 
  example: ./sort pacbio_hanchild_orig.vcf pacbio_hanchild_sorted.vcf

  contact data: Eric-Wubbo Lameijer, Xi'an Jiaotong University, eric_wubbo@hotmail.com
**/

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

#include <string>
#include <vector>

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

void transformFile(const std::string& nameOfInputFile, const std::string& nameOfOutputFile) {
  std::ifstream inputFile(nameOfInputFile.c_str());
  std::ofstream outputFile(nameOfOutputFile.c_str());
  std::vector<std::string> events;

  while (!inputFile.eof()) {
    std::string line;
    std::stringstream buffer_ss;
    getline(inputFile, line );
    if (line.length() == 0) {
      break;
    }

    // lines beginning with '#' are comment lines which need to be copied in their original sequence
    if (StringStartsWith(line,"#" )) {
      outputFile << line << "\n";
      continue;
    } else {
      events.push_back(line);
    }
  } // while not eof

  sort(events.begin(), events.end(), comesBefore);
  int numberOfEvents = events.size();
  for (int i = 0; i < events.size(); ++i) {
    outputFile << events[i] << std::endl;
  }
  
  inputFile.close();
  outputFile.close();
}

int main(int argc, char** argv) {
  if (argc == 1) {
    std::cout <<
      "sort\n"
      "\n"
      "Purpose: sorts a VCF file into the order chr1, chr2...chr22, chrX, chrY, chrM, as not "
      "all VCF files have this format (some have a format like chr1, chr11, chr12...chr19, chr2...)\n"
      "\n"
      "usage: ./sort original_vcf sorted_vcf \n"
      "example: ./sort pacbio_hanchild_orig.vcf pacbio_hanchild_sorted.vcf\n"
      "\n"
      "contact data: Eric-Wubbo Lameijer, Xi'an Jiaotong University, eric_wubbo@hotmail.com\n\n";
  } else if (argc < 3) {
    std::cout << "Invalid number of arguments. At least two arguments "
        "are needed, the name of the input file and the name of the "
        "output file.\n";
    return -1;
  } else {
    std::string nameOfInputFile = argv[1];
    std::string nameOfOutputFile = argv[2];
    transformFile(nameOfInputFile, nameOfOutputFile);
    return 0;
  }
}
