/**
  vcf_remove_events.cpp

  Purpose: removes certain events from an input VCF file. Note that the events list can either be of the format "chromosome:position"
  (like "chr1:1023492") or of the format "chromosome:position:reference:alt" (like "chr1:1023492:A:AT"). This affects the precision 
  of filtering

  usage: ./remove_events input_vcf events.txt output_vcf
  example: ./remove_events pindel_hanchild.vcf pindel_hanchild_multialts.txt pindel_hanchild_deduplicated.vcf

  contact data: Eric-Wubbo Lameijer, Xi'an Jiaotong University, eric_wubbo@hotmail.com
**/

#include <cctype> // isdigit
#include <fstream>
#include <iostream>
#include <set>
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

void loadEventsToBeRemoved(const std::string& nameOfFilterFile, std::set<std::string>& eventsToBeRemoved, bool& extensive) {
  std::ifstream filterFile(nameOfFilterFile.c_str());
  std::string line;
  extensive = true; // assume 'chr1:12893123:A:AT"
  while (!filterFile.eof()) {
    getline(filterFile, line);
    char lastChar = line[line.length()-1];
    if (isdigit(lastChar)) {
      extensive = false;
    }
    eventsToBeRemoved.insert(line);
  }
}

void transformFile(const std::string& nameOfInputFile, const std::string& nameOfFilterFile, const std::string& nameOfOutputFile) {
  std::ifstream inputFile(nameOfInputFile.c_str());
  std::ofstream outputFile(nameOfOutputFile.c_str());
  std::set<std::string> eventsToBeRemoved;
  bool extensiveDescriptor = false; // "chr1:1209231:A:AT" is extensive, "chr1:1209231" is regular
  loadEventsToBeRemoved(nameOfFilterFile, eventsToBeRemoved, extensiveDescriptor);

  while (!inputFile.eof()) {
    std::string line;
    std::stringstream buffer_ss;
    getline(inputFile, line );
    if (line.length() == 0) {
      break;
    }

    // just copy comment lines (lines beginning with '#') to the output file. Those don't need to be
    // filtered
    const char START_OF_COMMENT_CHAR = '#';
    if (line[0] == START_OF_COMMENT_CHAR ) {
      outputFile << line << "\n";
      continue;
    }

    buffer_ss << line;
    std::string chrom;
    buffer_ss >> chrom;
    std::string pos;
    buffer_ss >> pos;
    std::string dummy;
    buffer_ss >> dummy;
    std::string ref;
    buffer_ss >> ref;
    std::string alt;
    buffer_ss >> alt;

    std::string positionCode = chrom + ":" + pos;
    if (extensiveDescriptor) {
      positionCode += ":" + ref + ":" + alt;
    }

    if (eventsToBeRemoved.find(positionCode) != eventsToBeRemoved.end() ) {
      std::cout << "removed " << line << std::endl;
    } else {
      outputFile << line << "\n";
    }
  }     
  inputFile.close();
  outputFile.close();
}


int main(int argc, char** argv) {
  std::cout << "Converting the input VCF to output VCF";
  if (argc == 1) {
    std::cout << 
      "remove_events\n"
      "\n"
      "Purpose: removes certain events from an input VCF file. Note that the events list can either be of the format \"chromosome:position\" "
      "(like \"chr1:1023492\") or of the format \"chromosome:position:reference:alt\" (like \"chr1:1023492:A:AT\"). This affects the precision "
      "of filtering.\n"
      "\n"
      "usage: ./remove_events input_vcf events.txt output_vcf\n"
      "example: ./remove_events pindel_hanchild.vcf pindel_hanchild_multialts.txt pindel_hanchild_deduplicated.vcf\n"
      "\n"
      "contact data: Eric-Wubbo Lameijer, Xi'an Jiaotong University, eric_wubbo@hotmail.com\n\n";
      return 0;
  } else if (argc < 4) {
    std::cout << "Invalid number of arguments. At least two arguments "
        "are needed, the name of the input file and the name of the "
        "output file.";
    return -1;
  } else {
    std::string nameOfInputFile = argv[1];
    std::string nameOfFilterFile = argv[2];
    std::string nameOfOutputFile = argv[3];
    transformFile(nameOfInputFile, nameOfFilterFile, nameOfOutputFile);

    	
    return 0;
  }
}
