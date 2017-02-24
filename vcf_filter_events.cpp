/**
  vcf_filter_events.cpp

  Purpose: based on an existing VCF file, creates a new VCF file only containing those events that are given
  in a list of events (of the format chromosome:position:ref:alt, like "chr1:10:A:AT")

  usage: ./filter_events input_vcf event_list.txt output_vcf
  example: ./filter_events pindel_freebayes_merged_hanchild_del.vcf pfdel_shared_events.txt pindel_freebayes_merged_hanchild_shared_del.vcf

  contact data: Eric-Wubbo Lameijer, Xi'an Jiaotong University, eric_wubbo@hotmail.com
**/

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

void loadEvents(const std::string& nameOfFilterFile, std::set<std::string>& events) {
  std::ifstream filterFile(nameOfFilterFile.c_str());
  std::string line;
  while (!filterFile.eof()) {
    getline(filterFile, line);
    events.insert(line);
  }
}


void transformFile(const std::string& nameOfInputFile, const std::string& nameOfFilterFile, const std::string& nameOfOutputFile) {
  std::ifstream inputFile(nameOfInputFile.c_str());
  std::ofstream outputFile(nameOfOutputFile.c_str());
  std::set<std::string> events;
  loadEvents(nameOfFilterFile, events);

  while (!inputFile.eof()) {
    std::string line;
    std::stringstream buffer_ss;
    getline(inputFile, line );
    if (line.length() == 0) {
      break;
    }

    // lines beginning with '#' should simply be copied to the output
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

    std::string positionCode = chrom + ":" + pos + ":" + ref + ":" + alt;

    if (events.find(positionCode) == events.end() ) {
      std::cout << "removed, as not found in the list of filtered events " << line << std::endl;
    } else {
       events.erase(positionCode); // don't report an event twice
       outputFile << line << "\n";
    }    
  }     
  inputFile.close();
  outputFile.close();
}


int main(int argc, char** argv) {
  if (argc == 1) {
    std::cout <<
      "filter_events\n"
      "\n"
      "Purpose: based on an existing VCF file, creates a new VCF file only containing those events that are given "
      "in a list of events (of the format \"chromosome:position:ref:alt\", like \"chr1:10:A:AT\").\n"
      "\n"
      "usage: ./filter_events input_vcf event_list.txt output_vcf\n"
      "example: ./filter_events pindel_freebayes_merged_hanchild_del.vcf pfdel_shared_events.txt pindel_freebayes_merged_hanchild_shared_del.vcf\n"
      "\n"
      "contact data: Eric-Wubbo Lameijer, Xi'an Jiaotong University, eric_wubbo@hotmail.com\n\n";
  } else if (argc < 4) {
    std::cout << "Invalid number of arguments. At least three arguments "
        "are needed, the name of the input file, the name of the file containing the events that are to be maintained, and the name of the "
        "output file.";
    return -1;
  } else {
    std::string nameOfInputFile = argv[1];
    std::string nameOfFilterFile = argv[2];
    std::string nameOfOutputFile = argv[3];
    transformFile(nameOfInputFile, nameOfFilterFile, nameOfOutputFile);
  }
  return 0;
}
