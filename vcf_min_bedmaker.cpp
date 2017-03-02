/**
  vcf_min_bedmaker.cpp

  Purpose: makes a BED-file from a VCF, taking as the coordinates of each BED line
  the start coordinate of the event and the startposition + eventsize.
  Note that these assumptions can be problematic: it does not take into account that
  in repetitive regions an indel can be on many loci, and even the slightest error in
  reference or read can shift it greatly, or that an insertion should actually have
  a size of one, but that would not be helpful for establishing overlaps in noisy
  regions.

  usage: ./min_bedmaker input_vcf output_bed
  example: ./min_bedmaker pacbio_hanchild.vcf pacbio_hanchild.bed

  contact data: Eric-Wubbo Lameijer, Xi'an Jiaotong University, eric_wubbo@hotmail.com
**/

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

int eventLength(const std::string& ref, const std::string alt) {
  int length = ref.length() - alt.length();
  if (length < 0 ) {
    length = -length;
  }
  return length;
}

void transformFile(const std::string& nameOfInputFile, const std::string& nameOfOutputFile) {
  std::ifstream inputFile(nameOfInputFile.c_str());
  std::ofstream outputFile(nameOfOutputFile.c_str());

  std::string oldChrom = "";
  std::string oldPos = "";
  std::string oldRef = "";
  std::string oldAlt = "";

  while (!inputFile.eof()) {
    std::string line;
    std::stringstream buffer_ss;
    getline(inputFile, line );
    if (line.length() == 0) {
      return;
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
    std::string dummy;
    buffer_ss >> dummy;
    std::string ref;
    buffer_ss >> ref;
    std::string alt;
    buffer_ss >> alt;

    outputFile << chrom << "\t" << pos << "\t" << pos + eventLength(ref,alt) + 1 << std::endl;
  }     
  inputFile.close();
  outputFile.close();
}


int main(int argc, char** argv) {
  
  if (argc < 3) { std::cout <<
    "min_bedmaker\n"
    "\n"
    "Purpose: makes a BED-file from a VCF, taking as the coordinates of each BED line "
    "the start coordinate of the event and the startposition + eventsize.\n"
    "Note that these assumptions can be problematic: it does not take into account that "
    "in repetitive regions an indel can be on many loci, and even the slightest error in "
    "reference or read can shift it greatly, or that an insertion should actually have "
    "a size of one, but that would not be helpful for establishing overlaps in noisy "
    "regions.\n"
    "\n"
    "Usage: ./min_bedmaker input_vcf output_bed\n"
    "Example: ./min_bedmaker pacbio_hanchild.vcf pacbio_hanchild.bed\n"
    "\n"
    "contact data: Eric-Wubbo Lameijer, Xi'an Jiaotong University, eric_wubbo@hotmail.com\n\n";
  } else {
    std::cout << "Converting the input VCF to output BED.\n";
    std::string nameOfInputFile = argv[1];
    std::string nameOfOutputFile = argv[2];
    transformFile(nameOfInputFile, nameOfOutputFile);  
    std::cout << "Conversion completed.\n";	
  }
  return 0;
}
