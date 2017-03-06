/**
  vcf_find_mlma_events.cpp

  Purpose: finds all "multi-line, multi-alt events" (same chromosome, position, but dispersed over different lines in the VCF) 
  and writes a list containing them to standard output. This is especially handy for vcf files created from Pindel output,
  as in those files multiple alt alleles tend to be listed on different lines (in contrast to the same line, as is the case
  for GATK).

  usage: ./find_mlma input_vcf 
  example: ./find_mlma pindel_hanchild.vcf > pindel_hanchild_duplicates.txt

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

void transformFile(const std::string& nameOfInputFile) {
  std::ifstream inputFile(nameOfInputFile.c_str());

  std::string oldChrom = "";
  std::string oldPos = "";

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
    std::string pos;

    buffer_ss >> chrom;
    buffer_ss >> pos;

    if (chrom == oldChrom && pos == oldPos) {
       std::cout << chrom << ":" << pos << std::endl;
    }

    oldChrom = chrom;
    oldPos = pos;
  }     
  inputFile.close();
}


int main(int argc, char** argv) {
  //std::cout << "Converting the input VCF to output VCF";
  if (argc < 2) { 
    std::cout <<
      "find_mlma\n"
      "\n"
      "Purpose: finds all 'multi-line, multi-alt events' (same chromosome, position, but dispersed over different lines in the VCF) "
      "and writes a list containing them to standard output. This is especially handy for vcf files created from Pindel output, "
      "as in those files multiple alt alleles tend to be listed on different lines (in contrast to the same line, as is the case "
      "for GATK).\n"
      "\n"
      "usage: ./find_mlma input_vcf\n"
      "example: ./find_mlma pindel_hanchild.vcf > pindel_hanchild_duplicates.txt\n"
      "\n"
      "contact data: Eric-Wubbo Lameijer, Xi'an Jiaotong University, eric_wubbo@hotmail.com\n\n";
  } else {
    std::string nameOfInputFile = argv[1];
    transformFile(nameOfInputFile);  	
  }
  return 0;
}
