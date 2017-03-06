/**
  vcf_find_duplicates.cpp

  Purpose: finds all duplicate events (same chromosome, position, ref and alt) and writes a list containing them to standard output

  usage: ./find_duplicates input_vcf 
  example: ./find_duplicates gatk_hanchild.vcf > gatk_hanchild_duplicates.txt

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
  std::string oldRef = "";
  std::string oldAlt = "";

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
      //outputFile << line << "\n";
      continue;
    }

    buffer_ss << line;
    std::string dummy;
    std::string chrom;
    std::string pos;

    buffer_ss >> chrom;
    buffer_ss >> pos;
 
    for (int i = 0; i < 1; i++ ) {
       buffer_ss >> dummy;
    }
    std::string ref;
    buffer_ss >> ref;
    std::string alt;
    buffer_ss >> alt;

    if (chrom == oldChrom && pos == oldPos && ref == oldRef && alt == oldAlt) {
       std::cout << chrom << ":" << pos << ":" << ref << ":" << alt << std::endl;
    }

    if (alt.find_first_of(',') != std::string::npos) {
      std::cout << "Ref: " << ref << " alt " << alt << "\n";
      continue;
    }

    oldChrom = chrom;
    oldPos = pos;
    oldRef = ref;
    oldAlt = alt;

    for (int i = 0; i < 4; i++ ) {
       buffer_ss >> dummy;
    }

    std::string genotype;
    buffer_ss >> genotype;
  }     
  inputFile.close();
}


int main(int argc, char** argv) {
  //std::cout << "Converting the input VCF to output VCF";
  if (argc < 2) { std::cout <<
      "find_dup"
      "\n"
      "Purpose: finds all duplicate events (same chromosome, position, ref and alt) and writes a list containing them to standard output.\n"
      "\n"
      "usage: ./find_dup input_vcf\n"
      "example: ./find_dup gatk_hanchild.vcf > gatk_hanchild_duplicates.txt\n"
      "\n"
      "contact data: Eric-Wubbo Lameijer, Xi'an Jiaotong University, eric_wubbo@hotmail.com\n\n";
  } else {
    std::string nameOfInputFile = argv[1];
    transformFile(nameOfInputFile);  	
  }
  return 0;
}
