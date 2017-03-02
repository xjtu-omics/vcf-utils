/**
  vcf_uniquify_loci.cpp

  Purpose: finds all events that share a locus (same chromosome and position) and removes all of those events except
  the first one, leaving only one event with those data in the file (so "chr1 10 AT A Caller=Pindel ... chr1 10 A AGGCGGC Caller=Delly ... 
  chr1 10 AT T Caller=GATK ...chr1 17 G GCC Caller=GATK" 
  becomes "chr1 10 AT T Caller=Pindel ... chr1 17 G GCC Caller=GATK"

  usage: ./uniquify_loci input_vcf output_vcf
  example: ./uniquify_loci pacbio_hanchild.vcf pacbio_hanchild_unique_loci.vcf

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

void transformFile(const std::string& nameOfInputFile, const std::string& nameOfOutputFile) {
  std::ifstream inputFile(nameOfInputFile.c_str());
  std::ofstream outputFile(nameOfOutputFile.c_str());

  std::string oldChrom = "";
  std::string oldPos = "";

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
      outputFile << line << std::endl;
      continue;
    }

    buffer_ss << line;
    
    std::string chrom;
    buffer_ss >> chrom;
    std::string pos;
    buffer_ss >> pos;

    if (chrom == oldChrom && pos == oldPos) {
       std::cout << chrom << ":" << pos << std::endl;
    } else {
       outputFile << line << std::endl;
    }

    oldChrom = chrom;
    oldPos = pos;
  }     
  inputFile.close();
  outputFile.close();
}


int main(int argc, char** argv) {
  
  if (argc < 3) { std::cout <<
    "uniquify_loci\n"
    "\n"
    "Purpose: finds all events that share a locus (same chromosome and position) and removes all of those events except "
    "the first one, leaving only one event with those data in the file (so \"chr1 10 AT A Caller=Pindel ... chr1 10 A AGGCGGC Caller=Delly ... "
    "chr1 10 AT T Caller=GATK ...chr1 17 G GCC Caller=GATK\" " 
    "becomes \"chr1 10 AT T Caller=Pindel ... chr1 17 G GCC Caller=GATK\".\n"
    "\n"
    "usage: ./uniquify_loci input_vcf output_vcf\n"
    "example: ./uniquify_loci pacbio_hanchild.vcf pacbio_hanchild_unique_events.vcf\n"
    "\n"
    "contact data: Eric-Wubbo Lameijer, Xi'an Jiaotong University, eric_wubbo@hotmail.com\n\n";
  } else {
    std::cout << "Converting the input VCF to output VCF.\n";
    std::string nameOfInputFile = argv[1];
    std::string nameOfOutputFile = argv[2];
    transformFile(nameOfInputFile, nameOfOutputFile);  
    std::cout << "Conversion completed.\n";	
  }
  return 0;
}
