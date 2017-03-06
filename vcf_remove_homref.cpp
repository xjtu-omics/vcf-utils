/**
  vcf_remove_homref.cpp

  Purpose: from an input VCF file that has only one sample, removes all lines containing homref-
  only events (like 0/0, or ./., or .). While usually a VCF file would not have this structure,
  homref lines can come into being after creating a single sample VCF file from a multi-sample
  vcf file, by for example using vcftools.

  usage: ./remove_homref input_vcf output_vcf
  example: ./remove_homref gatk_hanchild.vcf gatk_hanchild_without_homref.vcf

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
    std::string chrom;
    std::string pos;

    buffer_ss >> chrom;
    buffer_ss >> pos;
 
    if (chrom != oldChrom) {
      std::cout << "Chromosome: " << chrom << std::endl;
    }
    
    if (chrom == oldChrom && pos == oldPos) {
       std::cout << chrom << ":" << pos << std::endl;
    }
    oldChrom = chrom;
    oldPos = pos;

    for (int i = 0; i < 1; i++ ) {
       buffer_ss >> dummy;
    }
    std::string ref;
    buffer_ss >> ref;
    std::string alt;
    buffer_ss >> alt;

    if (alt.find_first_of(',') != std::string::npos) {
      std::cout << "Ref: " << ref << " alt " << alt << "\n";
      continue;
    }
    for (int i = 0; i < 4; i++ ) {
       buffer_ss >> dummy;
    }

    std::string genotype;
    buffer_ss >> genotype;
    
    if (StringStartsWith(genotype,"0/0") || StringStartsWith(genotype,".")) {
      std::cout << genotype << "\n";
    } else {
      outputFile << line << "\n";
    }
  }     
  inputFile.close();
  outputFile.close();
}


int main(int argc, char** argv) {
  //std::cout << "Converting the input VCF to output VCF";
  if (argc == 1) {
    std::cout << ""
      "remove_homref\n"
      "\n"
      "Purpose: from an input VCF file that has only one sample, removes all lines containing homref-"
      "only events (like 0/0, or ./., or .). While usually a VCF file would not have this structure, "
      "homref lines can come into being after creating a single sample VCF file from a multi-sample "
      "VCF file, for example by using vcftools.\n"
      "\n"
      "usage: ./remove_homref input_vcf output_vcf\n"
      "example: ./remove_homref gatk_hanchild.vcf gatk_hanchild_without_homref.vcf\n"
      "\n"
      "contact data: Eric-Wubbo Lameijer, Xi'an Jiaotong University, eric_wubbo@hotmail.com\n\n";
  } else if (argc < 3) {
    std::cerr << "Invalid number of arguments. At least two arguments "
        "are needed, the name of the input file and the name of the "
        "output file.";
    return -1;
  } else {
    std::string nameOfInputFile = argv[1];
    std::string nameOfOutputFile = argv[2];
    transformFile(nameOfInputFile, nameOfOutputFile);  
  }
  return 0;
}
