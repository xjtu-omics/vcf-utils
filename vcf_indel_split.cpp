/**
  vcf_indel_split.cpp

  Purpose: from an input VCF file, creates two new files: one containing pure deletions (like ACTC -> A),
  another one containing pure insertions (like A -> ATTC). Note that SNPs and 'impure' indels 
  (like replacements, 'ACT -> AG') are not put into any output file.

  usage: ./indel_split input_vcf deletion_output_vcf insertion_output_vcf
  example: ./indel_split gatk_hanchild.vcf gatk_hanchild_deletions.vcf gatk_hanchild_insertions.vcf

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

void transformFile(const std::string& nameOfInputFile, const std::string& nameOfDeletionOutputFile,
  const std::string& nameOfInsertionOutputFile) {

  std::ifstream inputFile(nameOfInputFile.c_str());
  std::ofstream deletionOutputFile(nameOfDeletionOutputFile.c_str());
  std::ofstream insertionOutputFile(nameOfInsertionOutputFile.c_str());

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
      insertionOutputFile << line << std::endl;
      deletionOutputFile << line << std::endl;      
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
    
    /*if (StringStartsWith(genotype,"0/0") || StringStartsWith(genotype,".")) {
      std::cout << genotype << "\n";
    } else {
      outputFile << line << "\n";
    }*/
    if (isInsertion(ref,alt)) {
      insertionOutputFile << line << "\n";
    } else if (isDeletion(ref,alt)) {
      deletionOutputFile << line << "\n";
    } else {
      std::cout << "neither (pure) insertion nor (pure) deletion: " << line << "\n";
    }

  }     
  inputFile.close();
  deletionOutputFile.close();
  insertionOutputFile.close();
}


int main(int argc, char** argv) {
  //std::cout << "Converting the input VCF to output VCF";
  if (argc == 1) { std::cout << "indel_split\n"
    "\n"
    "Purpose: from an input VCF file, creates two new files: one containing pure deletions (like ACTC -> A),"
    "another one containing pure insertions (like A -> ATTC). Note that SNPs and 'impure' indels" 
    "(like replacements, 'ACT -> AG') are not put into any output file.\n"
    "\n"
    "usage: ./indel_split input_vcf deletion_output_vcf insertion_output_vcf\n"
    "example: ./indel_split gatk_hanchild.vcf gatk_hanchild_deletions.vcf gatk_hanchild_insertions.vcf\n";
    "\n"
    "contact data: Eric-Wubbo Lameijer, Xi'an Jiaotong University, eric_wubbo@hotmail.com\n\n"
    return 0;
  }
  if (argc < 4) {
    std::cerr << "Invalid number of arguments. At least three arguments "
        "are needed, the name of the input file, the name of the "
        "deletion output file, and the name of the insertion output file.\n";
    return -1;
  } else {
    std::string nameOfInputFile = argv[1];
    std::string nameOfDeletionOutputFile = argv[2];
    std::string nameOfInsertionOutputFile = argv[3];
    transformFile(nameOfInputFile, nameOfDeletionOutputFile, nameOfInsertionOutputFile);

    return 0;
  }
}
