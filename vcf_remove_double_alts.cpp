/** 
  vcf_remove_double_alts.cpp

  Purpose: removes events that have more than one alternative allele, as indicated by 
  comma-separation of alt alleles (for example "chr1 14053 A AT,AG"). This happens quite
  frequently in files produced by GATK, even though the alt calls themselves are
  (start of 2017) not necessarily very reliable from a Mendelian correctness point of view.
  Practically, multi-alt-calls also complicate further downstream data processing and
  analysis, so this tool can be used to remove them.
 
  Usage: ./remove_double_alts input_vcf output_vcf
  Example: ./remove_double_alts gatk_hanchild.vcf gatk_hanchild_wo_doublealts.vcf

  Contact data: Eric-Wubbo Lameijer, Xi'an Jiaotong University, eric_wubbo@hotmail.com
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
    
    for (int i = 0; i < 4; i++ ) {
       buffer_ss >> dummy;
    }

    std::string alt;
    buffer_ss >> alt;

    if (alt.find_first_of(',') != std::string::npos) {
      std::cout << "Alt " << alt << "\n";
    } else {
      outputFile << line << "\n";
    }
  }     
  inputFile.close();
  outputFile.close();
}


int main(int argc, char** argv) {
  if (argc == 1 ) {
    std::cout <<
      "remove_double_alts\n"
      "\n"
      "Purpose: removes events that have more than one alternative allele, as indicated by "
      "comma-separation of alt alleles (for example \"chr1 14053 A AT,AG\"). This happens quite "
      "frequently in files produced by GATK, even though the alt calls themselves are "
      "(start of 2017) not necessarily very reliable from a Mendelian correctness point of view. "
      "Practically, multi-alt-calls also complicate further downstream data processing and "
      "analysis, so this tool can be used to remove them.\n"
      "\n"
      "Usage: ./remove_double_alts input_vcf output_vcf\n"
      "Example: ./remove_double_alts gatk_hanchild.vcf gatk_hanchild_wo_doublealts.vcf\n"
      "\n"
      "Contact data: Eric-Wubbo Lameijer, Xi'an Jiaotong University, eric_wubbo@hotmail.com\n\n";
    return -1;
  }
  if (argc < 3) {
    std::cout << "Invalid number of arguments. At least two arguments "
        "are needed, the name of the input file and the name of the "
        "output file.";
    return -1;
  } else {
    std::string nameOfInputFile = argv[1];
    std::string nameOfOutputFile = argv[2];
    transformFile(nameOfInputFile, nameOfOutputFile);   	
    return 0;
  }
}
