/** 
  vcf_alt_unraveler.cpp

  Purpose: unravels events that have more than one alternative allele, as indicated by 
  comma-separation of alt alleles (for example "chr1 14053 A AT,AG"). This happens quite
  frequently in files produced by GATK, even though the alt calls themselves are
  (start of 2017) not necessarily very reliable from a Mendelian correctness point of view.
  Practically, multi-alt-calls also complicate further downstream data processing and
  analysis. While in many cases one may prefer to just get rid of such "multialt" calls by 
  using remove_double_alts instead, in cases where one really wants to keep all of the alts,
  this may be a better option for downstream compatibility (though it sacrifices the 
  higher factual correctness of the GATK calls, which in this regard could be considered
  superior to other SV callers, be it by complicating matters for VCF parsers)
 
  Usage: ./unravel_alts input_vcf output_vcf
  Example: ./unravel_alts gatk_hanchild.vcf gatk_hanchild_unraveled.vcf

  Contact data: Eric-Wubbo Lameijer, Xi'an Jiaotong University, eric_wubbo@hotmail.com
**/

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <set>
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

/** Splits a string pre character 'ch' into a vector of strings. So 4/5/3 split
 * on '/' would yield ["4","5","3"]. **/
std::vector<std::string> Split(const std::string& str, char separator) {

  std::vector<std::string> output;
  size_t startSearchPos = 0;
  while (true) {
    size_t separatorPos = str.find_first_of(separator, startSearchPos);
    //std::cout << "Getting " << str.substr(startSearchPos, separatorPos - startSearchPos) << std::endl;
    output.push_back(str.substr(startSearchPos, separatorPos- startSearchPos));
    if (separatorPos == std::string::npos) {
      break; // we're done
    }
    else {
      startSearchPos = separatorPos + 1;
    }
  };
  return output;
}

std::set<int> getAllUsedAlts(const std::string& line) {
  std::stringstream ss;
  ss << line;
  
  // skip items before genotypes
  for (int i = 0; i < 9; i++ ) {
    std::string dummy;
    ss >> dummy;
//std::cout << "Dummy " << dummy << std::endl;
  }

  std::set<int> usedAlts;

  while (ss) {
    std::string genotype;
    ss >> genotype;
    if (genotype == "") {
      // end of stream
      break;
    }
//std::cout << "GT " << genotype << std::endl;
    int colonPos = genotype.find_first_of(':');
    if ( colonPos != std::string::npos ) {
      genotype = genotype.substr(0,colonPos);
    }
//std::cout << "GT2 " << genotype << std::endl;
    std::vector<std::string> genotypeAsStrings = Split(genotype,'/');
    for (int i = 0; i < genotypeAsStrings.size(); ++i ) {
//std::cout << "GT3 " << genotype << std::endl;
      if (genotypeAsStrings[i] != "." && genotypeAsStrings[i] != "0") {
         usedAlts.insert( atoi(genotypeAsStrings[i].c_str()));
      }
    }
  }
//std::cout << "GT4 " << std::endl;
  return usedAlts;
}

std::string altify(const std::string& line, int altId, const std::string& correctAlt) {
  std::stringstream ss;
  std::stringstream outSs;
  ss << line;
  for (int i = 0; i < 4; i++ ) {
    // just output chrom, pos, id and ref normally
    std::string item;
    ss >> item;
    outSs << item << "\t";
  }
  std::string officialAlt;
  ss >> officialAlt;
  outSs << correctAlt << "\t"; // !important (of course) to replace the alt
  for (int i = 0; i < 4; i++ ) {
    // just output qual filter, infor and format normally
    std::string item;
    ss >> item;
    outSs << item << "\t";
  }
  while (ss) {
    // loop over the genotypes
    std::string genotype;
    ss >> genotype;
    int position;
    for (position = 0; position < genotype.length(); ++position) {
      // assumption: less than 9 alt alleles...
      char ch = genotype[position];
      if (ch == ':') {
         outSs << ch;
         break;
      }
      if (ch == altId + '0') {
        ch = '1';
      } else if (ch >= '0' && ch <= '9') {
        ch = '0';
      }
      outSs << ch;
    }
    ++position;
    while (position < genotype.length()) {
      outSs << genotype[position];
      ++position;
    } 
    if (ss) {
      outSs << '\t';
    }
  }
  return outSs.str();
}

void transformFile(const std::string& nameOfInputFile, const std::string& nameOfOutputFile) {
  std::ifstream inputFile(nameOfInputFile.c_str());
  std::ofstream outputFile(nameOfOutputFile.c_str());

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
      std::vector<std::string> alts = Split(alt,',');
      std::set<int> usedAltIds = getAllUsedAlts(line);
//std::cout << "X";
      for (std::set<int>::iterator usedAltIdIt = usedAltIds.begin(); usedAltIdIt != usedAltIds.end(); ++usedAltIdIt) {
         int altId = *usedAltIdIt;
         std::string correctAlt = alts[altId - 1]; // in 0/1, 1 refers to the first alt, so A T,C would be T, which is the 0th element of alts
         outputFile << altify(line, altId, correctAlt) << std::endl;
         //std::cout << altify(line, altId, correctAlt) << std::endl;
      }
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
      "unravel_alts\n"
      "\n"
      "Purpose: unravels events that have more than one alternative allele, as indicated by "
      "comma-separation of alt alleles (for example \"chr1 14053 A AT,AG\"). This happens quite "
      "frequently in files produced by GATK, even though the alt calls themselves are "
      "(start of 2017) not necessarily very reliable from a Mendelian correctness point of view. "
      "Practically, multi-alt-calls also complicate further downstream data processing and "
      "analysis. While in many cases one may prefer to just get rid of such \"multialt\" calls by "
      "using remove_double_alts instead, in cases where one really wants to keep all of the alts, "
      "this may be a better option for downstream compatibility (though it sacrifices the "
      "higher factual correctness of the GATK calls, which in this regard could be considered "
      "superior to many other SV callers, be it by complicating matters for VCF parsers).\n"
      "\n"
      "Usage: ./unravel_alts input_vcf output_vcf\n"
      "Example: ./unravel_alts gatk_hanchild.vcf gatk_hanchild_unraveled.vcf\n"
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
