/**
  vcf_eventizer.cpp

  Purpose: transforms a VCF file into a list of events, handy for later use in filtering.
  Events can be "wide" (chrom:pos:ref:alt, like "chr1:10:A:AT") or "narrow" (chrom:pos,
  like "chr1:10"). ("wide" or "narrow" need to be given as the second command line parameter).

  Usage: ./eventizer input_vcf wideness_flag
  Example: ./eventizer found_pacbio_events.vcf wide found_pacbio_events.txt

  Contact data: Eric-Wubbo Lameijer, Xi'an Jiaotong University, eric_wubbo@hotmail.com
**/

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

void transformFile(const std::string& nameOfInputFile, const std::string& format, 
    const std::string& nameOfOutputFile) {

  std::ifstream inputFile(nameOfInputFile.c_str());
  std::ofstream outputFile(nameOfOutputFile.c_str());
  bool isWide = (format == "wide");

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
    
    std::string dummy;
    buffer_ss >> dummy;

    std::string ref;
    buffer_ss >> ref;
    std::string alt;
    buffer_ss >> alt;
    
    outputFile << chrom << ":" << pos;
    if (isWide) {
      outputFile << ":" << ref << ":" << alt;
    }
    outputFile << std::endl;
  }     
  inputFile.close();
  outputFile.close();
}

bool parametersOkay(int argc, char** argv) {
  if (argc != 4) {
    return false;
  }
  std::string secondParameter = argv[2];
  return ((secondParameter == "wide") || (secondParameter == "narrow"));
}


int main(int argc, char** argv) {
  //std::cout << "Converting the input VCF to output VCF";
  if (!parametersOkay(argc, argv)) { 
    std::cout <<
      "eventizer\n"
      "\n"
      "Purpose: transforms a VCF file into a list of events, handy for later use in filtering. "
      "Events can be \"wide\" (chrom:pos:ref:alt, like \"chr1:10:A:AT\") or \"narrow\" (chrom:pos, "
      "like \"chr1:10\"). (\"wide\" or \"narrow\" need to be given as the second command line parameter).\n"
      "\n"
      "Usage: ./eventizer input_vcf wideness_flag output_txt\n"
      "Example: ./eventizer found_pacbio_events.vcf wide found_pacbio_events.txt\n"
      "\n"
      "Contact data: Eric-Wubbo Lameijer, Xi'an Jiaotong University, eric_wubbo@hotmail.com\n\n";
  } else {
    std::string nameOfInputFile = argv[1];
    std::string wideness = argv[2];
    std::string nameOfOutputFile = argv[3];
    transformFile(nameOfInputFile, wideness, nameOfOutputFile);  	
  }
  return 0;
}
