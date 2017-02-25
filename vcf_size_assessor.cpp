/**
  vcf_size_assessor.cpp

  Purpose: returns a table summarizing the occurence of events of certain
  sizes in a VCF file. Its output is like
  
  Size    Count   Caller  SvType
  1       107774  Pindel  deletion
  2       36753   Pindel  deletion
  3       14359   Pindel  deletion
  ...     ...     ...     ... 

  Note that it should only be used on a file that contains only one type
  of event (insertion or deletion)

  Usage: ./size_ass input_vcf output_txt name_of_caller name_of_sv_type
  Example: ./size_ass pindel_hanchild_del.vcf pindel_hanchild_del_sizes.txt Pindel deletion

  Contact data: Eric-Wubbo Lameijer, Xi'an Jiaotong University, eric_wubbo@hotmail.com
**/

#include <fstream>
#include <iostream>
#include <map>
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

int getEventSize(const std::string& ref, const std::string& alt) {
  int refLength = ref.size();
  int altLength = alt.size();
  int eventSize = refLength - altLength;
  if (eventSize < 0) { // deletion
    eventSize = -eventSize;
  }
  return eventSize;
}

void transformFile(const std::string& nameOfInputFile, const std::string& nameOfOutputFile,
  const std::string& nameOfCaller, const std::string& nameOfSvType) {
  std::ifstream inputFile(nameOfInputFile.c_str());
  std::ofstream outputFile(nameOfOutputFile.c_str());

  std::string oldChrom = "";
  std::string oldPos = "";

  std::map<int,int> sizeCounts;

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
    
    int size = getEventSize(ref,alt);
    std::cout << "Size: " << size << ", mapsize: " << sizeCounts.size() << std::endl;
    if (sizeCounts.count(size) == 0 ) {
      sizeCounts[size] = 1;
    } else {
      ++sizeCounts[size];
    }
  }     
  inputFile.close();
  outputFile << "Size\tCount\tCaller\tSvType" << std::endl;
  for (std::map<int,int>::iterator it = sizeCounts.begin(); it != sizeCounts.end(); ++it) {
    outputFile << it->first << "\t" << it->second << "\t" << nameOfCaller << "\t" << nameOfSvType << std::endl;
  }
  outputFile.close();
}


int main(int argc, char** argv) {
  if (argc == 1) {
    std::cout << 
      "size_ass\n"
      "\n"
      "Purpose: returns a table summarizing the occurence of events of certain "
      "sizes in a VCF file. Its output is like\n"
      "\n"
      "Size    Count   Caller  SvType\n"
      "1       107774  Pindel  deletion\n"
      "2       36753   Pindel  deletion\n"
      "3       14359   Pindel  deletion\n"
      "...     ...     ...     ...\n"
      "\n"
      "Note that it should only be used on a file that contains only one type "
      "of event (insertion or deletion).\n"
      "\n"
      "Usage: ./size_ass input_vcf output_txt name_of_caller name_of_sv_type\n"
      "Example: ./size_ass pindel_hanchild_del.vcf pindel_hanchild_del_sizes.txt Pindel deletion\n"
      "\n"
      "Contact data: Eric-Wubbo Lameijer, Xi'an Jiaotong University, eric_wubbo@hotmail.com\n\n";
    return -1;
  } else if (argc < 5) {
    std::cout << "Invalid number of arguments. At least four arguments "
        "are needed, the name of the input file and the name of the "
        "output file, the name of the caller and the name of the SV type (like deletion).";
    return -1;
  } else {
    std::string nameOfInputFile = argv[1];
    std::string nameOfOutputFile = argv[2];
    std::string nameOfCaller = argv[3];
    std::string nameOfSvType = argv[4];
    transformFile(nameOfInputFile, nameOfOutputFile, nameOfCaller, nameOfSvType);

    return 0;
  }
}
