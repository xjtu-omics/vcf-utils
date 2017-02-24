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
      //outputFile << line << "\n";
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
    int size = getEventSize(ref,alt);
    std::cout << "Size: " << size << ", mapsize: " << sizeCounts.size() << std::endl;
    if (sizeCounts.count(size) == 0 ) {
      sizeCounts[size] = 1;
    } else {
      ++sizeCounts[size];
    }
    /*if (isInsertion(ref,alt)) {
      //outputFile << line << "\n";
    } else if (isDeletion(ref,alt)) {
      
      outputFile << line << "\n";
    } else {
      std::cout << "alarm: " << line << "\n";
    }*/

  }     
  inputFile.close();
  outputFile << "Size\tCount\tCaller\tSvType" << std::endl;
  for (std::map<int,int>::iterator it = sizeCounts.begin(); it != sizeCounts.end(); ++it) {
    outputFile << it->first << "\t" << it->second << "\t" << nameOfCaller << "\t" << nameOfSvType << std::endl;
  }
  outputFile.close();
}


int main(int argc, char** argv) {
  //std::cout << "Converting the input VCF to output VCF";
  if (argc < 5) {
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