/** 
  read_reference
    
  Purpose: Finds a certain region in the reference genome.
  
  Usage: ./read_reference chromosome start_pos end_pos reference_fasta
  Example: ./read_reference chr1 1020 1040 hg38.fa

  Contact data: Eric-Wubbo Lameijer, Xi'an Jiaotong University, eric_wubbo@hotmail.com
**/
  
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

const std::string load(const std::string nameOfChromosome, std::ifstream& genomeFile, std::string& outputSequence);

void finalizeRegion(const std::string& nameOfCurrentChromosome, bool isNRegion, int lastSwapPosition, int newSwapPosition) {
  if (nameOfCurrentChromosome != "") {
    std::cout << nameOfCurrentChromosome << "\t" << lastSwapPosition << "\t" << newSwapPosition << "\t";
    std::cout << (isNRegion ? "unknown" : "genomic") << std::endl;
  }
}

std::string intToString(int i) {
  std::stringstream ss;
  ss << i;
  std::string str;
  ss >> str;
  return str;
}

std::string getChromosomeName(const std::string& line) {
  std::stringstream ss;
  ss << line;
  std::string chromosomeName;
  ss >> chromosomeName;
  chromosomeName = chromosomeName.substr(1); // skip '>'
  return chromosomeName;
}

bool isConsistent(char ch, bool isNRegion) {
  char upCh = toupper(ch);
  if (ch == 'N') {
    return isNRegion;
  } else if (ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T') {
    return !isNRegion;
  } else {
    //std::cout << "Unknown character " << ch << std::endl;
    return true; // no phase swap used
  } 
}

int main(int argc, char** argv) {
  if (argc == 1) {
    std::cout << 
      "read_reference\n"
      "\n"
      "Purpose: Finds a certain region in the reference genome.\n"
      "\n"
      "Usage: ./read_reference chromosome start_pos end_pos reference_fasta\n"
      "Example: ./read_reference chr1 1020 1040 hg38.fa\n"
      "\n"
      "Contact data: Eric-Wubbo Lameijer, Xi'an Jiaotong University, eric_wubbo@hotmail.com\n\n";
    return -1;
  } else if (argc != 5) {
    std::cout << "read_reference error: four arguments are required, the name of the chromosome, "
      "the start and end positions, and the name of the reference (fasta) file.";
    return -1;
  }
  
  const std::string nameOfTargetChromosome = argv[1];
  int startPos = atoi(argv[2]);
  int endPos = atoi(argv[3]);
  std::string nameOfReference = argv[4];
  std::string sequenceOfCurrentChromosome = "";
  std::string line = "";
  //std::cout << "S:" << startPos << ", E: " << endPos << std::endl;
  //char ch;
  //std::cin >> ch;
		
  std::ifstream genomeFile(nameOfReference.c_str());

  std::string nameOfCurrentChromosome = "";

  bool outputBase = false;
  std::string startPosAsString = argv[2];
  std::string bottomline = startPosAsString + ":";
  std::string topline = "";
  for (int i = 0; i < bottomline.length(); i++ ) {
    topline += " ";
  }
  while (!genomeFile.eof()) {
    getline(genomeFile, line);
   
    if (line.length() == 0) {
      break;
    }
    if (line[0] == '>') {
      // a new chromosome!
      //std::cout << line << std::endl;
      nameOfCurrentChromosome = getChromosomeName(line);
      //std::cout << "CC: [" << nameOfCurrentChromosome <<"], TC:[" << nameOfTargetChromosome <<"]" << (nameOfCurrentChromosome == nameOfTargetChromosome? "Y" : "N") << std::endl;
//char ch;
  //std::cin >> ch;
      sequenceOfCurrentChromosome = "N"; // to make C++ coordinates correspond to biological coordinates
    } else {
      // this is a regular line
      
      for (int i = 0; i < line.length(); ++i ) {
        if (nameOfCurrentChromosome == nameOfTargetChromosome) {
          //std::cout << sequenceOfCurrentChromosome.length() << ":" << line << std::endl;
          int position = sequenceOfCurrentChromosome.length();
          if (position == startPos) {
            outputBase = true;
          } else if (position == endPos) {
            outputBase = false;
            std::cout << topline << std::endl;
            std::cout << bottomline << std::endl;
            exit(0);
          }
          char newBase = toupper(line[i]);  
          sequenceOfCurrentChromosome += newBase;
          if (outputBase) {
            topline += intToString(position % 10);
            bottomline += newBase;
            if (position % 5 == 0) {
              topline += " ";
              bottomline += " ";
            }
          }  
        } 
      } // loop over all bases
    } // process a line
  } // while there is still data in the reference
}
