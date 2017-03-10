/**
 * vcf_standardizer.cpp
 *
 * Purpose: aligns indels in a VCF file to the leftmost position (not all SV-callers do so). Also, if
 * a VCF file is in a "PacBio format" (not TA T but T <DEL> SV=A it will convert it into the TA T format
 * to make comparison with other VCFs more straightforward)
 *
 * Usage: ./standardize input_vcf reference_fasta output_vcf
 * Example: ./standardize pacbio_hanchild.vcf hg38.fa pacbio_hanchild_leftaligned.vcf
 * 
 * Contact data: Eric-Wubbo Lameijer, Xi'an Jiaotong University, eric_wubbo@hotmail.com
 */

#include <cctype> // toupper
#include <cstdlib> // atoi
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

/** Utility function that halts/crashes the program, helps to catch bugs early.
**/
void Require(bool requirementMet, std::string errorMessage) {
  if (!requirementMet) {
    std::cerr << errorMessage << std::endl;
    exit(-1);
  }
}

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

std::string join(std::vector<std::string> stringVector, char ch) {
  std::string result = "";
  int numElements = stringVector.size();
  for (int i = 0; i < numElements - 1; i++ ) {
    result += stringVector[i] + ch;
  }
  if (numElements > 0 ) {
    result += stringVector[numElements - 1]; 
  }
  return result;
}

std::string charToString(char ch) {
  std::stringstream ss;
  ss << ch;
  std::string str;
  ss >> str;
  return str;
}

std::string intToString(int i) {
  std::stringstream ss;
  ss << i;
  std::string str;
  ss >> str;
  return str;
}

/**
 * 'Event' represents a genetic event (so insertion or deletion).
 * 
 * @author Eric-Wubbo Lameijer, Xi'an Jiaotong University, eric_wubbo@hotmail.com
 *
 */
class Event {

 private:	
  std::vector<std::string> allFields_; // saves all fields of the VCF line, handy for output
  std::string alternativeAllele_; // the alt allele/alternative allele: the sequence of DNA in the individual that is other than
		// expected by consulting the reference chromosome.
  std::string chromosomeName_; // the name of the chromosome in which this event takes place
  int position_; // position of the event
  std::string referenceAllele_; // the expected sequence at this place in this chromosome, as indicated in the reference genome
  bool isPacBio_; // PacBio VCF files have a different format, like T <DEL> instead of TTCG T

 public:	
  /**
  * Event constructor; creates an event out of a line of a VCF file.
  * 
  * @param vcfLine
  * 		The line of the VCF file to be turned into an event.
  */
  Event(const std::string& vcfLine, const std::string& reference) {
    allFields_ = Split(vcfLine,'\t');
    chromosomeName_ = allFields_[0];
    position_ = atoi(allFields_[1].c_str());
    referenceAllele_ = allFields_[3];
    alternativeAllele_ = allFields_[4];
    //std::cout << "\nRef: " << referenceAllele_ << ", alt: " << alternativeAllele_ << std::endl;
    isPacBio_ = ((alternativeAllele_ == "<INS>") || (alternativeAllele_ == "<DEL>"));
    if (isPacBio_) {
	
	   std::string info = allFields_[7];
      int svLenPosition = info.find("SVLEN=") + 6;
		int svLenEndCandidate = svLenPosition;
		while (info[svLenEndCandidate] >= '0' && info[svLenEndCandidate] <= '9') {
        ++svLenEndCandidate;
		}
		int svLenLength = svLenEndCandidate - svLenPosition;
		std::string svLenAsString = info.substr(svLenPosition, svLenLength);
		int svLength = atoi(svLenAsString.c_str());

		int sequencePosition = info.find("SEQ=") + 4;
	   char base = toupper(info[sequencePosition]);
  
		
		if (alternativeAllele_ == "<INS>") {
        std::cout << "Ref=" << referenceAllele_ << ":" << reference[position_] << std::endl;
        /*if (base != referenceAllele_[0] || base != reference[position_]) {
          std::cout << "Abort: weird base!" << "Base:" << base 
						<< "Ref: " << referenceAllele_[0] <<
						"chromosome: " << reference[position_] 
						<< "Position: " << position_ << std::endl;
			 exit(-1);
		  }*/
		  // sometimes the alt starts with a different base than the ref...
		  alternativeAllele_ = referenceAllele_;
	     while (isValidDnaBase(base)) {
	       alternativeAllele_ += toupper(base);
	       ++sequencePosition;
	       base = info[sequencePosition];
	     }
		  std::cout << "I_alt len=" << alternativeAllele_.length() <<
					"seq=" << referenceAllele_ << "/" << 
					alternativeAllele_ << "; info=" << info << std::endl << std::endl;
		  char ch;
       // std::cin >> ch;  

      } else if (alternativeAllele_ == "<DEL>") {
	     alternativeAllele_ = referenceAllele_;
	     for (int i = 1; i <= svLength; ++i ) {
	       referenceAllele_ += toupper(reference[position_ + i]);
	     }
		  std::cout << "D_ref len=" << referenceAllele_.length() << "; info=" << info << std::endl << std::endl;
		char ch;
      //std::cin >> ch;  
      } else {
        Require(false, "Event constructor error: event type " + alternativeAllele_ + " is unknown.");
      }
    } // if isPacBio
  }

  /**
  * Is this a valid base of DNA?
  * 
  * @param base the base (or character) to be checked
  * @return whether the base is a valid DNA base
  */
  bool isValidDnaBase(char base) {
    char baseInUpperCase = toupper(base);
    return (baseInUpperCase == 'A' || baseInUpperCase == 'C' || baseInUpperCase == 'G' || baseInUpperCase == 'T');
  }

  /**
   * Returns the name of the chromosome in which this event takes place.
   * 
   * @return the name of the chromosome in which this event takes place.
   */
  std::string getChromosome() {
    return chromosomeName_;
  }
	
  /**
   * Does this event have multiple alternative alleles?
   * 
   * @return whether the event has multiple alternative alleles.
   */
  bool hasMultipleAltAlleles() {
    //std::cout << "Alt allele: " << alternativeAllele_ << std::endl;
    return alternativeAllele_.find(",") != std::string::npos;
  }
	
  /**
   * Returns whether the event is a deletion.
   * 
   * @return whether the event is a deletion
   */
  bool isDeletion() {
    //std::cout << "Z3\n";
    return (referenceAllele_.length() > 1 && alternativeAllele_.length() == 1);
  }
	
  /**
   * Returns whether the event is an insertion.
   * 
   * @return whether the event is an insertion
  */
  bool isInsertion() {
    //std::cout << "Z2\n";
    return (referenceAllele_.length() == 1 && alternativeAllele_.length() > 1 && !hasMultipleAltAlleles());
  }


  /**
   * Left-aligns the event
   * 
   * @param sequenceOfCurrentChromosome
   */
  void leftAlign(const std::string& sequenceOfCurrentChromosome) {
    //std::cout << "Z1\n";
    if (isInsertion()) {
      //std::cout << "A1\n";
      Require(referenceAllele_[0] == sequenceOfCurrentChromosome[position_], 
	"Event.leftAlign error: reference chromosome/pos " + chromosomeName_ + ":" + intToString(position_) +
        "[" + sequenceOfCurrentChromosome[position_] + "] and reference sequence [" + referenceAllele_[0] + "] don't seem to match.");
      char referenceBase = referenceAllele_[0];
      //std::cout << "A2\n";
      while (alternativeAllele_[alternativeAllele_.length() - 1] == referenceBase && position_ > 1) {
	std::cout << "shifting " << chromosomeName_ << position_ << std::endl;
	--position_;
	referenceBase = sequenceOfCurrentChromosome[position_];
	referenceAllele_ = charToString(referenceBase);
	alternativeAllele_ = referenceAllele_ + alternativeAllele_.substr(0, alternativeAllele_.length() - 1);
      }
    } else if (isDeletion()){
      //std::cout << "B1\n";
      Require(referenceAllele_[0] == sequenceOfCurrentChromosome[position_], 
	"Event.leftAlign error: reference chromosome and reference sequence don't seem to match.");
      int eventLength = referenceAllele_.length() - alternativeAllele_.length();
      int lastPositionOfDeletion = position_ + eventLength;
      char referenceBase = referenceAllele_[0];
      //std::cout << "B2\n";
      while (sequenceOfCurrentChromosome[lastPositionOfDeletion] == referenceBase && position_ > 1) {
        std::cout << "shifting " << chromosomeName_ << position_ << std::endl;
	--position_;
	referenceBase = sequenceOfCurrentChromosome[position_];
	referenceAllele_ = referenceBase + referenceAllele_.substr(0, referenceAllele_.length() - 1);
        alternativeAllele_ = charToString(referenceBase);
	lastPositionOfDeletion = position_ + eventLength;			
      }		
    } else {
      std::cout << "Can't handle " << asLine();
    }
  }

  /** 
   * Returns the event as a line.
   * 
   * @return the event as a line/String.
   */
   std::string asLine() {
     //std::cout << "start\n";
     std::stringstream ss;
     ss << position_;
     ss >> allFields_[1];
//std::cout << "start21\n";
     allFields_[3] = referenceAllele_;
     allFields_[4] = alternativeAllele_;
     /*if (position_ == 90384281) {
       std::cout << "Ref: " << referenceAllele_ << "Alt: " <<alternativeAllele_ << std::endl;
       exit(0);
     }*/
//std::cout << "start2\n";
    /* if (isPacBio_) {
       if (referenceAllele_.length() > alternativeAllele_.length()) {
	 // deletion
	 allFields_[3] = alternativeAllele_;
	 allFields_[4]  = "<DEL>";
       } else {
	 // insertion
	 allFields_[3] = referenceAllele_;
	 allFields_[4] = "<INS>";
       }
     }*/
//std::cout << "star3t\n";
     return join(allFields_, '\t');
   }
};

const std::string load(const std::string nameOfChromosome, std::ifstream& genomeFile, std::string& outputSequence);

int main(int argc, char** argv) {
  if (argc != 4) {
    std::cout << 
      "standardize\n"
      "\n"
      "Purpose: aligns indels in a VCF file to the leftmost position (not all SV-callers do so). Also, if "
      "a VCF file is in a \"PacBio format\" (not TA T but T <DEL>) it will convert it into the TA T format "
      "to make comparison with other VCFs more straightforward).\n"
      "\n"
      "Usage: ./standardize input_vcf reference_fasta output_vcf\n"
      "Example: ./standardize pacbio_hanchild.vcf hg38.fa pacbio_hanchild_leftaligned.vcf\n"
      "\n"
      "contact data: Eric-Wubbo Lameijer, Xi'an Jiaotong University, eric_wubbo@hotmail.com\n\n";
    return -1;
  }   
  std::string nameOfInputVcf = argv[1];
  std::string nameOfReference = argv[2];
  std::string nameOfOutputVcf = argv[3];

  std::cout << "input VCF: " << nameOfInputVcf << std::endl;

  std::string nameOfCurrentChromosome = "";
  std::string sequenceOfCurrentChromosome = "";
  
  std::string line = "";
		
  std::ifstream inputVcf(nameOfInputVcf.c_str());
  std::ifstream referenceGenome(nameOfReference.c_str());
  std::ofstream outputVcf(nameOfOutputVcf.c_str());

  while (!inputVcf.eof()) {
    getline(inputVcf, line);
    //std::cout << line << std::endl;
    if (line.length() == 0) {
      break; // we're done
    }
    if (line[0] == '#') {
      // this is a comment line, copy comment lines directly to the output
      outputVcf << line << std::endl;

    } else {
		  // apparently, we've reached the first event
		std::stringstream ss;
		ss << line;
      std::string chromosomeOfEvent; 
		ss >> chromosomeOfEvent;

      if (chromosomeOfEvent != nameOfCurrentChromosome) {
        load(chromosomeOfEvent, referenceGenome, sequenceOfCurrentChromosome);
        nameOfCurrentChromosome = chromosomeOfEvent;

        // the below deals with mismatching reference genome - VCF chromosome orders
        if (sequenceOfCurrentChromosome.length() == 0) {
	  referenceGenome.close();
          referenceGenome.open(nameOfReference.c_str());
          sequenceOfCurrentChromosome = load(chromosomeOfEvent, referenceGenome, sequenceOfCurrentChromosome);
        }
      }
	   Event event(line, sequenceOfCurrentChromosome);
      //std::cout << "Ready with chromosome loading" << std::endl;
      event.leftAlign(sequenceOfCurrentChromosome);
      //std::cout << "Ready with aligning" << std::endl;
      outputVcf << event.asLine() << std::endl;
//std::cout << "Ready with outputting" << std::endl;
    }
  }
}


/**
 * Loads a chromosome with the specified name into memory.
 * @param nameOfChromosome
 * 		The name of the chromosome which is sought
 * @param genomeFile
 * 		The file which houses the reference genome
 * @return
 * 		The sequence of the sought chromosome.
 */
const std::string load(const std::string nameOfChromosome, std::ifstream& genomeFile, std::string& outputSequence) {
  outputSequence = "";
  //std::cout << "Loading chromosome " << nameOfChromosome << std::endl;
  std::string soughtStartSequence = ">" + nameOfChromosome;
  int startSequenceLength = soughtStartSequence.length();
  std::string line;
  int counter = 1;
  bool copyLines = false;
  std::streampos lastReadPosInFile;
  while (!genomeFile.eof()) {
    lastReadPosInFile = genomeFile.tellg();
    getline(genomeFile, line);
    //std::cout << line << std::endl;
    if (line.length() == 0) {
      break;
    }
    if (line[0] == '>') {
      std::cout << line << std::endl;
      // two cases: either turn copying on, or turn it off.
      if (StringStartsWith(line,soughtStartSequence) && isspace(line[startSequenceLength])) {
	copyLines = true;
	outputSequence += "N"; // for easy conversion of C++ coordinates to 1-based coordinates of most
					 // references
      } else { // >otherChromosomeName
	if (copyLines) {
	  genomeFile.seekg(lastReadPosInFile, genomeFile.beg);
	  return outputSequence; // done
	}
      }
    } else {
      // is just a normal chromosome sequence
      if (copyLines) {
	outputSequence += line;
      }
    }			
  }
	//Utilities.require(chromosomeSequence.length() > 0 , "VcfLeftAligner.load error: sought chromosome " + nameOfChromosome + " not found.");
  return outputSequence;
}

