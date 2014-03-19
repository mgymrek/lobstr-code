#include <sstream>
#include <string>

#include <stdlib.h>

#include "src/tests/DNATools.h"

namespace DNATools {
  char GetChar(int index){
    switch(index){
  case 0:
    return 'A';
    case 1:
      return 'C';
    case 2:
      return 'G';
    case 3:
      return 'T';
    default:
      return 'N';
    }
  }

  std::string RandDNA(const int length){
    std::stringstream seq;
    for (int j = 0; j < length; j++)
      seq << GetChar(rand()%4);
    return seq.str();
  }
}

