/**
 \file Implementation of class erhic::EventSimple. 
 \author    Barak Schmookler
 \date      2020-03-19
 */

#include "eicsmear/erhic/EventSimple.h"

#include <cmath>
#include <sstream>
#include <string>

namespace erhic {

EventSimple::EventSimple()
: numParticles(NAN){
}

bool EventSimple::Parse(const std::string& line) {
  static std::stringstream ss;
  ss.str("");
  ss.clear();
  ss << line;
  ss >>
  number >> number >>  // Skip first int in the line
    numParticles;
  // Protect against errors in the input file or the stream
  return !ss.fail();
}

}  // namespace erhic
