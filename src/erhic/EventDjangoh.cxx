/**
 \file
 Implementation of class erhic::EventDjangoh.
 
 \author    Thomas Burton
 \date      2011-07-07
 \copyright 2011 Brookhaven National Lab
 */

#include "eicsmear/erhic/EventDjangoh.h"

#include <sstream>
#include <string>

namespace erhic {

bool EventDjangoh::Parse(const std::string& line) {
  static std::stringstream ss;
  ss.str("");
  ss.clear();
  ss << line;
  ss >>
  number >> number >>  // Skip first int in the line
  IChannel >> dprocess >> process >> nucleon >> dstruckparton >>
  dpartontrck >> dY >> dQ2 >> dX >> dW2 >> dNu >>
  dtrueY >> dtrueQ2 >> dtrueX >> dtrueW2 >> dtrueNu >>
  sigTot >> sigTotErr >> D >> F1NC >> F3NC >> G1NC >> G3NC >>
  A1NC >> F1CC >> F3CC >> G1CC >> G5CC >>
  nTracks;
  // Protect against errors in the input file or the stream
  return !ss.fail();
}

}  // namespace erhic
