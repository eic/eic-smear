/**
 \file
 Implementation of class erhic::EventMilou.
 
 \author    Thomas Burton
 \date      2011-07-07
 \copyright 2011 Brookhaven National Lab
 */

#include "eicsmear/erhic/EventMilou.h"

#include <cmath>
#include <sstream>
#include <string>

namespace erhic {

EventMilou::EventMilou()
: radcorr(-1)
, weight(NAN)
, trueX(NAN)
, trueQ2(NAN)
, trueY(NAN)
, trueT(NAN)
, truePhi(NAN)
, phibelgen(NAN)
, phibelres(NAN)
, phibelrec(NAN) {
}

bool EventMilou::Parse(const std::string& line) {
  static std::stringstream ss;
  ss.str("");
  ss.clear();
  ss << line;
  ss >>
  number >> number >>  // Skip first int in the line
  nTracks >> weight >> process >> radcorr >> trueX >> trueQ2 >>
  trueY >> trueT >> truePhi >> phibelgen >> phibelres >> phibelrec;
  // Protect against errors in the input file or the stream
  return !ss.fail();
}

}  // namespace erhic
