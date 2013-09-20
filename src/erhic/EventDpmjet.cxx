/**
 \file
 Implementation of class erhic::EventDpmjet.
 
 \author    Thomas Burton
 \date      2011-08-19
 \copyright 2011 Brookhaven National Lab
 */

#include "eicsmear/erhic/EventDpmjet.h"

#include <sstream>
#include <string>

namespace erhic {

bool EventDpmjet::Parse(const std::string& line) {
  static std::stringstream ss;
  ss.str("");
  ss.clear();
  ss << line;
  ss >> I >> ievent >> process1 >> process2 >> IP >> dtrueW2
  >> dtrueNu >> dtrueQ2 >> dtrueX >> dtrueY >> theta_Evt >> photonFlux
  >> tgtparton >> prjparton >> xtgtparton >> xprjparton >> nucleon
  >> nTracks;
  number = ievent;
  return !ss.fail();
}

}  // namespace erhic
