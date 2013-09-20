/**
 \file
 Implementation of class erhic::EventPepsi.
 
 \author    Thomas Burton
 \date      2011-07-07
 \copyright 2011 Brookhaven National Lab
 */

#include "eicsmear/erhic/EventPepsi.h"

#include <sstream>
#include <string>

namespace erhic {

bool EventPepsi::Parse(const std::string& line) {
  static std::stringstream ss;
  ss.str("");
  ss.clear();
  ss << line;
  ss >>
  number >> number >>  // Skip first int in the line
  genevent >> process >> subprocess >> nucleon >> struckparton >>
  partontrck >> trueY >> trueQ2 >> trueX >> trueW2 >> trueNu >>
  FixedWeight >> Weight >> dxsec >> ExtraWeight >>
  Dilute >> F1 >> F2 >> A1 >> A2 >> R >> DePol >> D >> Eta >> Eps >> Chi >>
  gendilut >> genF1 >> genF2 >> genA1 >> genA2 >> genR >> genDepol >>
  gend >> geneta >> geneps >> genchi >> SigCorr >> radgamEnucl >> nTracks;
  // Protect against errors in the input file or the stream
  return !ss.fail();
}

}  // namespace erhic
