/**
 \file
 Implementation of class erhic::EventRapgap.

 \author    Thomas Burton
 \date      2011-07-07
 \copyright 2011 Brookhaven National Lab
 */

#include "eicsmear/erhic/EventRapgap.h"

#include <sstream>
#include <string>

namespace erhic {

bool EventRapgap::Parse(const std::string& line) {
  static std::stringstream ss;
  ss.str("");
  ss.clear();
  ss << line;
  ss >>
  number >> number >>  // Skip first int in the line
  genevent >> process >> idir >> idisdif >> cs >> sigma_cs >> s >> q2 >>
  y >> xgam >> xpr >> Pt_h >> pt2_hat >> sHat >> t >> x_pom >> sHat2 >> z >>
  x1 >> phi1 >> nTracks;
  // Protect against errors in the input file or the stream
  return !ss.fail();
}

}  // namespace erhic
