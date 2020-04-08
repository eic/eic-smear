/**
 \file Implementation of class erhic::EventSartre. 
 \author    Barak Schmookler
 \date      2020-04-08
 */

#include "eicsmear/erhic/EventSartre.h"

#include <cmath>
#include <sstream>
#include <string>

namespace erhic {

EventSartre::EventSartre()
: genevent(-1)
, trueT(NAN)
, trueQ2(NAN)
, trueX(NAN)
, trueY(NAN)
, trueW2(NAN)
, trueNu(NAN)
, trueXpom(NAN)
, s_cm(NAN)
, pol(-1)
, dmode(-1)
, bup(-1){
}

bool EventSartre::Parse(const std::string& line) {
  static std::stringstream ss;
  ss.str("");
  ss.clear();
  ss << line;
  ss >>
  number >> number >>  // Skip first int in the line
  genevent >>
  trueT >> trueQ2 >> trueX >> trueY >> trueW2 >> trueNu >>
  trueXpom >> s_cm >> pol >> dmode >> bup;
  // Protect against errors in the input file or the stream
  return !ss.fail();
}

}  // namespace erhic
