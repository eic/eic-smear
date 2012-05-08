//
// EventMilou.cxx
//
// Created by TB on 6/10/11.
// Copyright 2011 BNL. All rights reserved.
//

#include <cmath>
#include <sstream>
#include <string>

#include "eicsmear/erhic/EventMilou.h"

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
   
   bool EventMilou::Parse(const std::string& line ) {
      static std::stringstream ss;
      ss.str("" );
      ss.clear();
      ss << line;
      ss >>
      number >> number >> // Skip first int in the line
      nTracks >> weight >> process >> radcorr >> trueX >> trueQ2 >>
      trueY >> trueT >> truePhi >> phibelgen >> phibelres >> phibelrec;
      // Protect against errors in the input file or the stream
      return not ss.fail();
   }
} // namespace erhic
