//
// EventDjangoh.cxx
//
// Created by TB on 6/10/11.
// Copyright 2011 BNL. All rights reserved.
//

#include <sstream>
#include <string>

#include "eicsmear/erhic/EventDjangoh.h"

namespace erhic {
   
   bool EventDjangoh::Parse(const std::string& line ) {
      static std::stringstream ss;
      ss.str("" );
      ss.clear();
      ss << line;
      ss >>
      number >> number >> // Skip first int in the line
      IChannel >> dprocess >> process >> nucleon >> dstruckparton >>
      dpartontrck >> dY >> dQ2 >> dX >> dW2 >> dNu >>
      dtrueY >> dtrueQ2 >> dtrueX >> dtrueW2 >> dtrueNu >>
      nTracks;
      // Protect against errors in the input file or the stream
      return not ss.fail();
   }
} // namespace erhic