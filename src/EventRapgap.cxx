//
// EventRapgap.cxx
//
// Created by TB on 6/10/11.
// Copyright 2011 BNL. All rights reserved.
//

#include <sstream>
#include <string>

#include "EventRapgap.h"

bool
EventRapgap::Parse(const std::string& line ) {
   
   static std::stringstream ss;
   
   ss.str("" );
   ss.clear();
   
   ss << line;
   ss >>
   number >> number >> // Skip first int in the line
   genevent >> process >> idir >> idisdif >> cs >> sigma_cs >> s >> q2 >>
   y >> xgam >> xpr >> Pt_h >> pt2_hat >> sHat >> t >> x_pom >> sHat2 >> z >>
   x1 >> phi1 >> nTracks;
   
   // Protect against errors in the input file or the stream
   return not ss.fail();
}
