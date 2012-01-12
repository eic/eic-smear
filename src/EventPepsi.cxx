//
// EventPepsi.cxx
//
// Created by TB on 6/10/11.
// Copyright 2011 BNL. All rights reserved.
//

#include <sstream>
#include <string>

#include "EventPepsi.h"

bool
EventPepsi::Parse(const std::string& line ) {
   
   static std::stringstream ss;
   
   ss.str("" );
   ss.clear();
   
   ss << line;
   ss >>
   number >> number >> // Skip first int in the line
   genevent >> process >> subprocess >> nucleon >> struckparton >> partontrck >>
   trueY >> trueQ2 >> trueX >> trueW2 >> trueNu >>
   FixedWeight >> Weight >> dxsec >> ExtraWeight >>
   Dilute >> F1 >> F2 >> A1 >> A2 >> R >> DePol >> D >> Eta >> Eps >> Chi >>
   gendilut >> genF1 >> genF2 >> genA1 >> genA2 >> genR >> genDepol >> gend >> geneta >> geneps >> genchi >>
   SigCorr >> radgamEnucl >> nTracks;
   
   // Protect against errors in the input file or the stream
   return not ss.fail();
}
