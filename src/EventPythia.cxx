//
// EventPythia.cxx
//
// Created by TB on 6/10/11.
// Copyright 2011 BNL. All rights reserved.
//

#include <cmath>
#include <sstream>
#include <string>

#include "EventPythia.h"

namespace erhic {
   
   EventPythia::EventPythia(const std::string& /* Unused */ )
   : nucleon(std::numeric_limits<Int_t>::max() )
   , tgtparton(std::numeric_limits<Int_t>::max() )
   , beamparton(std::numeric_limits<Int_t>::max() )
   , genevent(-1 )
   , xtgtparton(NAN)
   , xbeamparton(NAN)
   , thetabeamparton(NAN)
   , leptonphi(NAN)
   , F1(NAN)
   , sigma_rad(NAN)
   , t_hat(NAN)
   , u_hat(NAN)
   , Q2_hat(NAN)
   , SigRadCor(NAN)
   , EBrems(NAN)
   , photonflux(NAN)
   , trueY(NAN)
   , trueQ2(NAN)
   , trueX(NAN)
   , trueW2(NAN)
   , trueNu(NAN)
   , F2(NAN)
   , R(NAN)
   , pt2_hat(NAN)
   , sHat(NAN)
   {
   }
   
   
   EventPythia::~EventPythia() { }
   
   
   bool
   EventPythia::Parse(const std::string& line ) {
      
      static std::stringstream ss;
      
      ss.str("" );
      ss.clear();
      
      ss << line;
      ss >>
      number >> number >> // Skip first int in the line
      genevent >> process >> nucleon >> tgtparton >> xtgtparton >>
      beamparton >> xbeamparton >> thetabeamparton >> trueY >> trueQ2 >>
      trueX >> trueW2 >> trueNu >> leptonphi >> sHat >> t_hat >> u_hat >>
      pt2_hat >> Q2_hat >> F2 >> F1 >> R >> sigma_rad >> SigRadCor >> EBrems >>
      photonflux >> nTracks;
      
      // Protect against errors in the input file or the stream
      return not ss.fail();
   }
   
} // namespace erhic
