//
// Pid.cxx
//
// Created by TB on 8/12/11.
// Copyright 2011 BNL. All rights reserved.
//

#include <exception>
#include <iostream>
#include <limits>

#include <TDatabasePDG.h>

#include "Pid.h"

namespace erhic {

   Pid::Pid(Int_t code)
   : mCode(code)
   { }

   Pid::~Pid()
   { }
   
   /**
    TODO May need to cache TParticlePDG* if performance proves too slow with
    database lookup upon each call.
    */
   TParticlePDG* Pid::Info() const {
      try {
         return TDatabasePDG::Instance()->GetParticle(Code());
      } // try
      catch(std::exception& e) {
         std::cerr
         << "Caught exception in Pid::Info(): " << e.what()
         << std::endl;
         return NULL;
      } // catch
   }
   
   // CINT can't handle <limits> in header files so we hide it away here.
   Int_t Pid::InvalidCode() {
      return std::numeric_limits<Int_t>::max();
   }
   
} // namespace
