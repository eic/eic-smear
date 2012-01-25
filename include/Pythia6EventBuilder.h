/**
 Pythia6EventBuilder.h

 \file
 Declaration of class Pythia6EventBuilder.

 Created by TB on 1/17/12.
 Copyright 2012 BNL. All rights reserved.
*/

#ifndef _Pythia6EventBuilder_H_
#define _Pythia6EventBuilder_H_

#include <Rtypes.h> // For ClassDef macro

#include "EventPythia.h"

namespace erhic {
   
   
   /**
    Interface to PYTHIA 6.
    Builds EventPythia objects directly from PYTHIA output, without
    an intermediate text file.
    Configure PYTHIA options via the ROOT TPythia6 class.
    \remark The following fields in EventPythia are not set:
    F1, F2, R, sigma_rad, SigRadCor, EBrems
    */
   class Pythia6EventBuilder {
      
   public:
      
      /** Constructor */
      Pythia6EventBuilder();
      
      /** Destructor */
      virtual ~Pythia6EventBuilder();
      
      /**
       Generates an event from the current state of TPythia6.
       The returned event is dynamically allocated and must
       be deleted by the user.
       */
      virtual EventPythia* Create();
      
      /** */
//      virtual int NEvents() const { return mNEvents; }
      
      /** */
//      virtual int NTrials() const { return mNTrials; }
      
   protected:
      
//      int mNEvents;
//      int mNTrials;
      
      ClassDef(erhic::Pythia6EventBuilder, 1)
   };
   
} // namespace erhic

#endif
