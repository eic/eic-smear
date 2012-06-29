/**
 Pythia6EventBuilder.h

 \file
 Declaration of class Pythia6EventBuilder.

 Created by TB on 1/17/12.
 Copyright 2012 BNL. All rights reserved.
*/

#ifndef _Pythia6EventBuilder_H_
#define _Pythia6EventBuilder_H_

#include <string>

#include <Rtypes.h> // For ClassDef macro

#include "eicsmear/erhic/EventFactory.h"

namespace erhic {

   class EventMCFilterABC;
   class EventPythia;
   
   /**
    Interface to PYTHIA 6.
    Builds EventPythia objects directly from PYTHIA output, without
    an intermediate text file.
    Configure PYTHIA options via the ROOT::TPythia6 class.
    \remark The following fields in EventPythia are not set:
    F1, F2, R, sigma_rad, SigRadCor, EBrems
    */
   class Pythia6EventBuilder : public VirtualEventFactory {

   public:
      
      /**
       Constructor.
       If a filter is provided, it will be applied so that all events
       yielded by Create() pass the filter.
       The filter should be allocated via new and is subsequently
       owned and deleted by the Pythia6EventBuilder.
       */
      Pythia6EventBuilder(EventMCFilterABC* = NULL);
      
      /** Destructor */
      virtual ~Pythia6EventBuilder();
      
      /**
       Generates an event from the current state of ROOT::TPythia6.
       The returned event is dynamically allocated and must
       be deleted by the user.
       */
      virtual EventPythia* Create();

      virtual std::string EventName() const;

      virtual TBranch* Branch(TTree&, const std::string&);
      virtual void Fill(TBranch&);
   protected:
      EventPythia* BuildEvent();
      Long64_t mNGenerated;
      Long64_t mNTrials;
      EventMCFilterABC* mFilter;
      EventPythia* mEvent;
      ClassDef(erhic::Pythia6EventBuilder, 1)
   };
   
} // namespace erhic

#endif
