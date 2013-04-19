/**
 EventDisFactory.h

 \file
 Declaration of class EventDisFactory.

 \author Thomas Burton 
 \date 5/11/12
 \copyright 2012 BNL. All rights reserved.
*/

#ifndef _EICSMEAR_SMEAR_EVENTDISBUILDER_H_
#define _EICSMEAR_SMEAR_EVENTDISBUILDER_H_

#include "eicsmear/smear/Detector.h"
#include "eicsmear/smear/EventSmear.h"
#include "eicsmear/smear/EventFactory.h"

class TBranch;

namespace erhic {
   class EventDis;
} // namespace erhic

namespace Smear {
   /**
    Factory class for smeared DIS events.
   */
   class EventDisFactory : public EventFactory<Smear::Event> {
   public:
      /** Destructor */
      virtual ~EventDisFactory();
      /**
       Constructor.
       Initialise with the Detector performing particle smearing and
       the tree branch providing DIS Monte Carlo events.
      */
      EventDisFactory(const Detector&, TBranch&);
      /**
       Create a smeared event corresponding to the current DIS Monte Carlo
       event in the input branch passed to the constructor.
       The user should call TTree::GetEntry() on the tree with the
       input branch between calls to Create().
      */
      virtual Event* Create();

      erhic::VirtualEvent* GetEvBufferPtr() { return mMcEvent;} ;
   protected:
      Detector mDetector;
      erhic::EventDis* mMcEvent;
   };
} // namespace Smear

#endif
