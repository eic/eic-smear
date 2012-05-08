/**
 smear/EventHadronic.h
 
 \file Declaration of class smear::EventHadronic
 
 \author TB
 \date 5/3/12
 \copyright 2012 BNL. All rights reserved.
 */

#ifndef _EICSMEAR_SMEAR_EVENTHADRONIC_H_
#define _EICSMEAR_SMEAR_EVENTHADRONIC_H_

#include <vector>

#include <Rtypes.h> // For ClassDef macro

#include <eicsmear/hadronic/EventMC.h>
#include "eicsmear/smear/ParticleMCS.h"
#include "eicsmear/erhic/VirtualEvent.h"
#include "eicsmear/erhic/EventFactory.h"

#include <TBranch.h>

namespace erhic {
   namespace hadronic {
      /**
       Realisation of hadronic::EventMC as an event with detector smearing.
       */
      class EventSmear : public VirtualEvent {
      public:
         virtual ~EventSmear();
         EventSmear();
         virtual const Smear::ParticleMCS* GetTrack(UInt_t) const;
         virtual Smear::ParticleMCS* GetTrack(UInt_t);
         virtual UInt_t GetNTracks() const;
         virtual void AddLast(Smear::ParticleMCS*);
      protected:
         std::vector<Smear::ParticleMCS*> particles;
         ClassDef(erhic::hadronic::EventSmear, 1)
      };
   } // namespace hadronic
} // namespace erhic

namespace Smear {
   /**
    Factory class for smeared hadronic events.
   */
   class HadronicEventBuilder :
   public EventFactory<erhic::hadronic::EventSmear> {
   public:
      /** Destructor */
      virtual ~HadronicEventBuilder() {
      }
      /**
       Constructor.
       Initialise with the Detector performing particle smearing and
       the tree branch providing DIS Monte Carlo events.
      */
      HadronicEventBuilder(const Detector& d, TBranch& mcBranch)
      : mDetector(d)
      , mMcEvent(NULL) {
         mcBranch.SetAddress(&mMcEvent);
      }
      /**
       Create a smeared event corresponding to the current DIS Monte Carlo
       event in the input branch passed to the constructor.
       The user should call TTree::GetEntry() themselves between calls to
       Create().
      */
      virtual EventType* Create() {
         EventType* event = new EventType;
         for(unsigned j(0); j < mMcEvent->GetNTracks(); j++) {
            const erhic::VirtualParticle* ptr = mMcEvent->GetTrack(j);
            if(ptr) {
               event->AddLast(mDetector.Smear(*ptr));
            } // if
         } // for
         return event;
      }
   protected:
      Detector mDetector;
      erhic::hadronic::EventMC* mMcEvent;
   };
} // namespace Smear

#endif
