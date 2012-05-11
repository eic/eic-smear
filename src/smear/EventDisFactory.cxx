/**
 EventDisFactory.cxx
 
 \file
 Implementation of class EventDisFactory.
 
 \author Thomas Burton 
 \date 5/11/12
 \copyright 2012 BNL. All rights reserved.
 */

#include "eicsmear/smear/EventDisFactory.h"

#include <TBranch.h>

#include "eicsmear/erhic/EventDis.h"
#include "eicsmear/erhic/ParticleIdentifier.h"
#include "eicsmear/erhic/VirtualParticle.h"
#include "eicsmear/smear/ParticleMCS.h"

namespace {
   Smear::ParticleMCS* mcToSmear(const erhic::VirtualParticle& mc) {
      return new Smear::ParticleMCS(mc.Get4Vector(), mc.Id(), mc.GetStatus());
   }
}

namespace Smear {
   EventDisFactory::~EventDisFactory() {
   }
   EventDisFactory::EventDisFactory(const Detector& d, TBranch& mcBranch)
   : mDetector(d)
   , mMcEvent(NULL) {
      mcBranch.SetAddress(&mMcEvent);
   }
   Event* EventDisFactory::Create() {
      Event* event = new Event;
      ParticleIdentifier pid(mMcEvent->BeamLepton()->Id());
      for(unsigned j(0); j < mMcEvent->GetNTracks(); j++) {
         const erhic::VirtualParticle* ptr = mMcEvent->GetTrack(j);
         if(not ptr) {
            continue;
         } // if
         // It's convenient to keep the initial beams, unsmeared, in the
         // smeared event record, so copy their properties exactly
         if(mMcEvent->BeamLepton() == ptr or mMcEvent->BeamHadron() == ptr) {
            event->AddLast(mcToSmear(*ptr));
         } // if
         else {
            event->AddLast(mDetector.Smear(*ptr));
            // If this is the scattered lepton, record the index.
            // Check that the scattered lepton hasn't been set yet so we
            // don't replace it with a subsequent match.
            if(pid.isScatteredLepton(*ptr) and not event->ScatteredLepton()) {
               event->mScatteredIndex = j;
            } // if
         } // else
      } // for
        // Fill the event-wise kinematic variables.
      mDetector.FillEventKinematics(*mMcEvent, event);
      return event;
   }
} // namespace Smear
