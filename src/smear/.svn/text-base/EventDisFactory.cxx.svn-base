/**
 \file
 Implementation of class Smear::EventDisFactory.
 
 \author    Michael Savastio
 \date      2011-08-19
 \copyright 2011 Brookhaven National Lab
 */

#include "eicsmear/smear/EventDisFactory.h"

#include <vector>

#include <TBranch.h>

#include "eicsmear/erhic/EventDis.h"
#include "eicsmear/erhic/ParticleIdentifier.h"
#include "eicsmear/erhic/VirtualParticle.h"
#include "eicsmear/smear/ParticleMCS.h"

namespace {

Smear::ParticleMCS* mcToSmear(const erhic::VirtualParticle& mc) {
  Smear::ParticleMCS* p = new Smear::ParticleMCS(mc.Get4Vector(),
                                                 mc.Id(), mc.GetStatus());
  p->SetStatus(mc.GetStatus());
  return p;
}

}  // anonymous namespace

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
  for (unsigned j(0); j < mMcEvent->GetNTracks(); j++) {
    const erhic::VirtualParticle* ptr = mMcEvent->GetTrack(j);
    if (!ptr) {
      continue;
    }  // if
    // If this is the scattered lepton, record the index.
    // Set the index even if the particle turns out to be outside the
    // acceptance (in which case it will just point to a NULL anyway).
    if (mMcEvent->ScatteredLepton() == ptr) {
      ParticleMCS* p = mDetector.Smear(*ptr);
      if (p) {
        p->SetStatus(ptr->GetStatus());
        event->SetScattered(j);
      }  // if
      event->AddLast(p);
      // Only set the index if the scattered electron is detected
    } else if (mMcEvent->BeamLepton() == ptr ||
            mMcEvent->BeamHadron() == ptr) {
      // It's convenient to keep the initial beams, unsmeared, in the
      // smeared event record, so copy their properties exactly
      event->AddLast(mcToSmear(*ptr));
    } else {
      ParticleMCS* p = mDetector.Smear(*ptr);
      if (p) {
        p->SetStatus(ptr->GetStatus());
      }  // if
      event->AddLast(p);
    }  // if
  }  // for
  // Fill the event-wise kinematic variables.
  mDetector.FillEventKinematics(event);
  return event;
}

}  // namespace Smear
