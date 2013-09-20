/**
 \file
 Implementation of class Smear::EventSmear.
 
 \author    Michael Savastio
 \date      2011-10-10
 \copyright 2011 Brookhaven National Lab
 */

#include "eicsmear/smear/EventSmear.h"

#include <iostream>
#include <vector>

namespace Smear {

Event::Event()
: nTracks(0)
, mScatteredIndex(-1) {
}

Event::~Event() {
  ClearParticles();
}

void Event::ClearParticles() {
  for (unsigned i(0); i < particles.size(); ++i) {
    if (GetTrack(i)) {
      delete GetTrack(i);
    }  // if
  }  // for
}

void Event::Reset() {
  ClearParticles();
  *this = Event();
}

void Event::AddLast(ParticleMCS* track) {
  particles.push_back(track);
}

// The scattered lepton should be the first non-NULL entry in the track list
const ParticleMCS* Event::ScatteredLepton() const {
  if (mScatteredIndex > -1 &&
      mScatteredIndex < static_cast<int>(GetNTracks())) {
    return GetTrack(mScatteredIndex);
  }  // if
  return NULL;
}

// Get the particles that belong to the hadronic final state.
// The stored Particle* are pointers to the original particles in the event
// so don't delete them!
void Event::HadronicFinalState(ParticlePtrList& final) const {
  // Skip the first two entries, as these are the incident beams
  for (unsigned i(2); i < GetNTracks(); ++i) {
    if (!GetTrack(i)) {
      continue;
    }  // if
    if (GetTrack(i) != ScatteredLepton()) {
      final.push_back(GetTrack(i));
    }  // if
  }  // for
}

std::vector<const erhic::VirtualParticle*> Event::GetTracks() const {
  std::vector<const erhic::VirtualParticle*> tracks;
  for (unsigned i(0); i < GetNTracks(); ++i) {
    tracks.push_back(GetTrack(i));
  }  // for
  return tracks;
}

void Event::SetScattered(int index) {
  if (index >= 0) {
    mScatteredIndex = index;
  }  // if
}
void Event::Print(Option_t* /* unused */) const {
  std::cout <<
  "x:  " << GetX() << std::endl <<
  "Q2: " << GetQ2() << std::endl <<
  "y:  " << GetY() << std::endl;
  for (unsigned i(0); i < GetNTracks(); ++i) {
    if (GetTrack(i)) {
      GetTrack(i)->Print();
    }  // if
  }  // for
}

}  // namespace Smear
