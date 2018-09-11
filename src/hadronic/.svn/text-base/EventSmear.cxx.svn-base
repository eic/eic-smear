/**
 \file Implementation of smearing claases for hadronic events.
 
 \author    Thomas Burton
 \date      2012-05-03
 \copyright 2012 Brookhaven National Lab
 */

#include <eicsmear/hadronic/EventSmear.h>

namespace erhic {
namespace hadronic {

EventSmear::~EventSmear() {
  for (unsigned i(0); i < GetNTracks(); ++i) {
    delete particles.at(i);
    particles.at(i) = NULL;
  }  // for
}

EventSmear::EventSmear() {
}

const Smear::ParticleMCS* EventSmear::GetTrack(UInt_t i) const {
  return particles.at(i);
}

Smear::ParticleMCS* EventSmear::GetTrack(UInt_t i) {
  return particles.at(i);
}

UInt_t EventSmear::GetNTracks() const {
  return particles.size();
}

void EventSmear::AddLast(Smear::ParticleMCS* p) {
  particles.push_back(p);
}

}  // namespace hadronic
}  // namespace erhic
