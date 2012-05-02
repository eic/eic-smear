#include "EventSmear.h"

namespace Smear {

   Event::Event()
   : nTracks(0)
   , mScatteredIndex(-1) {
   }

   Event::~Event() {
      ClearParticles();
   }

   void Event::ClearParticles() {
      for(unsigned i(0); i < particles.size(); ++i) {
         if(GetTrack(i)) {
            delete GetTrack(i);
         } // if
      } // for
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
      if(mScatteredIndex > -1 and mScatteredIndex < int(GetNTracks())) {
         return GetTrack(mScatteredIndex);
      } // if
      return NULL;
   }

   // Get the particles that belong to the hadronic final state.
   // The stored Particle* are pointers to the original particles in the event
   // so don't delete them!
   void Event::HadronicFinalState(ParticlePtrList& final) const {
      // Skip the first two entries, as these are the incident beams
      for(unsigned i(2); i < GetNTracks(); ++i) {
         if(not GetTrack(i)) {
            continue;
         } // if
         if(GetTrack(i) not_eq ScatteredLepton()) {
            final.push_back(GetTrack(i));
         } // if
      } // for
   }
} // namespace erhic
