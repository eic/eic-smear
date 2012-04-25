#include "EventSmear.h"

namespace Smear {

   Event::Event()
   : nTracks(0)
   , x(0.)
   , QSquared(0.)
   , y(0.)
   , WSquared(0.)
   , nu(0.)
   , yJB(0.)
   , QSquaredJB(0.)
   , xJB(0.)
   , WSquaredJB(0.)
   , yDA(0.)
   , QSquaredDA(0.)
   , xDA(0.)
   , WSquaredDA(0.) {
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

   void Event::AddLast(TrackType* track) {
      particles.push_back(track);
   }

} // namespace erhic
