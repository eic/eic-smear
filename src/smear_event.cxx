//
// Event.cxx
//
// Created by TB on 10/10/11.
// Copyright 2011 BNL. All rights reserved.
//

#include "EventSmear.h"

//#ifdef USE_NAMESPACE_ERHIC   
namespace Smear {
   //#endif   
   
   void Event::Reset() {
      for( unsigned i(0); i < smearedParticles.size(); ++i ) {
         if( GetTrack(i) ) {
            delete GetTrack(i);
         } // if
      } // for
      smearedParticles.clear();
   }
   
   void Event::Grow(unsigned n) {
      for(unsigned i(0); i < n; ++i) {
         smearedParticles.push_back(new ::Smear::ParticleMCS);
         smearedParticles.back()->id = i;
      }
   }
   
   void Event::AddLast(TrackType* track) {
      smearedParticles.push_back(track);
   }
   
   //#ifdef USE_NAMESPACE_ERHIC   
} // namespace erhic
  //#endif
