//
// Event.cxx
//
// Created by TB on 10/10/11.
// Copyright 2011 BNL. All rights reserved.
//


/**\todo Tidy up*/
#include "EventSmear.h"

//#ifdef USE_NAMESPACE_ERHIC   
namespace Smear {
   //#endif   
   
   void Event::Reset() {
      for( unsigned i(0); i < particles.size(); ++i ) {
         if( GetTrack(i) ) {
            delete GetTrack(i);
         } // if
      } // for
      particles.clear();
   }
   
/*   void Event::Grow(unsigned n) {
      for(unsigned i(0); i < n; ++i) {
         particles.push_back(new ::Smear::ParticleMCS);
         particles.back()->id = i;
      }
   }
   */
   void Event::AddLast(TrackType* track) {
      particles.push_back(track);
   }
   
   //#ifdef USE_NAMESPACE_ERHIC   
} // namespace erhic
  //#endif
