/**
 VirtualEvent.h
 
 \file
 Declaration of class VirtualEvent.
 
 \author TB
 \date 8/19/11
 \copyright 2011 BNL. All rights reserved.
 */

#ifndef _ERHIC_EVENT_H_
#define _ERHIC_EVENT_H_

#include <vector>

#include <TObject.h>

namespace erhic {

   class VirtualParticle;

   /**
    Abstract base class for a lepton-hadron event.
    An "event" is defined here as a collection of tracks.
    */
   class VirtualEvent : public TObject {
      
   public:

      virtual ~VirtualEvent() { }

      /**
       Returns the nth track from the event.
       Indices run from 0 to (n-1).
       */
      virtual const VirtualParticle* GetTrack(UInt_t) const = 0;
      
      /**
       Returns the nth track from the event.
       Indices run from 0 to (n-1).
       */
      virtual VirtualParticle* GetTrack(UInt_t) = 0;
      
      /**
       Returns the number of tracks in the event.
       */
      virtual UInt_t GetNTracks() const = 0;

      typedef std::vector<const erhic::VirtualParticle*> ParticlePtrList;
      virtual void HadronicFinalState(ParticlePtrList&) const { }

      ClassDef(VirtualEvent, 1)
   };
} // namespace erhic

#endif
