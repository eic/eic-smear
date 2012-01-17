#ifndef _ERHIC_EVENT_H_
#define _ERHIC_EVENT_H_

// erhic_event.h
// BuildTree
//
// Created by TB on 8/19/11.
// Copyright 2011 BNL. All rights reserved.

#include <TObject.h>

namespace erhic {
   
   /**
    Abstract base class for a lepton-hadron event.
    An event is a collection of tracks plus (possibly) other quantities.
    The template argument should be the name of a "track" type,
    whatever that may be for your needs.
    The event stores a list of these track (in principle the
    track can be whatever you need a list of).
    */
   template<typename T>
   class VirtualEvent : public TObject {
      
   public:
      
      typedef T TrackType; // The type of track in the event
      
      virtual ~VirtualEvent() { }
      
      /**
       Returns Bjorken-x of the event.
       x<sub>B</sub> = Q<sup>2</sup>/(2p.q)
       */
      virtual Double_t GetX() const = 0;
      
      /**
       Returns the four-momentum transfer (exchange boson mass) Q<sup>2</sup>.
       Q<sup>2</sup> = 2EE`(1+cos(theta)) = (e-e`)<sup>2</sup>
       */
      virtual Double_t GetQ2() const = 0;
      
      /**
       Returns the event inelasticity.
       y = (p.q)/(p.e)
       */
      virtual Double_t GetY() const = 0;
      
      /**
       Returns the invariant mass of the hadronic final state.
       W<sup>2</sup> = M<sup>2</sup> + Q<sup>2</sup>(1-x)/x
       */
      virtual Double_t GetW2() const = 0;
      
      /**
       Returns the exchange boson energy in the beam hadron rest frame.
       nu = q.p/M
       */
      virtual Double_t GetNu() const = 0;
      
      /**
       Returns the nth track from the event.
       Indices run from 0 to (n-1).
       */
      virtual const TrackType* GetTrack(UInt_t) const = 0;
      
      /**
       Returns the nth track from the event.
       Indices run from 0 to (n-1).
       */
      virtual TrackType* GetTrack(UInt_t) = 0;
      
      /**
       Returns the number of tracks in the event.
       */
      virtual UInt_t GetNTracks() const = 0;
      
      /**
       Add a new track to the end of the track list.
       */
      virtual void AddLast(TrackType*) = 0;
      
      ClassDef(VirtualEvent, 1)
   };
   
} // namespace erhic

#endif
