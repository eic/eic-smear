/**
 Event.h
 
 \file
 Declaration of class Event.
 
 \author Thomas Burton 
 \date 5/2/12
 \copyright 2012 BNL. All rights reserved.
 */

#ifndef _EICSMEAR_HADRONIC_EventMC_H_
#define _EICSMEAR_HADRONIC_EventMC_H_

#include <vector>

#include <TClonesArray.h>

#include <eicsmear/hadronic/ParticleMC.h>
#include "eicsmear/erhic/VirtualEvent.h"

namespace erhic {

namespace hadronic {

class ParticleMC;

/**
 A realisation of erhic::VirtualEvent for a hadron-hadron
 Monte Carlo event.
 */
class EventMC : public erhic::VirtualEvent {
public:

   /** Destructor */
   virtual ~EventMC();

   /** Constructor */
   EventMC();

   /**
    Returns the nth track from the event.
    Indices run from 0 to (n-1).
    */
   virtual const ParticleMC* GetTrack(UInt_t) const;

   /**
    Returns the nth track from the event.
    Indices run from 0 to (n-1).
    */
   virtual ParticleMC* GetTrack(UInt_t);

   /**
    Returns the number of tracks in the event.
    */
   virtual UInt_t GetNTracks() const;

   /**
    Add a track to the list.
    The track must be dynamically allocated and is owned by the event.
    Returns the new list size.
   */
   virtual UInt_t Add(ParticleMC*);

   /**
    Returns the centre-of-mass energy of the event (GeV)
    */
   virtual Double_t GetCentreOfMassEnergy() const;

   virtual void Clear(Option_t* = "");

protected:

   TClonesArray mTracks;

   ClassDef(erhic::hadronic::EventMC, 1)
};

} // namespace hadronic

} // namespace erhic

#endif
