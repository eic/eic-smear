/**
 \file
 Declaration of class erhic::hadronic::Event.
 
 \author    Thomas Burton
 \date      2012-05-02
 \copyright 2012 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_HADRONIC_EVENTMC_H_
#define INCLUDE_EICSMEAR_HADRONIC_EVENTMC_H_

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
  /**
   Destructor.
   */
  virtual ~EventMC();

  /**
   Constructor.
   */
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
  virtual UInt_t Add(ParticleMC* particle);

  /**
   Returns the centre-of-mass energy of the event (GeV)
   */
  virtual Double_t GetCentreOfMassEnergy() const;

  /**
   Clear the track list.
   */
  virtual void Clear(Option_t* = "");

 protected:
  TClonesArray mTracks;

  ClassDef(erhic::hadronic::EventMC, 1)
};

}  // namespace hadronic
}  // namespace erhic

#endif  // INCLUDE_EICSMEAR_HADRONIC_EVENTMC_H_
