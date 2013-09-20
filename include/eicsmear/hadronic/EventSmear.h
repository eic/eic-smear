/**
 \file Declaration of smearing claases for hadronic events.
 
 \author    Thomas Burton
 \date      2012-05-03
 \copyright 2012 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_HADRONIC_EVENTSMEAR_H_
#define INCLUDE_EICSMEAR_HADRONIC_EVENTSMEAR_H_

#include <vector>

#include <Rtypes.h>  // For ClassDef macro
#include <TBranch.h>

#include "eicsmear/erhic/VirtualEvent.h"
#include "eicsmear/hadronic/EventMC.h"
#include "eicsmear/smear/EventFactory.h"
#include "eicsmear/smear/ParticleMCS.h"

namespace erhic {
namespace hadronic {

/**
 Realisation of hadronic::EventMC as an event with detector smearing.
 */
class EventSmear : public VirtualEvent {
 public:
  /**
   Destructor.
   */
  virtual ~EventSmear();

  /**
   Constructor.
   */
  EventSmear();

  /**
   Returns the numbered track for the event.
   */
  virtual const Smear::ParticleMCS* GetTrack(UInt_t) const;

  /**
   \overload
   */
  virtual Smear::ParticleMCS* GetTrack(UInt_t);

  /**
   Returns the number of tracks in the event.
   */
  virtual UInt_t GetNTracks() const;

  /**
   Add a particle to the end of the list.
   */
  virtual void AddLast(Smear::ParticleMCS*);

 protected:
  std::vector<Smear::ParticleMCS*> particles;

  ClassDef(erhic::hadronic::EventSmear, 1)
};

}  // namespace hadronic
}  // namespace erhic

namespace Smear {

/**
 Factory class for smeared hadronic events.
 */
class HadronicEventBuilder : public EventFactory<erhic::hadronic::EventSmear> {
 public:
  /**
   Destructor.
   */
  virtual ~HadronicEventBuilder() { }

  /**
   Constructor.
   Initialise with the Detector performing particle smearing and
   the tree branch providing DIS Monte Carlo events.
   */
  HadronicEventBuilder(const Detector& d, TBranch& mcBranch)
  : mDetector(d)
  , mMcEvent(NULL) {
    mcBranch.SetAddress(&mMcEvent);
  }

  /**
   Create a smeared event corresponding to the current DIS Monte Carlo
   event in the input branch passed to the constructor.
   The user should call TTree::GetEntry() themselves between calls to
   Create().
   */
  virtual EventType* Create() {
    EventType* event = new EventType;
    for (unsigned j(0); j < mMcEvent->GetNTracks(); j++) {
      const erhic::VirtualParticle* ptr = mMcEvent->GetTrack(j);
      if (ptr) {
        event->AddLast(mDetector.Smear(*ptr));
      }  // if
    }  // for
    return event;
  }

 protected:
  Detector mDetector;
  erhic::hadronic::EventMC* mMcEvent;
};

}  // namespace Smear

#endif  // INCLUDE_EICSMEAR_HADRONIC_EVENTSMEAR_H_
