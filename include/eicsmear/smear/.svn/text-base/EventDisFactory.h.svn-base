/**
 \file
 Declaration of class Smear::EventDisFactory.
 
 \author    Michael Savastio
 \date      2011-08-19
 \copyright 2011 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_SMEAR_EVENTDISFACTORY_H_
#define INCLUDE_EICSMEAR_SMEAR_EVENTDISFACTORY_H_

#include "eicsmear/smear/Detector.h"
#include "eicsmear/smear/EventSmear.h"
#include "eicsmear/smear/EventFactory.h"

class TBranch;

namespace erhic {

class EventDis;

}  // namespace erhic

namespace Smear {

/**
 Factory class for smeared DIS events.
 */
class EventDisFactory : public EventFactory<Smear::Event> {
 public:
  /**
   Destructor.
   */
  virtual ~EventDisFactory();

  /**
   Constructor.
   Initialise with the Detector performing particle smearing and
   the tree branch providing DIS Monte Carlo events.
   */
  EventDisFactory(const Detector&, TBranch&);

  /**
   Create a smeared event corresponding to the current DIS Monte Carlo
   event in the input branch passed to the constructor.
   The user should call TTree::GetEntry() on the tree with the
   input branch between calls to Create().
   */
  virtual Event* Create();

  erhic::VirtualEvent* GetEvBufferPtr();

 protected:
  Detector mDetector;
  erhic::EventDis* mMcEvent;
};

inline erhic::VirtualEvent* EventDisFactory::GetEvBufferPtr() {
  return mMcEvent;
}

}  // namespace Smear

#endif  // INCLUDE_EICSMEAR_SMEAR_EVENTDISFACTORY_H_
