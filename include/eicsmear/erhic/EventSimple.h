/**
 \file      Declaration of class erhic::EventSimple. 
 \author    Barak Schmookler
 \date      2020-03-19
 */

#ifndef INCLUDE_EICSMEAR_ERHIC_EVENTSIMPLE_H_
#define INCLUDE_EICSMEAR_ERHIC_EVENTSIMPLE_H_

#include <string>
#include <Rtypes.h>
#include "eicsmear/erhic/EventMC.h"

namespace erhic {

/**
 Describes an event from the generator SIMPLE.
 */
class EventSimple : public EventMC {
 public:
  /**
   Constructor.
   */
  EventSimple();

  /**
   Parses the event information from a text string.
   
   The string must have the following format (no newlines):
   \verbatim
   "0 eventnumber numParticles
   \endverbatim
   Returns true in the event of a successful read operation,
   false in case of an error.
   */
  virtual bool Parse(const std::string&);
  const ParticleMC* ScatteredLepton() const;

  Double32_t numParticles;

  ClassDef(erhic::EventSimple, 1)
};

}  // namespace erhic

#endif  // INCLUDE_EICSMEAR_ERHIC_EVENTSIMPLE_H_
