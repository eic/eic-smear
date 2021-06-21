/**
 \file Implementation of class erhic::EventDEMP. 
 \author    Wenliang Li
 \date      2021-06-21
 \email     wenliang.billlee@gmail.com
 */


#ifndef INCLUDE_EICSMEAR_ERHIC_EVENTDEMP_H_
#define INCLUDE_EICSMEAR_ERHIC_EVENTDEMP_H_

#include <string>
#include <Rtypes.h>
#include "eicsmear/erhic/EventMC.h"

namespace erhic {

/**
 Describes an event from the generator DEMP.
 */
class EventDEMP : public EventMC {
 public:
  /**
   Constructor.
   */
  EventDEMP();

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
  Double32_t weight;

  ClassDef(erhic::EventDEMP, 1)
};

}  // namespace erhic

#endif  // INCLUDE_EICSMEAR_ERHIC_EVENTDEMP_H_
