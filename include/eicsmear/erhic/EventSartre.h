/**
 \file      Declaration of class erhic::EventSartre. 
 \author    Barak Schmookler
 \date      2020-04-20
 */

#ifndef INCLUDE_EICSMEAR_ERHIC_EVENTSARTRE_H_
#define INCLUDE_EICSMEAR_ERHIC_EVENTSARTRE_H_

#include <string>
#include <Rtypes.h>
#include "eicsmear/erhic/EventMC.h"

namespace erhic {

/**
 Describes an event from the generator SIMPLE.
 */
class EventSartre : public EventMC {
 public:
  /**
   Constructor.
   */
  EventSartre();

  /**
   Parses the event information from a text string.
   
   The string must have the following format (no newlines):
   \verbatim
   "0 ievent genevent truet trueQ2 truex truey trueW2
    truenu truexpom s_cm pol dmode bup"
   \endverbatim
   Returns true in the event of a successful read operation,
   false in case of an error.
   */
  virtual bool Parse(const std::string&);

  Int_t       genevent;         ///< Trials required for this event (dummy right now)
  Double32_t  trueT; 
  Double32_t  trueQ2;
  Double32_t  trueX;
  Double32_t  trueY;
  Double32_t  trueW2;
  Double32_t  trueNu;
  Double32_t  trueXpom;
  Double32_t  s_cm;
  Int_t       pol;
  Int_t       dmode;
  Int_t       bup; 

  ClassDef(erhic::EventSartre, 1)
};

}  // namespace erhic

#endif  // INCLUDE_EICSMEAR_ERHIC_EVENTSARTRE_H_
