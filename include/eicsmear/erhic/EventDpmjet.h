/**
 \file
 Declaration of class erhic::EventDpmjet.

 \author    Thomas Burton
 \date      2011-08-19
 \copyright 2011 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_ERHIC_EVENTDPMJET_H_
#define INCLUDE_EICSMEAR_ERHIC_EVENTDPMJET_H_

#include <string>

#include <Rtypes.h>

#include "eicsmear/erhic/EventMC.h"

namespace erhic {

/**
 Describles the event from DPMJET.
 */
class EventDpmjet : public EventMC {
 public:
  /**
   Parses the event information from a text string.
   
   The string must have the following format (no newlines):
   \verbatim
   "0 eventnumber subprocessId hardProcessId particleCombination y Q2 x W2 nu
   theta photonFlux targetParton projectileParton xTargetParton
   xProjectileParton nucleon numTracks"
   \endverbatim
   Returns true in the event of a successful read operation,
   false in case of an error.
   */
  virtual bool Parse(const std::string&);

  Int_t ievent;
  Int_t I;
  Int_t process1;
  Int_t process2;
  Int_t IP;
  Int_t tgtparton;
  Int_t prjparton;
  Int_t nucleon;
  Double32_t xtgtparton;
  Double32_t xprjparton;
  Double32_t dtrueW2;
  Double32_t dtrueNu;
  Double32_t dtrueQ2;
  Double32_t dtrueY;
  Double32_t dtrueX;
  Double32_t theta_Evt;
  Double32_t photonFlux;

  ClassDef(erhic::EventDpmjet, 1)
};

}  // namespace erhic

#endif  // INCLUDE_EICSMEAR_ERHIC_EVENTDPMJET_H_
