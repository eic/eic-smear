/**
 \file
 Declaration of class erhic::EventRapgap.

 \author    Thomas Burton
 \date      2011-07-07
 \copyright 2011 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_ERHIC_EVENTRAPGAP_H_
#define INCLUDE_EICSMEAR_ERHIC_EVENTRAPGAP_H_

#include <string>

#include <Rtypes.h>

#include "eicsmear/erhic/EventMC.h"

namespace erhic {

/**
 Describes an event from the generator RAPGAP.
 */
class EventRapgap : public EventMC {
 public:
  virtual bool Parse(const std::string&);

  Int_t idir;
  Int_t idisdif;
  Int_t genevent;

  Double32_t cs;
  Double32_t sigma_cs;
  Double32_t s;
  Double32_t q2;
  Double32_t xgam;
  Double32_t xpr;
  Double32_t Pt_h;
  Double32_t t;
  Double32_t x_pom;
  Double32_t sHat2;
  Double32_t z;
  Double32_t x1;
  Double32_t phi1;
  Double32_t pt2_hat;
  Double32_t sHat;  // Partonic centre-of-mass energy

  ClassDef(erhic::EventRapgap, 1)
};

}  // namespace erhic

#endif  // INCLUDE_EICSMEAR_ERHIC_EVENTRAPGAP_H_
