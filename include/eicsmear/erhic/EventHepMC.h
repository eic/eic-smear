/**
 \file
 Declaration of class erhic::EventMilou.
 
 \author    Thomas Burton
 \date      2011-07-07
 \copyright 2011 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_ERHIC_EVENTHEPMC_H_
#define INCLUDE_EICSMEAR_ERHIC_EVENTHEPMC_H_

#include <string>

#include <Rtypes.h>

#include "eicsmear/erhic/EventMC.h"

namespace erhic {

/**
 Describes an event from the generator HEPMC.
 */
class EventHepMC : public EventMC {
 public:
  /**
   Constructor.
   */
  EventHepMC();

  /**
   dummy - the reading is done by a hepmc reader
   */
  virtual bool Parse(const std::string&) {return true;}

  /* Bool_t radcorr; */
  /* Double32_t weight; */
  /* Double32_t trueX; */
  /* Double32_t trueQ2; */
  /* Double32_t trueY; */
  /* Double32_t trueT; */
  /* Double32_t truePhi; */
  /* Double32_t phibelgen;  ///> the azimuthal angle between the production */
  /*                        ///> and the scattering plane */
  /* Double32_t phibelres;  ///> the resolution of the previous angle */
  /*                        ///> according to H1 */
  /* Double32_t phibelrec;  ///> the reconstructed phi angle */

  ClassDef(erhic::EventHepMC, 1)
};

}  // namespace erhic

#endif  // INCLUDE_EICSMEAR_ERHIC_EVENTHEPMC_H_
