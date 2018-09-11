/**
 \file
 Declaration of class erhic::hadronic::EventPythiaPP.
 
 \author    Thomas Burton
 \date      2012-05-02
 \copyright 2012 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_HADRONIC_EVENTPYTHIA_H_
#define INCLUDE_EICSMEAR_HADRONIC_EVENTPYTHIA_H_

#include <Rtypes.h>  // For ClassDef macro

#include <eicsmear/hadronic/EventMC.h>

namespace erhic {
namespace hadronic {

class EventPythiaPP : public EventMC {
 public:
  /**
   Destructor.
   */
  virtual ~EventPythiaPP();

  /**
   Default constructor.
   */
  EventPythiaPP();

  /**
   Initialise event kinematics.
   */
  EventPythiaPP(double Q2, double x1, double x2);

  /**
   Copy constructor.
   */
  EventPythiaPP(const EventPythiaPP&);

  /**
   Assignment operator.
   */
  EventPythiaPP& operator=(const EventPythiaPP&);

  /**
   Q-squared of the interaction.
   */
  virtual Double_t GetQ2() const;

  /**
   x of the parton from the +z beam.
   */
  virtual Double_t GetX1() const;

  /**
   x of the parton from the -z beam.
   */
  virtual Double_t GetX2() const;

 protected:
  Double32_t QSquared;
  Double32_t x1;
  Double32_t x2;

  ClassDef(erhic::hadronic::EventPythiaPP, 1)
};

}  // namespace hadronic
}  // namespace erhic

#endif  // INCLUDE_EICSMEAR_HADRONIC_EVENTPYTHIA_H_
