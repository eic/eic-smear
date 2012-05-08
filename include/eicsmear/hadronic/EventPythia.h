/**
 EventPythiaPP.h
 
 \file
 Declaration of class EventPythiaPP.
 
 \author Thomas Burton 
 \date 5/2/12
 \copyright 2012 BNL. All rights reserved.
 */

#ifndef _EICSMEAR_HADRONIC_EventPythiaPP_H_
#define _EICSMEAR_HADRONIC_EventPythiaPP_H_

#include <Rtypes.h> // For ClassDef macro

#include <eicsmear/hadronic/EventMC.h>

namespace erhic {
   namespace hadronic {
      class EventPythiaPP : public EventMC {
      public:
         virtual ~EventPythiaPP();
         EventPythiaPP();
         EventPythiaPP(double Q2, double x1, double x2);
         EventPythiaPP(const EventPythiaPP&);
         EventPythiaPP& operator=(const EventPythiaPP&);
         /** Q-squared of the interaction */
         virtual Double_t GetQ2() const;
         /** x of the parton from the +z beam */
         virtual Double_t GetX1() const;
         /** x of the parton from the -z beam */
         virtual Double_t GetX2() const;
      protected:
         Double32_t QSquared;
         Double32_t x1;
         Double32_t x2;
      private:
         ClassDef(erhic::hadronic::EventPythiaPP, 1)
      };
   } // namespace hadronic
} // namespace erhic

#endif
