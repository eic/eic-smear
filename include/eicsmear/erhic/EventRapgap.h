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


  /**  Returns a pointer to the scattered lepton in the event record.
       This is the first (only?) particle that matches the following:
       1) pdg code equals that of incident lepton beam.
       2) status code is 1 i.e. it's a stable/final-state particle.
       3) the parent is track 1 or 2
  */
  const ParticleMC* ScatteredLepton() const;

  /**
   Returns a pointer to the exchanged boson.   
   It would probably be the third track, but we'll go with the first status=21 boson
   that has particle 1 or 2 as parent
   */
  virtual const ParticleMC* ExchangeBoson() const;
  
  Double32_t Get_cs() const {return cs;}
  Double32_t Get_sigma_cs() const {return sigma_cs;}
  Double32_t Get_s() const {return s;}
  Double32_t Get_q2() const {return q2;}
  Double32_t Get_xgam() const {return xgam;}
  Double32_t Get_xpr() const {return xpr;}
  Double32_t Get_Pt_h() const {return Pt_h;}
  Double32_t Get_t() const {return t;}
  Double32_t Get_x_pom() const {return x_pom;}
  Double32_t Get_sHat2() const {return sHat2;}
  Double32_t Get_z() const {return z;}
  Double32_t Get_x1() const {return x1;}
  Double32_t Get_phi1() const {return phi1;}
  Double32_t Get_pt2_hat() const {return pt2_hat;}
  Double32_t Get_sHat() const {return sHat;}

 protected:

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
      
  ClassDef(erhic::EventRapgap, 2)

};

}  // namespace erhic

#endif  // INCLUDE_EICSMEAR_ERHIC_EVENTRAPGAP_H_
