/**
 EventPythia.h
 
 \file
 Declaration of class EventPythia.
 
 \author Thomas Burton
 \date 8/31/11
 \copyright 2011 BNL. All rights reserved.
 */

#ifndef _ERHIC_BUILDTREE_EVENTPYTHIA_
#define _ERHIC_BUILDTREE_EVENTPYTHIA_

#include <string>

#include <Rtypes.h>

#include "EventMC.h"

namespace erhic {
   
   /**
    Describes an event from the EIC PYTHIA6 implementation.
    \todo Add accessor methods
    \todo Complete method comments
    */
   class EventPythia : public EventMC {
      
   public:
      
      /**
       Constructor.
       @param [in] str A text string setting event-wise quantities. The string
       format should be (no newlines):
       "I ievent genevent subprocess nucleon targetparton xtargparton
       beamparton xbeamparton thetabeamprtn truey trueQ2 truex trueW2 trueNu
       leptonphi s_hat t_hat u_hat pt2_hat Q2_hat F2 F1 R sigma_rad SigRadCor
       EBrems photonflux nrTracks"
       */
      EventPythia(const std::string& str = "");
      
      /** Destructor */
      virtual ~EventPythia();
      
      /**
       Parses the event information from a text string.
       See the constructor for the string format.
       Returns true in the event of a successful read operation,
       false in case of an error.
       */
      virtual bool Parse(const std::string&);
      
      /**
       Sets the nucleon species.
       @param n [in] PDG code of the hadron beam, see MSTI(12)
       */
      virtual void SetNucleon(int n);
      
      /**
       Sets the target parton species.
       @param n [in] PDG code of the struck parton in the hadron beam,
       see MSTI(16)
       */
      virtual void SetTargetParton(int n);
      
      /**
       Sets the beam parton species.
       @param n [in] PDG code of the parton interacting with the hadron beam in
       the case of resolved photon processes and soft VMD, see MSTI(15)
       */
      virtual void SetBeamParton(int n);
      
      /**
       Sets the number of trials required to generate this event
       @param n [in] The number of trials
       */
      virtual void SetGenEvent(int n);
      virtual void SetTargetPartonX(double xB);
      virtual void SetBeamPartonX(double xB);
      virtual void SetBeamPartonTheta(double radians);
      virtual void SetLeptonPhi(double radians);
      virtual void SetF1(double f1);
      virtual void SetF2(double f2);
      virtual void SetSigmaRad(double sr);
      virtual void SetHardS(double s);
      virtual void SetHardT(double t);
      virtual void SetHardU(double u);
      virtual void SetHardQ2(double Q2);
      virtual void SetHardPt2(double pt2);
      virtual void SetSigRadCor(double s);
      virtual void SetEBrems(double e);
      virtual void SetPhotonFlux(double f);
      virtual void SetTrueY(double inelasticity);
      virtual void SetTrueQ2(double Q2);
      virtual void SetTrueX(double xB);
      virtual void SetTrueW2(double W2);
      virtual void SetTrueNu(double Nu);
      virtual void SetR(double r);
      
//   protected:
      
      // Inline comments after field names will appear in ROOT
      // when EventPythia::Dump() is called.
      
      Int_t       nucleon;          ///< PDG code of the hadron beam,
                                    ///< see MSTI(12)
      Int_t       tgtparton;        ///< PDG code of the struck parton
                                    ///< in the hadron beam, see MSTI(16)
      Int_t       beamparton;       ///< Parton interacting with the hadron
                                    ///< beam in the case of resolved photon
                                    ///< processes and soft VMD, see MSTI(15)
      Int_t       genevent;         ///< Trials required for this event
      Double32_t  xtgtparton;       ///< Momentum fraction taken by the
                                    ///< target parton, see PARI(34)
      Double32_t  xbeamparton;      ///< Momentum fraction taken by the
                                    ///< beam parton, see PARI(33)
      Double32_t  thetabeamparton;  ///< Polar angle of the beam parton in
                                    ///< the cm frame, between 0 and pi
                                    ///< radians, see PARI(53)
      Double32_t  leptonphi;        ///< Azimuthal angle of the scattered
                                    ///< lepton in the cm frame
      Double32_t  F1;               ///< Value used for radiative corrections
      Double32_t  sigma_rad;        ///< Value used for radiative corrections
      Double32_t  t_hat;            ///< Mandelstam t of the hard subprocess,
                                    ///< see PARI(15)
      Double32_t  u_hat;            ///< Mandelstam u of the hard subprocess,
                                    ///< see PARI(16)
      Double32_t  Q2_hat;           ///< Q<sup>2</sup> of the hard subprocess,
                                    ///< see PARI(22)
      Double32_t  SigRadCor;        ///< Value used for radiative corrections
      Double32_t  EBrems;           ///< Energy per radiative photon in the
                                    ///< nuclear rest frame.
      Double32_t  photonflux;       ///< Flux factor, see VINT(319)
      Double32_t  trueY;            ///< Generated y of the event,
                                    ///< see VINT(309)
      Double32_t  trueQ2;           ///< Generated Q<sup>2</sup> of the event,
                                    ///< see VINT(307)
      Double32_t  trueX;            ///< Generated x of the event
      Double32_t  trueW2;           ///< Generated W<sup>2</sup> of the event
      Double32_t  trueNu;           ///< Generated nu of the event
      Double32_t  F2;               ///< Value used for radiative corrections
      Double32_t  R;                ///< Value used for radiative corrections
      Double32_t  pt2_hat;          ///< Squared p<sub>T</sub> of the hard
                                    ///< subprocess, see PARI(18)
      Double32_t  sHat;             ///< Mandelstam s of the hard subprocess,
                                    ///< see PARI(14)
      
      ClassDef(erhic::EventPythia, 1)
   };
   
   inline void EventPythia::SetNucleon(int n) {
      nucleon = n;
   }
   inline void EventPythia::SetTargetParton(int n) {
      tgtparton = n;
   }
   inline void EventPythia::SetBeamParton(int n) {
      beamparton = n;
   }
   inline void EventPythia::SetGenEvent(int n) {
      genevent = n;
   }
   inline void EventPythia::SetTargetPartonX(double xB) {
      xtgtparton = xB;
   }
   inline void EventPythia::SetBeamPartonX(double xB) {
      xbeamparton = xB;
   }
   inline void EventPythia::SetBeamPartonTheta(double radians) {
      thetabeamparton = radians;
   }
   inline void EventPythia::SetLeptonPhi(double radians) {
      leptonphi = radians;
   }
   inline void EventPythia::SetF1(double f1) {
      F1 = f1;
   }
   inline void EventPythia::SetF2(double f2) {
      F2 = f2;
   }
   inline void EventPythia::SetSigmaRad(double sr) {
      sigma_rad = sr;
   }
   inline void EventPythia::SetHardS(double s) {
      sHat = s;
   }
   inline void EventPythia::SetHardT(double t) {
      t_hat = t;
   }
   inline void EventPythia::SetHardU(double u) {
      u_hat = u;
   }
   inline void EventPythia::SetHardQ2(double Q2) {
      Q2_hat = Q2;
   }
   inline void EventPythia::SetHardPt2(double pt2) {
      pt2_hat = pt2;
   }
   inline void EventPythia::SetSigRadCor(double s) {
      SigRadCor = s;
   }
   inline void EventPythia::SetEBrems(double e) {
      EBrems = e;
   }
   inline void EventPythia::SetPhotonFlux(double f) {
      photonflux = f;
   }
   inline void EventPythia::SetTrueY(double inelasticity) {
      trueY = inelasticity;
   }
   inline void EventPythia::SetTrueQ2(double Q2) {
      trueQ2 = Q2;
   }
   inline void EventPythia::SetTrueX(double xB) {
      trueX = xB;
   }
   inline void EventPythia::SetTrueW2(double W2) {
      trueW2 = W2;
   }
   inline void EventPythia::SetTrueNu(double Nu) {
      trueNu = Nu;
   }
   inline void EventPythia::SetR(double r) {
      R = r;
   }
   
} // namespace erhic

#endif
