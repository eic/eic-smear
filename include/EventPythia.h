#ifndef _ERHIC_BUILDTREE_EVENTPYTHIA_
#define _ERHIC_BUILDTREE_EVENTPYTHIA_

#include <string>

#include <Rtypes.h>

#include "EventBase.h"

/**
 Describes an event from the generator PYTHIA
 */
struct EventPythia : public EventBase {
   
   EventPythia(const std::string& = "" );
   
   /**
    Parses the event information from a text string.
    Returns true in the event of a successful read operation,
    false in case of an error.
    */
   virtual bool Parse(const std::string& );
   
	Int_t       nucleon;          ///< PDG code of the hadron beam
   Int_t       tgtparton;        ///< PDG code of the struck parton in the hadron beam
   Int_t       beamparton;       ///< Parton interacting with the hadron beam in the case of resolved photon processes and soft VMD
   Int_t       genevent;
	Double32_t  xtgtparton;
   Double32_t  xbeamparton;
	Double32_t  thetabeamparton;
   Double32_t  leptonphi;
   Double32_t  F1;
   Double32_t  sigma_rad;
	Double32_t  t_hat;
   Double32_t  u_hat;
   Double32_t  Q2_hat;
   Double32_t  SigRadCor;
   Double32_t  EBrems;
   Double32_t  photonflux;
	Double32_t  trueY;
   Double32_t  trueQ2;
   Double32_t  trueX;
   Double32_t  trueW2;
   Double32_t  trueNu;
   Double32_t  F2;
   Double32_t  R;
   Double32_t  pt2_hat;
   Double32_t  sHat;             ///< Partonic centre-of-mass energy
   
   virtual void SetNucleon(int n) { nucleon = n; }
   virtual void SetTargetParton(int n) { tgtparton = n; }
   virtual void SetBeamParton(int n) { beamparton = n; }
   virtual void SetGenEvent(int n) { genevent = n; }
   virtual void SetTargetPartonX(double xB) { xtgtparton = xB; }
   virtual void SetBeamPartonX(double xB) { xbeamparton = xB; }
   virtual void SetBeamPartonTheta(double radians) { thetabeamparton = radians; }
   virtual void SetLeptonPhi(double radians) { leptonphi = radians; }
   virtual void SetF1(double f1) { F1 = f1; }
   virtual void SetF2(double f2) { F2 = f2; }
   virtual void SetSigmaRad(double sr) { sigma_rad = sr; }
   virtual void SetHardS(double s) { sHat = s; }
   virtual void SetHardT(double t) { t_hat = t; }
   virtual void SetHardU(double u) { u_hat = u; }
   virtual void SetHardQ2(double Q2) { Q2_hat = Q2; }
   virtual void SetHardPt2(double pt2) { pt2_hat = pt2; }
   virtual void SetSigRadCor(double s) { SigRadCor = s; }
   virtual void SetEBrems(double e) { EBrems = e; }
   virtual void SetPhotonFlux(double f) { photonflux = f; }
   virtual void SetTrueY(double inelasticity) { trueY = inelasticity; }
   virtual void SetTrueQ2(double Q2) { trueQ2 = Q2; }
   virtual void SetTrueX(double xB) { trueX = xB; }
   virtual void SetTrueW2(double W2) { trueW2 = W2; }
   virtual void SetTrueNu(double Nu) { trueNu = Nu; }
   virtual void SetR(double r) { R = r; }
   
   ClassDef(EventPythia, 1 )
};

#endif