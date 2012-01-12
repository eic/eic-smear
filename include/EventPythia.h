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
   
   ClassDef(EventPythia, 1 )
};

#endif