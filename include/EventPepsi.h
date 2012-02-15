#ifndef _ERHIC_BUILDTREE_EVENTPEPSI_
#define _ERHIC_BUILDTREE_EVENTPEPSI_

#include <string>

#include <Rtypes.h>

#include "EventBase.h"

/**
 Describes an event from the generator PEPSI
 \todo Add accessor and setter methods
 \todo Add comments describing each field
 */
struct EventPepsi : public EventBase {
   
   virtual bool Parse(const std::string& );
   
   Int_t nucleon;
   Int_t struckparton;
   Int_t partontrck;
   Int_t genevent;
   Int_t subprocess;
   Double32_t trueY;
   Double32_t trueQ2;
   Double32_t trueX;
   Double32_t trueW2;
   Double32_t trueNu;
   Double32_t FixedWeight;
   Double32_t Weight;
   Double32_t dxsec;
   Double32_t ExtraWeight;
   Double32_t Dilute;
   Double32_t F1;
   Double32_t F2;
   Double32_t A1;
   Double32_t A2;
   Double32_t R;
   Double32_t DePol;
   Double32_t D;
   Double32_t Eta;
   Double32_t Eps;
   Double32_t Chi;
   Double32_t gendilut;
   Double32_t genF1;
   Double32_t genF2;
   Double32_t genA1;
   Double32_t genA2;
   Double32_t genR;
   Double32_t genDepol;
   Double32_t gend;
   Double32_t geneta;
   Double32_t geneps;
   Double32_t genchi;
   Double32_t SigCorr;
   Double32_t radgamEnucl;
   
   ClassDef(EventPepsi, 1 )
};

#endif