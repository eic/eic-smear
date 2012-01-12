#ifndef _ERHIC_BUILDTREE_EVENTRAPGAP_
#define _ERHIC_BUILDTREE_EVENTRAPGAP_

#include <string>

#include <Rtypes.h>

#include "EventBase.h"

/**
 Describes an event from the generator RAPGAP
 */
struct EventRapgap : public EventBase {
   
   virtual bool Parse(const std::string& );
   
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
   Double32_t sHat; // Partonic centre-of-mass energy

   ClassDef(EventRapgap, 1 )
};

#endif
