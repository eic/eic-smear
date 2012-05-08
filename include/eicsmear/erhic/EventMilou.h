#ifndef _EICSMEAR_EVENTMILOU_
#define _EICSMEAR_EVENTMILOU_

#include <string>

#include <Rtypes.h>

#include "eicsmear/erhic/EventMC.h"

namespace erhic {
   
   /**
    Describes an event from the generator MILOU
    */
   class EventMilou : public EventMC {
   public:
      EventMilou();
      
      virtual bool Parse(const std::string&);
      
      Double_t GetPhiBelGen() const;
      Double_t GetPhiBelRes() const;
      Double_t GetPhiBelRec() const;
      
      Bool_t radcorr;
      Double32_t weight;
      Double32_t trueX;
      Double32_t trueQ2;
      Double32_t trueY;
      Double32_t trueT;
      Double32_t truePhi;
      Double32_t phibelgen; ///> the azimuthal angle between the production and the scattering plane
      Double32_t phibelres; ///> the resolution of the previous angle according to H1
      Double32_t phibelrec; ///> the reconstructed phi angle
      
      ClassDef(erhic::EventMilou, 1 )
   };
   
   inline Double_t EventMilou::GetPhiBelGen() const {
      return phibelgen;
   }
   
   inline Double_t EventMilou::GetPhiBelRes() const {
      return phibelres;
   }
   
   inline Double_t EventMilou::GetPhiBelRec() const {
      return phibelrec;
   }
} // namespace erhic

#endif
