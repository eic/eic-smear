#ifndef _EICSMEAR_EVENTDJANGOH_
#define _EICSMEAR_EVENTDJANGOH_

#include <string>

#include <Rtypes.h>

#include "eicsmear/erhic/EventMC.h"

namespace erhic {
   
   /**
    Describes an event from the generator DJANGOH.
    */
   class EventDjangoh : public EventMC {
   public:
      EventDjangoh() { }
      virtual bool Parse(const std::string& );
      
      Int_t nucleon;
      Int_t IChannel;
      Int_t dprocess;
      Int_t dstruckparton;
      Int_t dpartontrck;
      Double32_t dY;
      Double32_t dQ2;
      Double32_t dX;
      Double32_t dW2;
      Double32_t dNu;
      Double32_t dtrueY;
      Double32_t dtrueQ2;
      Double32_t dtrueX;
      Double32_t dtrueW2;
      Double32_t dtrueNu;   
      Double32_t sigTot;
      Double32_t sigTotErr;
      Double32_t D;
      Double32_t F1NC;
      Double32_t F3NC;
      Double32_t G1NC;
      Double32_t G3NC;
      Double32_t A1NC;
      Double32_t F1CC;
      Double32_t F3CC;
      Double32_t G1CC;
      Double32_t G5CC;
      
      ClassDef(erhic::EventDjangoh, 1 )
   };
} // namespace erhic

#endif
