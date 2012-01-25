//
// EventFactory.cxx
//
// Created by TB on 11/1/11.
// Copyright 2011 BNL. All rights reserved.
//

#include <memory>

#include "EventFactory.h"
#include "EventPythia.h"
#include "EventMilou.h"
#include "EventDjangoh.h"
#include "EventDpmjet.h"
#include "EventRapgap.h"
#include "EventPepsi.h"
#include "functions.h" // For getFirstNonBlank()
#include "ParticleMC.h"

namespace erhic {
   
   
   template<typename T>
   bool EventFromAsciiFactory<T>::AtEndOfEvent() const {
      return mLine.find("finished") not_eq std::string::npos;
   }
   
   
   template<typename T>
   T* EventFromAsciiFactory<T>::Create() {
      
      mEvent.reset(new T);
      
      while(std::getline(*mInput, mLine).good()) {
//         std::cout << mLine << std::endl;
         if(AtEndOfEvent()) {
            FinishEvent();
            break;
         } // if
         else if('0' == getFirstNonBlank(mLine)) {
            mEvent->Parse(mLine);
         }
         else if('=' not_eq getFirstNonBlank(mLine)) {
            AddParticle();
         }
      } // if
      
      if(not *mInput) {
         mEvent.reset(NULL);
      } // if
      
      return mEvent.release();
   }
   
   
   template<typename T>
   Int_t
   EventFromAsciiFactory<T>::FinishEvent() {
      
      if(not mEvent->Compute()) {
         std::cerr << "Event computation of event kinematics failed"
         << std::endl;
      } // if
      
      return 0;
   }
   
   
   template<typename T>
   bool EventFromAsciiFactory<T>::AddParticle() {
      
      std::auto_ptr<ParticleMC> particle(new ParticleMC(mLine));
//      particle->ComputeDerivedQuantities();
      
      // That's a keeper!
      particle->SetEvent(mEvent.get());
      mEvent->AddLast(particle.release());
      
      return true;
   }
   
} // namespace erhic

namespace {
   
   // Need this to generate the code for each version
   erhic::EventFromAsciiFactory<EventDjangoh> ed;
   erhic::EventFromAsciiFactory<EventDpmjet> ej;
   erhic::EventFromAsciiFactory<EventPepsi> ee;
   erhic::EventFromAsciiFactory<EventMilou> em;
   erhic::EventFromAsciiFactory<EventRapgap> er;
   erhic::EventFromAsciiFactory<erhic::EventPythia> ep;
   
}
