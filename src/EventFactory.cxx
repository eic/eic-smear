//
// EventFactory.cxx
//
// Created by TB on 11/1/11.
// Copyright 2011 BNL. All rights reserved.
//

#include <memory>

#include <TProcessID.h>

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
      // Save current object count
      int objectNumber = TProcessID::GetObjectCount();

      mEvent.reset(new T);
      while(std::getline(*mInput, mLine).good()) {
         if(AtEndOfEvent()) {
            FinishEvent();
            break;
         } // if
         else if('0' == getFirstNonBlank(mLine)) {
            mEvent->Parse(mLine);
         } // else if
         else if('=' not_eq getFirstNonBlank(mLine)) {
            AddParticle();
         } // else if
      } // if
      if(not *mInput) {
         mEvent.reset(NULL);
      } // if

      // Restore Object count 
      // See example in $ROOTSYS/test/Event.cxx
      // To save space in the table keeping track of all referenced objects
      // we assume that our events do not address each other. We reset the 
      // object count to what it was at the beginning of the event.
      TProcessID::SetObjectCount(objectNumber);

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
