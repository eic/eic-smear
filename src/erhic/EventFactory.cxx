//
// EventFactory.cxx
//
// Created by TB on 11/1/11.
// Copyright 2011 BNL. All rights reserved.
//

#include <memory>

#include <TClass.h>
#include <TProcessID.h>

#include "eicsmear/erhic/BeamParticles.h"
#include "eicsmear/erhic/EventFactory.h"
#include "eicsmear/erhic/EventPythia.h"
#include "eicsmear/erhic/EventMilou.h"
#include "eicsmear/erhic/EventDjangoh.h"
#include "eicsmear/erhic/EventDpmjet.h"
#include "eicsmear/erhic/EventRapgap.h"
#include "eicsmear/erhic/EventPepsi.h"
#include "eicsmear/erhic/EventGmcTrans.h"
#include "eicsmear/functions.h" // For getFirstNonBlank()
#include "eicsmear/erhic/Kinematics.h"
#include "eicsmear/erhic/ParticleIdentifier.h"
#include "eicsmear/erhic/ParticleMC.h"

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
      DisKinematics* nm = LeptonKinematicsComputer(*mEvent).Calculate();
      DisKinematics* jb = JacquetBlondelComputer(*mEvent).Calculate();
      DisKinematics* da = DoubleAngleComputer(*mEvent).Calculate();
      if(nm) {
         mEvent->SetLeptonKinematics(*nm);
      } // if
      for(unsigned n(0); n < mEvent->GetNTracks(); ++n) {
         mEvent->GetTrack(n)->ComputeEventDependentQuantities(*mEvent);
      } // for
      if(jb) {
         mEvent->SetJacquetBlondelKinematics(*jb);
      } // if
      if(da) {
         mEvent->SetDoubleAngleKinematics(*da);
      } // if
      // We also have to set the remaining variables not taken care of
      // by the general DIS event kinematic computations.
      // Find the beams, exchange boson, scattered lepton.
      BeamParticles beams;
      if(not ParticleIdentifier::IdentifyBeams(*mEvent, beams)) {
         std::cerr <<
         "EventFromAsciiFactory::FinishEvent(): failed to find beams"
         << std::endl;
         return 1;
      } // if
      const TLorentzVector h = beams.BeamHadron();
      TLorentzVector l = beams.BeamLepton();
      TLorentzVector s = beams.ScatteredLepton();
      TVector3 boost = -h.BoostVector();
      l.Boost(boost);
      s.Boost(boost);
      mEvent->SetELeptonInNuclearFrame(l.E());
      mEvent->SetEScatteredInNuclearFrame(s.E());
      return 0;
   }

   template<typename T>
   bool EventFromAsciiFactory<T>::AddParticle() {
      ParticleMC particle(mLine);
      particle.SetEvent(mEvent.get());
      mEvent->AddLast(&particle);
      return true;
   }

   template<typename T>
   std::string EventFromAsciiFactory<T>::EventName() const {
      return T::Class()->GetName();
   }
} // namespace erhic

namespace {
   // Need this to generate the CINT code for each version
   erhic::EventFromAsciiFactory<erhic::EventDjangoh> ed;
   erhic::EventFromAsciiFactory<erhic::EventDpmjet> ej;
   erhic::EventFromAsciiFactory<erhic::EventPepsi> ee;
   erhic::EventFromAsciiFactory<erhic::EventMilou> em;
   erhic::EventFromAsciiFactory<erhic::EventRapgap> er;
   erhic::EventFromAsciiFactory<erhic::EventPythia> ep;
   erhic::EventFromAsciiFactory<erhic::EventGmcTrans> eg;
}
