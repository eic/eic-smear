//
// EventFactory.cxx
//
// Created by TB on 11/1/11.
// Copyright 2011 BNL. All rights reserved.
//

#include <memory>
#include <stdexcept>

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
   return not mLine.empty()
          and mLine.find("finished") not_eq std::string::npos;
}

template<typename T>
T* EventFromAsciiFactory<T>::Create() {
   // Save current object count
   int objectNumber = TProcessID::GetObjectCount();
   mEvent.reset(new T);
   // We use this flag to check input doesn't end mid-event.
   // Initialised finished flag to "success" in case of no input.
   int finished(0);
   while(std::getline(*mInput, mLine).good()) {
      if(AtEndOfEvent()) {
         finished = FinishEvent(); // Zero upon success
         break;
      } // if
      else if('0' == getFirstNonBlank(mLine)) {
         // An event started, set finished flag to "unfinished".
         finished = -1;
         // Parse string and check for validity.
         if(not mEvent->Parse(mLine)) {
            mEvent.reset(NULL);
            throw std::runtime_error("Bad event input: " + mLine);
         } // if
      } // else if
      else if('=' not_eq getFirstNonBlank(mLine)) {
         AddParticle(); // Will throw upon bad input
      } // else if
   } // if
   // Check for events that started but did not finish.
   if(finished not_eq 0) {
      mEvent.reset(NULL);
      throw std::runtime_error("Ended mid-event");
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
   std::auto_ptr<DisKinematics> nm(
      LeptonKinematicsComputer(*mEvent).Calculate());
   std::auto_ptr<DisKinematics> jb(
      JacquetBlondelComputer(*mEvent).Calculate());
   std::auto_ptr<DisKinematics> da(
      DoubleAngleComputer(*mEvent).Calculate());
   if(nm.get()) {
      mEvent->SetLeptonKinematics(*nm);
   } // if
   for(unsigned n(0); n < mEvent->GetNTracks(); ++n) {
      mEvent->GetTrack(n)->ComputeEventDependentQuantities(*mEvent);
   } // for
   if(jb.get()) {
      mEvent->SetJacquetBlondelKinematics(*jb);
   } // if
   if(da.get()) {
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
   ParticleMC particle(mLine); // Throws if the string is bad
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

} // namespace
