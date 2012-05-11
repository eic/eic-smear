//
// Detector.cxx
//
// Created by TB on 8/15/11.
// Copyright 2011 BNL. All rights reserved.
//

#include "eicsmear/smear/Detector.h"

#include <algorithm>
#include <functional>
#include <memory>

#include "eicsmear/erhic/EventDis.h"
#include "eicsmear/smear/EventSmear.h"
#include "eicsmear/erhic/Kinematics.h"
#include "eicsmear/smear/ParticleMCS.h"
#include "eicsmear/smear/Smearer.h"
#include "eicsmear/erhic/VirtualParticle.h"

namespace Smear {
   Detector::Detector()
   : useNM(false)
   , useJB(false)
   , useDA(false) {
   }
   Detector::Detector(const Detector& other)
   : TObject(other) {
      useNM = other.useNM;
      useJB = other.useJB;
      useDA = other.useDA;
      Devices = other.CopyDevices();
   }
   Detector& Detector::operator=(const Detector& that) {
      if(this not_eq &that) {
         useNM = that.useNM;
         useJB = that.useJB;
         useDA = that.useDA;
         Devices = that.CopyDevices();
      } // if
      return *this;
   }
   Detector::~Detector() {
      DeleteAllDevices();
   }
   void Detector::DeleteAllDevices() {
      for(unsigned i(0); i < GetNDevices(); i++) {
         delete Devices.at(i);
         Devices.at(i) = NULL;
      } // for
      Devices.clear();
   }
   void Detector::AddDevice(Smearer& dev) {
      Devices.push_back(dev.Clone());
   }
   void Detector::SetEventKinematicsCalculator(TString s) {
      s.ToLower();
      useNM = s.Contains("nm") or s.Contains("null");
      useJB = s.Contains("jb") or s.Contains("jacquet");
      useDA = s.Contains("da") or s.Contains("double");
   }
   Smearer* Detector::GetDevice(int n) {
      Smearer* smearer(NULL);
      if(unsigned(n) < Devices.size()) {
         smearer = Devices.at(n);
      } // if
      return smearer;
   }
   void Detector::FillEventKinematics(const erhic::EventDis& event,
                                      Event* eventS) {
      if(not (useNM or useJB or useDA)) {
         return;
      } // if
      // Need a bit of jiggery-pokery here, as the incident beam info isn't
      // associated with the smeared event.
      // So, get the beam info from the MC event, but replace the scattered
      // electron with the smeared version.
      // Then we can use the standard JB/DA algorithms on the smeared event.
      BeamParticles beams;
      ParticleIdentifier::IdentifyBeams(event, beams);
      const ParticleMCS* scattered = eventS->ScatteredLepton();
      if(scattered) {
         beams.SetScatteredLepton(scattered->Get4Vector());
      } // if
      typedef std::auto_ptr<erhic::DisKinematics> KinPtr;
      if(useNM and scattered) {
         KinPtr kin(erhic::LeptonKinematicsComputer(beams).Calculate()); 
         if(kin.get()) {
            eventS->SetLeptonKinematics(*kin);
         } // if
      } // if
      else {
        eventS->SetLeptonKinematics(
           erhic::DisKinematics(-1., -1., -1., -1., -1.));
      } // else
      if(useJB) {
         KinPtr kin(erhic::JacquetBlondelComputer(*eventS, &beams).Calculate());
         if(kin.get()) {
            eventS->SetJacquetBlondelKinematics(*kin);
         } // if
      } // if
      if(useDA and scattered) {
         KinPtr kin(erhic::DoubleAngleComputer(*eventS, &beams).Calculate());
         if(kin.get()) {
            eventS->SetDoubleAngleKinematics(*kin);
         } // if
      } // if
   }
   ParticleMCS* Detector::Smear(const erhic::VirtualParticle& prt) const {
      // Does the particle fall in the acceptance of any device?
      // If so, we smear it, if not, we skip it (store a NULL pointer).
      bool accept(false);
      for(unsigned i(0); i < Devices.size(); ++i) {
         if(Devices.at(i)->Accept.Is(prt)) {
            accept = true;
            break;
         } // if
      } // for
      ParticleMCS* prtOut(NULL);
      if(prt.GetStatus() == 1 and accept) {
         // It passes through at least one device, so smear it.
         // Devices in which it doesn't pass won't smear it.
         prtOut = new ParticleMCS();
         for(unsigned i=0; i<GetNDevices(); i++) {
            Devices.at(i)->Smear(prt,*prtOut);
         } // for
         prtOut->px = prtOut->p * sin(prtOut->theta) * cos(prtOut->phi);
         prtOut->py = prtOut->p * sin(prtOut->theta) * sin(prtOut->phi);
         prtOut->pt = sqrt(pow(prtOut->px, 2.) + pow(prtOut->py, 2.));
         prtOut->pz = prtOut->p * cos(prtOut->theta);
      } // if
      return prtOut;
   }
   std::vector<Smearer*> Detector::CopyDevices() const {
      std::vector<Smearer*> copies;
      std::transform(Devices.begin(), Devices.end(),
                     std::back_inserter(copies),
                     std::bind2nd(std::mem_fun(&Smearer::Clone), ""));
      return copies;
   }
} // namespace Smear
