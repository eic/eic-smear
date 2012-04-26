//
// Detector.cxx
//
// Created by TB on 8/15/11.
// Copyright 2011 BNL. All rights reserved.
//

#include "Detector.h"

#include <algorithm>
#include <functional>

#include "EventMC.h"
#include "EventSmear.h"
#include "ParticleMC.h"
#include "ParticleMCS.h"
#include "Smearer.h"

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
      EventKinComp = other.EventKinComp;
      Devices = other.CopyDevices();
   }

   Detector::~Detector() {
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

   void Detector::SetPDGLeptonCode(int n) {
      EventKinComp.SetPDGLeptonCode(n);
   }

   void Detector::SetMissingEnergyTolerance(double E) {
      EventKinComp.SetMissingEnergyTolerance(E);
   }
   
   void Detector::TolerateBadEvents(double b) {
      EventKinComp.TolerateBadEvents(b);
   }
   
   void Detector::SetSupressEventWarnings(bool b) {
      EventKinComp.SetSupressWarnings(b);
   }
   
   void Detector::SetEventKinematicsCalculator(TString s) {
      s.ToLower();
      useNM = s.Contains("nm") or s.Contains("null");
      useJB = s.Contains("jb") or s.Contains("jacquet");
      useDA = s.Contains("da") or s.Contains("double");
   }
   
   void Detector::RemoveDevice(int n) {
      if(unsigned(n) < Devices.size()) {
         delete Devices.at(n);
         Devices.erase(Devices.begin()+n);
      } // if
   }

   Smearer* Detector::GetDevice(int n) {
      Smearer* smearer(NULL);
      if(unsigned(n) < Devices.size()) {
         smearer = Devices.at(n);
      } // if
      return smearer;
   }

   void Detector::FillEventKinematics(const erhic::EventMC* event,
                                      Event* eventS) {
      if(not (useNM or useJB or useDA)) {
         return;
      } // if

      if(useJB or useDA) {
         EventKinComp.SetNeedBackground(true);
      } // if
      else {
         EventKinComp.SetNeedBackground(false);
      } // else
      EventKinComp.ReadEvent(event, eventS);

      if(useNM) {
         EventKinComp.ComputeUsingNM();
         eventS->QSquared = EventKinComp.QSquared();
         eventS->y = EventKinComp.y();
         eventS->x = EventKinComp.x();;
         eventS->WSquared = EventKinComp.WSquared();
      } // if
      
      if(useJB) {
         EventKinComp.ComputeUsingJB();
         eventS->QSquaredJB = EventKinComp.QSquared();
         eventS->yJB = EventKinComp.y();
         eventS->xJB = EventKinComp.x();
         eventS->WSquaredJB = EventKinComp.WSquared();
      } // if
      
      if(useDA) {
         EventKinComp.ComputeUsingDA();
         eventS->QSquaredDA = EventKinComp.QSquared();
         eventS->yDA = EventKinComp.y();
         eventS->xDA = EventKinComp.x();
         eventS->WSquaredDA = EventKinComp.WSquared();
      } // if
   }

   ParticleS* Detector::DetSmear(const Particle& prt) {
      bool accept(false);
      for(unsigned i(0); i < Devices.size(); ++i) {
         if(Devices.at(i)->Accept.Is(prt)) {
            accept = true;
            break;
         } // if
      } // for
      ParticleS* prtOut(NULL);
      if(prt.GetStatus() == 1 and accept) {
         prtOut = new ParticleS();
         for(unsigned i=0; i<GetNDevices(); i++) {
            Devices.at(i)->DevSmear(prt,*prtOut);
         } // for
         prtOut->px = prtOut->p * sin(prtOut->theta) * cos(prtOut->phi);
         prtOut->py = prtOut->p * sin(prtOut->theta) * sin(prtOut->phi);
         prtOut->pz = prtOut->p * cos(prtOut->theta);
      } // if
      return prtOut;
   }

   std::vector<Smearer*> Detector::CopyDevices() const {
      using namespace std;
      vector<Smearer*> copies;
      transform(Devices.begin(), Devices.end(), back_inserter(copies),
                bind2nd(mem_fun(&Smearer::Clone), ""));
      return copies;
   }
} // namespace Smear
