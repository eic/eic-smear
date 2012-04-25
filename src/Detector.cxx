//
// Detector.cxx
//
// Created by TB on 8/15/11.
// Copyright 2011 BNL. All rights reserved.
//

#include "Detector.h"

#include <algorithm>
#include <functional>

#include "Device.h"
#include "EventBase.h"
#include "ParticleID.h"

namespace Smear {

   Detector::Detector()
   : useNM(false)
   , useJB(false)
   , useDA(false)
   , NIdentifiers(0) {
   }

   Detector::Detector(const Detector& other)
   : TObject(other) {
      useNM = other.useNM;
      useJB = other.useJB;
      useDA = other.useDA;
      EventKinComp = other.EventKinComp;
      NIdentifiers = other.NIdentifiers;
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

   void Detector::DeleteAllIdentifiers() {
      for(int i(0); i < NIdentifiers; i++) {
         delete Identifiers.at(i);
         Identifiers.at(i) = NULL;
      } // for
      Identifiers.clear();
      NIdentifiers=Identifiers.size();
   }

   void Detector::AddDevice(Smearer& dev) {
      Devices.push_back(dev.Clone());
   }

   void Detector::AddDevice(ParticleID& ident) {
      Identifiers.push_back(ident.Clone());
      NIdentifiers = Identifiers.size();
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
      if(s.Contains("NM") or s.Contains("null momentum") or s.Contains("zero mass")) {
         useNM = true;
      } else useNM=false;
      if(s.Contains("JB") or s.Contains("Jacquet") or s.Contains("Blondel")) {
         useJB = true;
      } else useJB=false;
      if(s.Contains("DA") or s.Contains("double angle")) {
         useDA = true;
      } else useDA=false;
   }
   
   void Detector::RemoveDevice(int n) {
      // TODO Protect against out-of-range
      delete Devices.at(n);
      Devices.erase(Devices.begin()+n);
      //NDevices = Devices.size();
   }
   
   void Detector::RemoveIdentifier(int n) {
      // TODO Protect against out-of-range
      delete Identifiers.at(n);
      Identifiers.erase(Identifiers.begin()+n);
      NIdentifiers = Identifiers.size();
   }
   
   Smearer* Detector::GetDevice(int n) {
      // TODO Protect against out-of-range
      return Devices.at(n);
   }
   
   ParticleID* Detector::GetIdentifier(int n) {
      // TODO Protect against out-of-range
      return Identifiers.at(n);
   }

   void Detector::FillEventKinematics(const EventBase *event, EventS *eventS) {
      if(not (useNM or useJB or useDA)) {
         return;
      } // if
      
      if (not (useJB or useDA)) {
         EventKinComp.SetNeedBackground(false);
      } // if
      else {
         EventKinComp.SetNeedBackground(true);
      } // else
      
      EventKinComp.ReadEvent(event,eventS);
      
      if(useNM) {
         EventKinComp.ComputeUsingNM();
         eventS->QSquared = EventKinComp.QSquared();
         eventS->y = EventKinComp.y();
         eventS->x = EventKinComp.x();;
         eventS->WSquared = EventKinComp.WSquared();
      } // if
      
      if(useJB) {
         EventKinComp.ComputeUsingJB(/*eventS*/);
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
      ParticleS *prtOut = NULL;
      //check if particle is stable
      // TODO Particles should also only be created if they lie within the
      // detector acceptance.
      bool accept(false);
      for(unsigned i(0); i < Devices.size(); ++i) {
         if(Devices.at(i)->Accept.Is(prt)) {
            accept = true;
            break;
         } // if
      } // for
      if(NIdentifiers not_eq 0 and not accept) {
         for(unsigned i(0); i < Identifiers.size(); ++i) {
            if(Identifiers.at(i)->Accept.Is(prt)) {
               accept = true;
               break;
            } // if
         } // for
      } // if
      
      if(prt.GetStatus() == 1 and accept) {
         
         prtOut = new ParticleS();
         
         for(unsigned i=0; i<GetNDevices(); i++) {
            Devices.at(i)->DevSmear(prt,*prtOut);
         } // for
         
         for(int i=0; i<NIdentifiers; i++) {
            Identifiers.at(i)->DevSmear(prt,*prtOut);
         } // for
         prtOut->px = prtOut->p*sin(prtOut->theta)*cos(prtOut->phi);
         prtOut->py = prtOut->p*sin(prtOut->theta)*sin(prtOut->phi);
         prtOut->pz = prtOut->p*cos(prtOut->theta);
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
