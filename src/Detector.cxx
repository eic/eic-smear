//
// Detector.cxx
//
// Created by TB on 8/15/11.
// Copyright 2011 BNL. All rights reserved.
//

#include "Device.h"
#include "EventBase.h"
#include "ParticleID.h"

#include "Detector.h"

namespace Smear {
   
   Detector::Detector()
//   : NDevices(0)
   : useNM(false)
   , useJB(false)
   , useDA(false)
   , NIdentifiers(0)
   {}
   
   void Detector::DeleteAllDevices() {
      for(unsigned i=0; i<GetNDevices(); i++) {
         delete Devices.at(i);
      } // for
      Devices.clear();
      //NDevices=Devices.size();
   }
   
   void Detector::RemoveAllDevices() {
      Devices.clear();
      //NDevices=Devices.size();
   }
   
   void Detector::DeleteAllIdentifiers() {
      for(int i=0; i<NIdentifiers; i++) {
         delete Identifiers.at(i);
      }
      Identifiers.clear();
      NIdentifiers=Identifiers.size();
   }
   
   void Detector::RemoveAllIdentifiers() {
      Identifiers.clear();
      NIdentifiers=Identifiers.size();
   }
   
   void Detector::AddDevice(Device* dev) {
      Devices.push_back(dev);
      //NDevices = Devices.size();
   }
   
   void Detector::AddDevice(Device& dev) {
      Devices.push_back(dev.Clone());
      //NDevices = Devices.size();
   }
   
   void Detector::AddDevice(ParticleID* ident) {
      Identifiers.push_back(ident);
      NIdentifiers = Identifiers.size();
   }
   
   void Detector::AddDevice(ParticleID& ident) {
      Identifiers.push_back(ident.Clone());
      NIdentifiers = Identifiers.size();
   }
   
   void Detector::SetPDGLeptonCode(int n) {
      EventKinComp.SetPDGLeptonCode(n);
   }
   
   void Detector::SetMasterAccept(KinType type, double min, double max) {
      for(unsigned i=0; i<GetNDevices(); i++) {
         Devices.at(i)->Accept.Set(type,min,max);
      }
   }
   
   void Detector::SetMasterParametrization(TString param) {
      for(unsigned i=0; i<GetNDevices(); i++) {
         Devices.at(i)->SetParametrization(param);
      }
   }
   
   void Detector::SetMasterRanSeed(int n) {
      for(unsigned i=0; i<GetNDevices(); i++) {
         Devices.at(i)->SetRanSeed(n);
      }
   }
   
   void Detector::SetMasterParams(int n, double xi) {
      for(unsigned i=0; i<GetNDevices(); i++) {
         Devices.at(i)->SetParams(n,xi);
      }
   }
   
   void Detector::SetMasterDistribution(TString f) {
      for(unsigned i=0; i<GetNDevices(); i++) {
         Devices.at(i)->SetDistribution(f);
      }
   }
   
   void Detector::SetMasterDistributionRange(double min, double max) {
      for(unsigned i=0; i<GetNDevices(); i++) {
         Devices.at(i)->SetDistributionRange(min,max);
      }
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
   
   Device* Detector::GetDevice(int n) {
      // TODO Protect against out-of-range
      return Devices.at(n);
   }
   
   ParticleID* Detector::GetIdentifier(int n) {
      // TODO Protect against out-of-range
      return Identifiers.at(n);
   }
   
   Device* Detector::GetDeviceUsedFor(Particle prt, KinType kin) {
      Device* dev(NULL);
      for(unsigned i = 0; i < GetNDevices(); i++) {
         if(GetDevice(i)->Accept.Is(prt) and GetDevice(i)->GetTargetVariable()==kin) {
            dev = GetDevice(i);
         } // if
      } // for
      return dev;
   }
   
   int Detector::WhichDeviceUsedFor(Particle prt, KinType kin) {
      int l(-1);
      for(unsigned i=0; i<GetNDevices(); i++) {
         if(GetDevice(i)->Accept.Is(prt) and GetDevice(i)->GetTargetVariable()==kin) {
            l=i;
         } // if
      } // for
      return l;
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
//      std::cout << "Detector::DetSmear()" << std::endl;
      //check if particle is stable
      // TODO Particles should also only be created if they lie within the
      // detector acceptance.
      bool accept(false);
      for(unsigned i(0); i < Devices.size(); ++i) {
         if(Devices.at(i)->Accept.Is(prt)) {
//            std::cout << "\tIn acceptance of Device " << Devices.at(i)->name << std::endl;
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
//            std::cout << "\tSmearing with device " << i << " " << Devices.at(i)->name << std::endl;
            Devices.at(i)->DevSmear(prt,*prtOut);
         } // for
         
         if(NIdentifiers not_eq 0) {
            for(int i=0; i<NIdentifiers; i++) {
//               std::cout << "\tSmearing with device " << Devices.at(i)->name << std::endl;
               Identifiers.at(i)->DevSmear(prt,*prtOut);
            } // for
         } // if
         prtOut->px = prtOut->p*sin(prtOut->theta)*cos(prtOut->phi);
         prtOut->py = prtOut->p*sin(prtOut->theta)*sin(prtOut->phi);
         prtOut->pz = prtOut->p*cos(prtOut->theta);
      } // if
      
      return prtOut;
   }
   
} // namespace Smear
