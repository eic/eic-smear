//
// Acceptance.cxx
//
// Created by TB on 8/15/11.
// Copyright 2011 BNL. All rights reserved.
//

#include "Acceptance.h"

#include <TLorentzVector.h>

namespace Smear {

   Acceptance::Acceptance(int genre)
   : Genre(genre) {
      mZones.push_back(Zone());
   }

   void Acceptance::AddZone(const Zone& z) {
      mZones.push_back(z);
   }

   void Acceptance::RemoveZone(unsigned n) {
      if(n < mZones.size()) {
         mZones.erase(mZones.begin() + n);
      } // if
   }

   void Acceptance::SetTheta(double min, double max, int n) {
      mZones.at(n).thetaMin = FixTopologyTheta(min);
      mZones.at(n).thetaMax = FixTopologyTheta(max);
   }

   void Acceptance::SetPhi(double min, double max, int n) {
      mZones.at(n).phiMin = FixTopologyPhi(min);
      mZones.at(n).phiMax = FixTopologyPhi(max);
   }

   void Acceptance::SetE(double min, double max, int n) {
      mZones.at(n).EMin = min;
      mZones.at(n).EMax = max;
   }

   void Acceptance::SetP(double min, double max, int n) {
      mZones.at(n).PMin = min;
      mZones.at(n).PMax = max;	
   }

   void Acceptance::Set(KinType type, double min, double max, int n) {
      switch(type) {
         case kE: 
            SetE(min,max,n); break;
         case kP:
            SetP(min,max,n); break;
         case kTheta:
            SetTheta(min,max,n); break;
         case kPhi:
            SetPhi(min,max,n); break;
         case kPz:
            SetPz(min, max, n); break;
         case kPt:
            SetPt(min, max, n); break;
            break; // Do nothing
      } // switch
   }

   void Acceptance::SetPt(double min, double max, int n) {
      mZones.at(n).pTMin = min;
      mZones.at(n).pTMax = max;
   }

   void Acceptance::SetPz(double min, double max, int n) {			
      mZones.at(n).pZMin = min;
      mZones.at(n).pZMax = max;
   }

   void Acceptance::SetGenre(int n) {
      if(n > 0 and n < 3) {
         Genre = n;
      } // if
      else {
         Genre = 0;
      } // else
   }

   void Acceptance::AddCustomAcceptance(TString s, double min, double max, int n) {
      KinType kin1=kP;
      KinType kin2=kTheta;
      CustomCut C;
      C.dim = ParseInputFunction(s, kin1, kin2);
      if(not IsCoreType(kin1) or not IsCoreType(kin2)) {
         std::cerr <<
         "ERROR! Custom acceptance is not a function of E, p, theta, phi"
         << std::endl;
         return;
      } // if
      
      double min1=0.;
      double max1=1e12;
      double min2=0.;
      double max2=1e12;
      
      if(kin1==kTheta) {
         max1 = pi;
      } // if
      else if(kin1==kPhi) {
         max1 = 2.*pi;
      } // else if
      
      if(kin2==kTheta) {
         max2 = pi;
      } // if
      else if(kin2==kPhi) {
         max2 = 2.*pi;
      } // else if
      
      if(1 == C.dim) {
         C.F1 = TF1("f1",s,min1,max1);
      } // if
      else if(2 == C.dim) {
         C.F2 = TF2("f2",s,min1,min2,max1,max2);
      } // else if
      else {
         std::cerr <<
         "ERROR! Provided custom acceptance is not of dimension 1 or 2."
         << std::endl;
         return;
      } // else if
      
      C.Min = min;
      C.Max = max;
      
      C.Kin1 = kin1;
      C.Kin2 = kin2;
      
      mZones.at(n).CustomCuts.push_back(C);
   }

   void Acceptance::ClearCustomAcceptance(int n) {
      for(unsigned i(0); i < mZones.at(n).CustomCuts.size(); ++i) {
         delete &(mZones.at(n).CustomCuts.at(i));
      } // for
      mZones.at(n).CustomCuts.clear();
   }

   void Acceptance::AddParticle(int n) {
      Particles.insert(n);
   }

   void Acceptance::ClearParticles() {
      Particles.clear();
   }

   void Acceptance::RemoveParticle(int n) {
      Particles.erase(n);
   }

   bool Acceptance::IsCustomAccepted(CustomCut& C, const erhic::ParticleMC& prt) const {
      bool accept(false);
      
      double x = SwitchKinGetFromParticle(prt, C.Kin1);
      double y(NAN), z(NAN);
      
      switch(C.dim) {
         case 1:
            z = C.F1.Eval(x);
            accept = z > C.Min and z < C.Max;
            break;
         case 2:
            y = SwitchKinGetFromParticle(prt, C.Kin2);
            z = C.F2.Eval(x,y);
            accept = z > C.Min and z < C.Max;
            break;
         default:
            break;
      } // switch
      
      return accept;
   }

   bool Acceptance::Zone::IsCustomAccepted(Acceptance::CustomCut C, const erhic::ParticleMC& prt) const {
      bool accept(false);
      double x = SwitchKinGetFromParticle(prt, C.Kin1);
      double y(NAN), z(NAN);
      switch(C.dim) {
         case 1:
            z = C.F1.Eval(x);
            accept = z > C.Min and z < C.Max;
            break;
         case 2:
            y = SwitchKinGetFromParticle(prt, C.Kin2);
            z = C.F2.Eval(x,y);
            accept = z > C.Min and z < C.Max;
            break;
         default:
            break;
      } // switch
      
      return accept;
   }

   Bool_t Acceptance::Zone::Contains(erhic::ParticleMC& prt) const {
      bool accept(true);
      
      const double theta = FixTopologyTheta(prt.GetTheta());
      const double phi = FixTopologyPhi(prt.GetPhi());
      
      if(theta < thetaMin or theta > thetaMax) {
         accept = false;
      } // if...
      else if(phi < phiMin or phi > phiMax) {
         accept = false;
      } // else if...
      else if(prt.GetE() < EMin or prt.GetE() > EMax) {
         accept = false;
      } // else if...
      else if(prt.GetP() < PMin or prt.GetP() > PMax) {
         accept = false;
      } // else if
      
      for(unsigned j(0); j < CustomCuts.size(); ++j) {
         if(not IsCustomAccepted(CustomCuts.at(j), prt)) {
            accept = false;
            break;
         } // if
      } // for
      
      return accept;
   }

   bool Acceptance::Is(const erhic::ParticleMC& prt) {
      bool b = false;
      bool intheta;
      bool inphi;
      bool inE;
      bool inp; 
      bool inpt;
      bool inpz;
      bool inCustom=true;
      if(PGenre(prt)==0 or (Genre not_eq 0 and PGenre(prt) not_eq Genre)) {
         return b;
      } // if
      for(unsigned i(0); i < mZones.size(); i++) {
         double theta = prt.GetTheta();
         double phi = prt.GetPhi();
         intheta =
            theta >= mZones.at(i).thetaMin and theta <= mZones.at(i).thetaMax;
         inphi   =
            phi >= mZones.at(i).phiMin and phi <= mZones.at(i).phiMax;
         inE     =
            prt.GetE() >= mZones.at(i).EMin and prt.GetE() <= mZones.at(i).EMax;
         inp     =
            prt.GetP() >= mZones.at(i).PMin and prt.GetP() <= mZones.at(i).PMax;
         TLorentzVector p = prt.PxPyPzE();
         inpz = p.Pz() >= mZones.at(i).pZMin and p.Pz() <= mZones.at(i).pZMax;
         inpt = p.Pt() >= mZones.at(i).pTMin and p.Pt() <= mZones.at(i).pTMax;
         // Loop through custom cuts in zone, if there are any
         if(mZones.at(i).CustomCuts.size() not_eq 0) {
            for(unsigned j=0; j<mZones.at(i).CustomCuts.size(); j++) {
               inCustom = IsCustomAccepted(mZones.at(i).CustomCuts.at(j),prt);
            } // for
         } // if
         if(intheta and inphi and inE and inp and inpz and inpt and inCustom) {
            b = true;
         } // if
      } // for
      // Check against exclusive particle list
      if(not Particles.empty()) {
         b = b and Particles.count(prt.Id()) > 0;
      } // if
      return b;
   }

   Acceptance::CustomCut::CustomCut()
   : dim(0)
   , Kin1(kP)
   , Kin2(kTheta) {
   }

   Acceptance::Zone::Zone(double thMin, double thMax,
                          double phMin, double phMax,
                          double eMin, double eMax,
                          double pMin, double pMax,
                          double ptmin, double ptmax,
                          double pzmin, double pzmax)
   : thetaMin(thMin)
   , thetaMax(thMax)
   , phiMin(phMin)
   , phiMax(phMax)
   , EMin(eMin)
   , EMax(eMax)
   , PMin(pMin)
   , PMax(pMax)
   , pTMin(ptmin)
   , pTMax(ptmax)
   , pZMin(pzmin)
   , pZMax(pzmax) {
   }
}
