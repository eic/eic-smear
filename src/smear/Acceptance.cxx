//
// Acceptance.cxx
//
// Created by TB on 8/15/11.
// Copyright 2011 BNL. All rights reserved.
//

#include "eicsmear/smear/Acceptance.h"

#include <TLorentzVector.h>

namespace Smear {
   Acceptance::Acceptance(int genre)
   : mGenre(genre) {
   }
   void Acceptance::AddZone(const Zone& z) {
      mZones.push_back(z);
   }
   void Acceptance::SetGenre(int n) {
      if(n > 0 and n < 3) {
         mGenre = n;
      } // if
      else {
         mGenre = 0;
      } // else
   }
   void Acceptance::AddCustomAcceptance(TString s, double min, double max, int n) {
      KinType kin1 = kP;
      KinType kin2 = kTheta;
      CustomCut C;
      C.dim = ParseInputFunction(s, kin1, kin2);
      std::cout << "Added custom cut " << s << std::endl;
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
   void Acceptance::AddParticle(int n) {
      mParticles.insert(n);
   }
   void Acceptance::ClearParticles() {
      mParticles.clear();
   }
   void Acceptance::RemoveParticle(int n) {
      mParticles.erase(n);
   }
   bool Acceptance::IsCustomAccepted(CustomCut& C, const erhic::VirtualParticle& prt) const {
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
   bool Acceptance::Is(const erhic::VirtualParticle& prt) {
      bool b = false;
      // Check for genre first (em, hadronic, any)
      if(PGenre(prt) == 0 or (mGenre not_eq 0 and PGenre(prt) not_eq mGenre)) {
         return b;
      } // if
      // If there are no Zones, accept everything that passed genre check
      if(mZones.empty()) {
         return true;
      } // if
      // We have zones to check, so look at each in turn.
      // Accept particles if they fall in any zone.
      bool intheta;
      bool inphi;
      bool inE;
      bool inp; 
      bool inpt;
      bool inpz;
      bool inCustom=true;
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
         TLorentzVector p = prt.Get4Vector();
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
      if(not mParticles.empty()) {
         b = b and mParticles.count(prt.Id()) > 0;
      } // if
      return b;
   }
   //
   // class Acceptance::CustomCut
   //
   Acceptance::CustomCut::CustomCut()
   : dim(0)
   , Kin1(kP)
   , Kin2(kTheta) {
   }
   //
   // class Acceptance::Zone
   //
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
   bool Acceptance::Zone::IsCustomAccepted(Acceptance::CustomCut C, const erhic::VirtualParticle& prt) const {
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
   Bool_t Acceptance::Zone::Contains(erhic::VirtualParticle& prt) const {
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
}
