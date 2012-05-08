//
// Acceptance.cxx
//
// Created by TB on 8/15/11.
// Copyright 2011 BNL. All rights reserved.
//

#include "eicsmear/smear/Acceptance.h"

#include <TLorentzVector.h>
#include <TString.h>

namespace Smear {

   Acceptance::~Acceptance() {
   }

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

   void Acceptance::AddParticle(int n) {
      mParticles.insert(n);
   }

   bool Acceptance::Is(const erhic::VirtualParticle& prt) const {
      bool b = false;
      // Check for genre first (em, hadronic, any)
      if(PGenre(prt) == 0 or (mGenre not_eq 0 and PGenre(prt) not_eq mGenre)) {
         return b;
      } // if
      // If there are no Zones, accept everything that passed genre check
      if(mZones.empty()) {
         return true;
      } // if
      for(unsigned i(0); i < mZones.size(); i++) {
         if(mZones.at(i).Contains(prt)) {
            b = true;
            break;
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

   Acceptance::CustomCut::~CustomCut() {
   }

   Acceptance::CustomCut::CustomCut()
   : mFormula("CustomCutFormula", "0")
   , dim(0)
   , Kin1(kP)
   , Kin2(kTheta)
   , Min(-TMath::Infinity())
   , Max(TMath::Infinity()) {
   }

   Acceptance::CustomCut::CustomCut(const TString& formula,
                                    double min, double max)
   : mFormula("CustomCutFormula", "0")
   , dim(0)
   , Kin1(kP)
   , Kin2(kTheta)
   , Min(min)
   , Max(max) {
      if(not IsCoreType(Kin1) or not IsCoreType(Kin2)) {
         std::cerr <<
         "ERROR! Custom acceptance is not a function of E, p, theta, phi"
         << std::endl;
      } // if
      TString s(formula);
      dim = ParseInputFunction(s, Kin1, Kin2);
      if(1 == dim or 2 == dim) {
         mFormula = TFormula("CustomCutFormula", s);
      } // if
      else {
         std::cerr <<
         "ERROR! Provided custom acceptance is not of dimension 1 or 2."
         << std::endl;
         return;
      } // else if
      std::cout << "Added custom cut " << formula << std::endl;
   }

   bool Acceptance::CustomCut::Contains(
                                 const erhic::VirtualParticle& prt) const {
      double x = SwitchKinGetFromParticle(prt, Kin1);
      double y(0.);
      if(2 == dim) {
         y = SwitchKinGetFromParticle(prt, Kin2);
      } // if
      double z = mFormula.Eval(x, y);
      return z >= Min and z < Max;
   }

   //
   // class Acceptance::Zone
   //

   Acceptance::Zone::~Zone() {
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

   void Acceptance::Zone::Add(const CustomCut& cut) {
      CustomCuts.push_back(cut);
   }

   Bool_t Acceptance::Zone::Contains(const erhic::VirtualParticle& prt) const {
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
      else if(prt.GetPz() < pZMin or prt.GetPz() > pZMax) {
         accept = false;
      } // else if
      else if(prt.GetPt() < pTMin or prt.GetPt() > pTMax) {
         accept = false;
      } // else if
      // If it made it this far, test the custom cut(s)
      if(accept) {
         for(unsigned j(0); j < CustomCuts.size(); ++j) {
            if(not CustomCuts.at(j).Contains(prt)) {
               accept = false;
               break;
            } // if
         } // for
      } // if
      return accept;
   }
}
