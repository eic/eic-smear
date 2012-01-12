//
// Acceptance.cxx
//
// Created by TB on 8/15/11.
// Copyright 2011 BNL. All rights reserved.
//

#include "Acceptance.h"

namespace Smear {
   
   Acceptance::Acceptance(int genre)
   : Genre(genre) {
      mZones.push_back(Zone());
   }
   
//   void Acceptance::AddZone() {
//      mZones.push_back(Zone());
//   }
   
   void Acceptance::AddZone(double thetamin, double thetamax,
                            double phimin, double phimax,
                            double Emin, double Emax,
                            double pmin, double pmax) {
      mZones.push_back(Zone(thetamin, thetamax, phimin, phimax,
                           Emin, Emax, pmin, pmax));
   }
   
   void Acceptance::RemoveZone(unsigned n) {
//      std::cout << "Removing zone " << n << " of " << mZones.size() << std::endl;
//      std::vector<Zone>::iterator i = (mZones.begin() + n);
//      if(i not_eq mZones.end()) {
//         std::cout << i->thetaMin << std::endl;
//      } // if
//      else {
//         std::cerr << "No such zone" << std::endl;
//      } // if
      if(n < mZones.size()) {
         mZones.erase(mZones.begin() + n);
      } // if
   }
   
   void Acceptance::SetTheta(double min, double max, int n) {
//      try {
         mZones.at(n).thetaMin = FixTopologyTheta(min);
         mZones.at(n).thetaMax = FixTopologyTheta(max);
//      } // try
//      catch(std::exception& e) {
//         std::cerr << "Exception in Smear::Acceptance::SetTheta() - " <<
//         e.what() << std::endl;
//      } // catch
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
         case kPz: // Note deliberate fall-through
         case kPt:
            break; // Do nothing
      } // switch
   }
   
   void Acceptance::SetPt(double min, double max, int n) {
      AddCustomAcceptance("P*sin(theta)",min,max,n);
   }
   
   void Acceptance::SetPz(double min, double max, int n) {			
      AddCustomAcceptance("P*cos(theta)",min,max,n);
   }
   
   void Acceptance::SetGenre(int n) {
      if(n > 0 and n < 3) {
         Genre = n;
      } // if...
      else {
         Genre = 0;
      } // else
   }
   /*
   void Acceptance::SetGenre(TString genre) {
      if(genre.Contains("EM")) {
         Genre = 1;
      }
      else if(genre.Contains("Had") or genre.Contains("had")) {
         Genre = 2;
      }
      else {
         Genre = 0;
      }
   }
   */
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
      Particles.push_back(n);
   }
   
   void Acceptance::ClearParticles() {
      Particles.clear();
   }
   
   void Acceptance::RemoveParticle(int n) {
      for(unsigned i=0; i<Particles.size(); i++) {
         if(Particles.at(i)==n) {
            Particles.erase(Particles.begin()+i);
         } // if
         break;
      } // for
   }
   
   bool Acceptance::IsCustomAccepted(CustomCut& C, const Particle& prt) const {
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
   
   bool Acceptance::Zone::IsCustomAccepted(Acceptance::CustomCut C, const Particle& prt) const {
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
   
   Bool_t Acceptance::Zone::Contains(Particle& prt) const {
      bool accept(true);
      
      const double theta = FixTopologyTheta(prt.Theta());
      const double phi = FixTopologyPhi(prt.Phi());
      
      if(theta < thetaMin or theta > thetaMax) {
         accept = false;
      } // if...
      else if(phi < phiMin or phi > phiMax) {
         accept = false;
      } // else if...
      else if(prt.GetE() < EMin or prt.GetE() > EMax) {
         accept = false;
      } // else if...
      else if(prt.P() < PMin or prt.P() > PMax) {
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

   bool Acceptance::Is(const Particle& prt) {
      bool b = false;
      bool intheta;
      bool inphi;
      bool inE;
      bool inp; 
      bool inCustom=true;
//      if(Genre==2)std::cout<<"In HCal"<<std::endl;
//      if(Genre==1)std::cout<<"In ECal"<<std::endl;
//      bool isHadron = prt.KS==1 and abs(prt.id)>110 and Genre==2;
      
      if(PGenre(prt)==0 or (Genre not_eq 0 and PGenre(prt) not_eq Genre)) {
         return b;
      } // if
      
//      if(isHadron)std::cout<<"Hadron passed genre"<<std::endl;
      
      // Loop through mZones to see if accepted
      //#if 0
      //      std::cout<<"Looping over "<<mZones.size()<<" zones"<<std::endl;
      for(unsigned i(0); i < mZones.size(); i++) {
         double theta = prt.Theta();
         double phi = prt.Phi();
         intheta =
            theta >= mZones.at(i).thetaMin and theta <= mZones.at(i).thetaMax;
         inphi   =
            phi >= mZones.at(i).phiMin and phi <= mZones.at(i).phiMax;
         inE     =
            prt.GetE() >= mZones.at(i).EMin and prt.GetE() <= mZones.at(i).EMax;
         inp     =
            prt.P() >= mZones.at(i).PMin and prt.P() <= mZones.at(i).PMax;
//         std::cout<<"\tLooping through "<<mZones.at(i).CustomCuts.size()<<" custom cuts"<<std::endl;
         // Loop through custom cuts in zone, if there are any
         if(mZones.at(i).CustomCuts.size() not_eq 0) {
            for(unsigned j=0; j<mZones.at(i).CustomCuts.size(); j++) {
               inCustom = IsCustomAccepted(mZones.at(i).CustomCuts.at(j),prt);
            } // for
         } // if
           //         std::cout<<"Acceptances: "<<intheta<<" "<<inphi<<" "<<inE<<" "<<inp<<" "<<inCustom<<std::endl;
         if(intheta and inphi and inE and inp and inCustom) {
            b = true;
         } // if
           //         if(mZones.at(i).Contains(prt)) {
           //            b = true;
           //            break;
           //         } // if
      } // for
        //#endif
#if 0
      // Alternative implementation
      unsigned i(0);
      do {
         b = mZones.at(i++).Contains(prt);
      } while(not b and i < mZones.size());
#endif
//      std::cout<<"Looping through "<<Particles.size()<<" exclusive particles"<<std::endl;
      // Loop through exclusive particle list
      if(b and Particles.size() not_eq 0) {
         b = false;
         for(unsigned i=0; i<Particles.size(); i++) {
            if(prt.Id()==Particles.at(i)) {
               b = true;
            } // if
         } // for
      } // if
        //      std::cout<<"Accepted? "<<b<<std::endl;
      return b;
   }
   
   Acceptance::CustomCut::CustomCut()
   : dim(0)
   , Kin1(kP)
   , Kin2(kTheta) { }
   
   Acceptance::Zone::Zone(double thMin, double thMax,
                          double phMin, double phMax,
                          double eMin, double eMax,
                          double pMin, double pMax)
   : thetaMin(thMin)
   , thetaMax(thMax)
   , phiMin(phMin)
   , phiMax(phMax)
   , EMin(eMin)
   , EMax(eMax)
   , PMin(pMin)
   , PMax(pMax) { }
}
