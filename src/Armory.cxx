/**
 Armory.cxx

 \file
 Implementation of class Armory.

 \author Thomas Burton 
 \date 4/25/12
 \copyright 2012 BNL. All rights reserved.
*/

#include "Armory.h"

namespace Smear {

   Tracker::Tracker()
   : B(2.)
   , Nrl(0.03)
   , SigmaRPhi(0.001)
   , n(25.)
   , R1(0.1)
   , R2(1.)
   , l(6.)
   , fault(-999.)
   , bMoreCrit(false)
   , P(NULL)
   , OutKin(kP) {
      FixAcceptance();
   }
   
   Tracker::Tracker(double inner, double outer, double L, double Bb,
                    double nrl, double sigmaRPhi, double N) {
      SetVals(inner, outer, L, Bb, nrl, sigmaRPhi, N);
   }

   Tracker::Tracker(double inner, double outer, double L) {
      SetDimensions(inner, outer, L);
   }

   Tracker::~Tracker() {
   }

   double Tracker::MultipleScatteringContribution() {
      if(not P) {
         return fault;
      } // if
      double c1 = 0.0136 / 0.3;
      return c1 * P->p * sqrt(Nrl) / (L() * beta() * B * cosSquaredgamma());
   }

   double Tracker::IntrinsicContribution() {
      if(not P) {
         return fault;
      } // if
      double c1 = sqrt(720.) / 0.3;
      double over = c1 * P->p * P->p * SigmaRPhi;
      return over / (B * pow(LPrime(), 2.) * sqrt(n + 4.));
   }

   double Tracker::EvaluateRes() {
      if(P) {
         return MultipleScatteringContribution() + IntrinsicContribution();
      } // if
      else {
         return fault;
      } // else
   }

   void Tracker::SetRadii(double inner, double outer) {
      if(inner < outer) {
         R1 = inner;
         R2 = outer;
         FixAcceptance();
      } // if
      else {
         std::cerr << "ERROR! Outer radius must exceed inner." << std::endl;
      } // else
   }

   double Tracker::L() {
      if(not P) {
         return fault;
      } // if
      if(bMoreCrit) {
         return LPrime() / sin(P->theta);
      } // if
      else {
         return sqrt(pow(LPrime(), 2.) + pow((l / 2.) - (R1 / tanTheta()), .2));
      } // else
   }

   double Tracker::LPrime() {
      if(not P) {
         return fault;
      } // if
      if(bMoreCrit) {
         return R2 - R1;
      } // if
      else {
         return (l / 2.) * tanTheta() - R1;
      } // else
   }

   void Tracker::SetParticle(const Particle &prt) {
      P = &prt;
      if(P->theta >= ThetaCrit() and P->theta < pi-ThetaCrit()) {
         bMoreCrit = true;
      } // if
      else {
         bMoreCrit = false;
      } // else
   }

   void Tracker::SetVals(double inner, double outer, double L, double Bb,
                         double nrl, double sigmaRPhi, double N) {
      SetDimensions(inner, outer, L);
      SetMagneticField(Bb);
      SetNumberOfRadiationLengths(nrl);
      SetPositionResolution(sigmaRPhi);
      SetNumberOfMeasurements(N);
   }

   void Tracker::DevSmear(const Particle& prt, ParticleS& prtOut) {
      Double32_t y(0.);
      if(Accept.Is(prt)) {
         SetParticle(prt);
         y = SwitchKinGetFromParticle(prt, OutKin);
         y = Distribution.Generate(y, EvaluateRes());
         SwitchKinStoreToParticle(prtOut, y, OutKin);
         //make sure angular coordinates live in S^2
         if(OutKin == kTheta) {
            prtOut.theta = FixTopologyTheta(prtOut.theta);
         } // if
         if(OutKin == kPhi) {
            prtOut.phi = FixTopologyPhi(prtOut.phi);
         } // if
         //make sure E, p are positive definite
         HandleBogusValues(prtOut, OutKin);
      } //if
   }
} // namespace Smear
