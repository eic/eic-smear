/**
 Tracker.cxx

 \file
 Implementation of class Tracker.

 \author Thomas Burton 
 \date 5/8/12
 \copyright 2012 BNL. All rights reserved.
*/

#include "eicsmear/smear/Tracker.h"

#include <TMath.h>

namespace Smear {

   Tracker::Tracker()
   : mMagField(2.)
   , mNRadLengths(0.03)
   , mSigmaRPhi(0.001)
   , mNFitPoints(25.)
   , mInnerRadius(0.1)
   , mOuterRadius(1.)
   , mLength(6.) {
      Accept.AddZone(Acceptance::Zone(GetThetaMin(),
                                      TMath::Pi() - GetThetaMin()));
   }

   Tracker::Tracker(double inner, double outer, double L, double Bb,
                    double nrl, double sigmaRPhi, double N)
   : mMagField(Bb)
   , mNRadLengths(nrl)
   , mSigmaRPhi(sigmaRPhi)
   , mNFitPoints(N)
   , mInnerRadius(inner)
   , mOuterRadius(outer)
   , mLength(L) {
      Accept.AddZone(Acceptance::Zone(GetThetaMin(),
                                      TMath::Pi() - GetThetaMin()));
   }

   Tracker::~Tracker() {
   }

   double Tracker::MultipleScatteringContribution(
                                       const erhic::VirtualParticle& p) {
      double theta = p.GetTheta();
      double c1 = 0.0136 / 0.3;
      double beta = p.Get4Vector().Beta();
      double cos2gamma = pow(cos(TMath::PiOver2() - theta), 2.);
      return c1 * p.GetP() * sqrt(mNRadLengths) /
             (L(p) * beta * mMagField * cos2gamma);
   }

   double Tracker::IntrinsicContribution(const erhic::VirtualParticle& p) {
      double c1 = sqrt(720.) / 0.3;
      double over = c1 * p.GetP() * p.GetP() * mSigmaRPhi;
      return over / (mMagField * pow(LPrime(p), 2.) * sqrt(mNFitPoints + 4.));
   }

   double Tracker::EvaluateRes(const erhic::VirtualParticle& p) {
         return MultipleScatteringContribution(p) + IntrinsicContribution(p);
   }

   double Tracker::L(const erhic::VirtualParticle& p) {
      if(p.GetTheta() >= ThetaCrit() and
         p.GetTheta() < TMath::Pi() - ThetaCrit()) {
         return LPrime(p) / sin(p.GetTheta());
      } // if
      else {
         return sqrt(pow(LPrime(p), 2.) +
                     pow((mLength / 2.) - (mInnerRadius / tanTheta(p)), .2));
      } // else
   }

   double Tracker::LPrime(const erhic::VirtualParticle& p) {
      if(p.GetTheta() >= ThetaCrit() and
         p.GetTheta() < TMath::Pi() - ThetaCrit()) {
         return mOuterRadius - mInnerRadius;
      } // if
      else {
         return (mLength / 2.) * tanTheta(p) - mInnerRadius;
      } // else
   }

   void Tracker::Smear(const erhic::VirtualParticle& prt,
                          ParticleMCS& prtOut) {
      if(Accept.Is(prt)) {
         double y = GetVariable(prt, kP);
         y = Distribution.Generate(y, EvaluateRes(prt));
         SetVariable(prtOut, y, kP);
         //make sure E, p are positive definite
         HandleBogusValues(prtOut, kP);
      } //if
   }
} // namespace Smear
