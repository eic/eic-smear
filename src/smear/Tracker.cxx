/**
 Tracker.cxx

 \file
 Implementation of class Tracker.

 \author Thomas Burton 
 \date 5/8/12
 \copyright 2012 BNL. All rights reserved.
*/

#include "eicsmear/smear/Tracker.h"

#include <algorithm>

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
//      Accept.AddZone(Acceptance::Zone(GetThetaMin(),
//                                      TMath::Pi() - GetThetaMin()));
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
//      Accept.AddZone(Acceptance::Zone(GetThetaMin(),
//                                      TMath::Pi() - GetThetaMin()));
   }

   Tracker::~Tracker() {
   }

   double Tracker::MultipleScatteringContribution(
                                       const erhic::VirtualParticle& p) const {
      // Technically should be a factor of particle charge in the numerator
      // but this is effectively always one.
      double c1 = 0.0136 / 0.3; // * z
      double beta = p.Get4Vector().Beta();
      // gamma is the angle to the vertical, theta is angle to the horizontal.
      return c1 * p.GetP() * sqrt(mNRadLengths) / (L(p) * beta * mMagField);
   }

   double Tracker::IntrinsicContribution(
                      const erhic::VirtualParticle& p) const {
      double c1 = sqrt(720.) / 0.3;
      double over = c1 * p.GetP() * p.GetP() * mSigmaRPhi;
      return over / (mMagField * pow(LPrime(p), 2.) * sqrt(mNFitPoints + 4.));
   }

   double Tracker::Evaluate(const erhic::VirtualParticle& p) const {
      return sqrt(pow(MultipleScatteringContribution(p), 2.) +
                  pow(IntrinsicContribution(p), 2.));
   }

   double Tracker::L(const erhic::VirtualParticle& p) const {
      double l(0.);
      if(sin(p.GetTheta()) > 0.) {
         l = LPrime(p) / sin(p.GetTheta());
      } // if
      return l;
   }

   double Tracker::LPrime(const erhic::VirtualParticle& p) const {
      // Compute the transverse track length.
      // Use abs(tan(theta)) as we don't distinguish forward and backward angles
      double l = mLength / 2. * fabs(tan(p.GetTheta())) - mInnerRadius;
      // Make sure to cap it at both min and max values
      // so it lies in the range [0, outer radius - inner radius]
      return std::min(std::max(l, 0.), mOuterRadius - mInnerRadius);
   }

   void Tracker::Smear(const erhic::VirtualParticle& prt,
                          ParticleMCS& prtOut) {
      bool accept = prt.GetTheta() > GetThetaMin() and 
                    prt.GetTheta() < TMath::Pi() - GetThetaMin();
      if(Accept.GetNZones() > 0) {
         accept = accept and Accept.Is(prt);
      } // if
      if(accept) {
         double y = GetVariable(prt, kP);
         y = Distribution.Generate(y, Evaluate(prt));
         SetVariable(prtOut, y, kP);
         //make sure E, p are positive definite
         HandleBogusValues(prtOut, kP);
      } //if
   }
   void Tracker::Print(Option_t*) const {
      std::cout << "Tracker with:" << std::endl <<
      "\tinner radius " << mInnerRadius << " m\n" <<
      "\touter radius " << mOuterRadius << " m\n" <<
      "\tlenth " << mLength << " m\n" <<
      "\tmagnetic field " << mMagField << " Tesla\n" <<
      "\t" << mNRadLengths << " radiation lengths\n" <<
      "\tpoint resolution " << mSigmaRPhi * 1.e6 << " microns\n" <<
      "\t" << mNFitPoints << " fit points" << std::endl;
   }
} // namespace Smear
