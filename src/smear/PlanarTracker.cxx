/**
 PlanarTracker.cxx

 \file
 Implementation of class PlanarTracker.

 \author Thomas Burton 
 \date 6/1/12
 \copyright 2012 BNL. All rights reserved.
*/

#include "eicsmear/smear/PlanarTracker.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <list>

#include <TMath.h>

namespace Smear {

   PlanarTracker::PlanarTracker()
   : Tracker(2., 0.03, 0.001)
   , mNPlanes(25.)
   , mInnerRadius(0.1)
   , mOuterRadius(1.)
   , mZMin(-1.5)
   , mZMax(1.5) {
   }
   
   PlanarTracker::PlanarTracker(double inner, double outer,
                                double zmin, double zmax,
                                double Bb, double nrl,
                                double sigmaRPhi, double N)
   : Tracker(Bb, nrl, sigmaRPhi)
   , mNPlanes(N)
   , mInnerRadius(inner)
   , mOuterRadius(outer)
   // Require zmin < zmax
   , mZMin(std::min(zmin, zmax))
   , mZMax(std::max(zmin, zmax)) {
   }
      
   PlanarTracker::~PlanarTracker() {
   }

   void PlanarTracker::Print(Option_t*) const {
      std::cout << ClassName() << " with:" << std::endl <<
      "\tinner radius " << mInnerRadius << " m\n" <<
      "\touter radius " << mOuterRadius << " m\n" <<
      "\tlength " << mZMin << " to " << mZMax <<
      " (= " << mZMax - mZMin << ") m\n" <<
      "\tmagnetic field " << mMagField << " Tesla\n" <<
      "\t" << mNRadLengths << " radiation lengths\n" <<
      "\tpoint resolution " << mSigmaRPhi * 1.e6 << " microns\n" <<
      "\t" << mNPlanes << " planes" << std::endl;
   }

   double PlanarTracker::GetThetaMin() {
      if(mZMax > 0.) {
         return atan2(mInnerRadius, mZMax);
      } // if
      else {
         return atan2(mOuterRadius, mZMax);
      } // else
   }

   double PlanarTracker::GetThetaMax() {
      if(mZMin > 0.) {
         return atan2(mOuterRadius, mZMin);
      } // if
      else {
         return atan2(mInnerRadius, mZMin);
      } // else
   }

   //
   // The following need implementing:
   //
   
   double PlanarTracker::L(const erhic::VirtualParticle& p) const {
      return 0.;
   }

   double PlanarTracker::LPrime(const erhic::VirtualParticle& p) const {
      return 0.;
   }

   int PlanarTracker::NPoints(const erhic::VirtualParticle&) const {
      return mNPlanes;
   }

   bool PlanarTracker::Accepts(const erhic::VirtualParticle& p) const {
      return true;
   }

} // namespace Smear
