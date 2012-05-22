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
#include <cmath>
#include <limits>
#include <list>

#include <TMath.h>

namespace {
   // Functor for testing non-existent intersection points,
   struct NoIntersection {
      bool operator()(const TVector3& v) const {
         return TMath::IsNaN(v.z());
      }
   };
} // namespace

namespace Smear {

   Tracker::Tracker()
   : mMagField(2.)
   , mNRadLengths(0.03)
   , mSigmaRPhi(0.001)
   , mNFitPoints(25.)
   , mInnerRadius(0.1)
   , mOuterRadius(1.)
   , mZMin(-1.5)
   , mZMax(1.5) {
   }

   Tracker::Tracker(double inner, double outer, double L, double Bb,
                    double nrl, double sigmaRPhi, double N)
   : mMagField(Bb)
   , mNRadLengths(nrl)
   , mSigmaRPhi(sigmaRPhi)
   , mNFitPoints(N)
   , mInnerRadius(inner)
   , mOuterRadius(outer)
   , mZMin(-L / 2.)
   , mZMax(L / 2.) {
   }

   Tracker::Tracker(double inner, double outer, double zmin, double zmax,
                    double Bb, double nrl, double sigmaRPhi, double N)
   : mMagField(Bb)
   , mNRadLengths(nrl)
   , mSigmaRPhi(sigmaRPhi)
   , mNFitPoints(N)
   , mInnerRadius(inner)
   , mOuterRadius(outer)
     // Require zmin < zmax
   , mZMin(std::min(zmin, zmax))
   , mZMax(std::max(zmin, zmax)) {
   }

   Tracker::~Tracker() {
   }

   double Tracker::MultipleScatteringContribution(
                                       const erhic::VirtualParticle& p) const {
      // Technically should be a factor of particle charge in the numerator
      // but this is effectively always one.
      return 0.0136 / 0.3 * p.GetP() * sqrt(mNRadLengths) / L(p) /
             p.Get4Vector().Beta() / mMagField;
   }

   double Tracker::IntrinsicContribution(
                      const erhic::VirtualParticle& p) const {
      return sqrt(720.) / 0.3 * pow(p.GetP(), 2.) * mSigmaRPhi /
             mMagField / pow(LPrime(p), 2.) / sqrt(mNFitPoints + 4.);
   }

   double Tracker::Resolution(const erhic::VirtualParticle& p) const {
      // Don't compute resolution for very small path length
      // as L in the denominator causes it to blow up.
      if(not Accepts(p)) {
         return 0.;
      } // if
      return sqrt(pow(MultipleScatteringContribution(p), 2.) +
                  pow(IntrinsicContribution(p), 2.));
   }

   TVector3 Tracker::ComputeIntersectionWithRadius(
      const erhic::VirtualParticle& p, double radius) const {
      // Compute the longitudinal position at which the intersection occurs.
      // Adjust for any z offset in the vertex of the particle.
      const double z = radius / tan(p.GetTheta()) + p.GetVertex().z();
      TVector3 intersection(0., 0., std::numeric_limits<double>::quiet_NaN());
      if(z > mZMin and z < mZMax) {
         // Don't care about phi angle, so set x and y arbitrarily.
         intersection.SetXYZ(radius, 0., z);
      } // if
      return intersection;
   }

   TVector3 Tracker::ComputeIntersectionWithPlane(
      const erhic::VirtualParticle& p, double z) const {
      // Compute the radius of intersection, adusting for any z offset
      // in the particle vertex.
      const double r = (z - p.GetVertex().z()) * tan(p.GetTheta());
      TVector3 intersection(0., 0., std::numeric_limits<double>::quiet_NaN());
      if(r > mInnerRadius and r < mOuterRadius) {
         // Don't care about phi angle, so set x and y arbitrarily.
         intersection.SetXYZ(r, 0., z);
      } // if
      return intersection;
   }

   TVector3 Tracker::ComputePath(const erhic::VirtualParticle& p) const {
      // There are 4 valid ways for the particle to intersect the cylinder.
      // It can intersect:
      // 1) nothing
      // 2) both faces
      // 3) both radii
      // 4) one face and one radius
      // Finding anything else indicates an error - if it
      // enters the (finite) volume, it must exit it.
      // So, first compute all 4 possible intersection points
      std::list<TVector3> xyz;
      xyz.push_back(ComputeIntersectionWithRadius(p, mInnerRadius));
      xyz.push_back(ComputeIntersectionWithRadius(p, mOuterRadius));
      xyz.push_back(ComputeIntersectionWithPlane(p, mZMin));
      xyz.push_back(ComputeIntersectionWithPlane(p, mZMax));
      // Remove those where there was no intersection
      xyz.erase(std::remove_if(xyz.begin(), xyz.end(), NoIntersection()),
                xyz.end());
      // We should have either exactly 2 intersection points, or none, in
      // which case the particle missed the detector and we return zero.
      // Anything else indicates an error, so return zero path length for
      // those also.
      TVector3 path(0., 0., 0.);
      if(2 == xyz.size()) {
         // Arrange the points so that the path vector is going outward
         // from the origin
         if(xyz.front().Mag() > xyz.back().Mag()) {
            path = xyz.front() - xyz.back();
         } // if
         else {
            path = xyz.back() - xyz.front();
         } // else 
      } // if
      // 
      return path;
   }

   double Tracker::L(const erhic::VirtualParticle& p) const {
      return ComputePath(p).Mag();
   }

   double Tracker::LPrime(const erhic::VirtualParticle& p) const {
      return ComputePath(p).Perp();
   }

   bool Tracker::Accepts(const erhic::VirtualParticle& p) const {
      // Require the transverse path length to exceed half of the
      // radial width per fit point, otherwise the detector essentially
      // doesn't "see" the particle.
      const double l = LPrime(p);
      return l > 0.5 * (mOuterRadius - mInnerRadius) / mNFitPoints;
   }

   void Tracker::Smear(const erhic::VirtualParticle& pIn,
                       ParticleMCS& pOut) {
      if(Accept.Is(pIn)) {
         double y = GetVariable(pIn, kP);
         // Randomly generate a smeared value from the resolution
         // and set it in the smeared particle.
         SetVariable(pOut, Distribution.Generate(y, Resolution(pIn)), kP);
         // Ensure E, p are positive definite
         HandleBogusValues(pOut, kP);
      } //if
   }

   void Tracker::Print(Option_t*) const {
      std::cout << "Tracker with:" << std::endl <<
      "\tinner radius " << mInnerRadius << " m\n" <<
      "\touter radius " << mOuterRadius << " m\n" <<
      "\tlength " << mZMin << " to " << mZMax <<
      " (= " << mZMax - mZMin << ") m\n" <<
      "\tmagnetic field " << mMagField << " Tesla\n" <<
      "\t" << mNRadLengths << " radiation lengths\n" <<
      "\tpoint resolution " << mSigmaRPhi * 1.e6 << " microns\n" <<
      "\t" << mNFitPoints << " fit points" << std::endl;
   }

   double Tracker::GetThetaMin() {
      if(mZMax > 0.) {
         return atan2(mInnerRadius, mZMax);
      } // if
      else {
         return atan2(mOuterRadius, mZMax);
      } // else
   }

   double Tracker::GetThetaMax() {
      if(mZMin > 0.) {
         return atan2(mOuterRadius, mZMin);
      } // if
      else {
         return atan2(mInnerRadius, mZMin);
      } // else
   }
} // namespace Smear
