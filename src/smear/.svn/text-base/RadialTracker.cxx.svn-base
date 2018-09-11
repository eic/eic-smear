/**
 \file
 Implementation of class Smear::RadialTracker.
 
 \author    Will Foreman
 \date      2011-08-19
 \copyright 2011 Brookhaven National Lab
 */

#include "eicsmear/smear/RadialTracker.h"

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

}  // anonymous namespace

namespace Smear {

RadialTracker::RadialTracker()
: Tracker(2., 0.03, 0.001)
, mNFitPoints(25.)
, mInnerRadius(0.1)
, mOuterRadius(1.)
, mZMin(-1.5)
, mZMax(1.5) {
}

RadialTracker::RadialTracker(double inner, double outer,
                             double zmin, double zmax,
                             double Bb, double nrl,
                             double sigmaRPhi, double N)
: Tracker(Bb, nrl, sigmaRPhi)
, mNFitPoints(N)
, mInnerRadius(inner)
, mOuterRadius(outer)
// Require zmin < zmax
, mZMin(std::min(zmin, zmax))
, mZMax(std::max(zmin, zmax)) {
}

RadialTracker::~RadialTracker() {
}

void RadialTracker::Print(Option_t* /* option */) const {
  std::cout << ClassName() << " with:" << std::endl <<
  "\tinner radius " << mInnerRadius << " m\n" <<
  "\touter radius " << mOuterRadius << " m\n" <<
  "\tlength " << mZMin << " to " << mZMax <<
  " (= " << mZMax - mZMin << ") m\n" <<
  "\tmagnetic field " << mMagField << " Tesla\n" <<
  "\t" << mNRadLengths << " radiation lengths\n" <<
  "\tpoint resolution " << mSigmaRPhi * 1.e6 << " microns\n" <<
  "\t" << mNFitPoints << " fit points" << std::endl;
}

TVector3 RadialTracker::ComputeIntersectionWithRadius(
    const erhic::VirtualParticle& p, double radius) const {
  // Compute the longitudinal position at which the intersection occurs.
  // Adjust for any z offset in the vertex of the particle.
  const double z = radius / tan(p.GetTheta()) + p.GetVertex().z();
  TVector3 intersection(0., 0., std::numeric_limits<double>::quiet_NaN());
  if (z > mZMin && z < mZMax) {
    // Don't care about phi angle, so set x and y arbitrarily.
    intersection.SetXYZ(radius, 0., z);
  }  // if
  return intersection;
}

TVector3 RadialTracker::ComputeIntersectionWithPlane(
    const erhic::VirtualParticle& p, double z) const {
  // Compute the radius of intersection, adusting for any z offset
  // in the particle vertex.
  const double r = (z - p.GetVertex().z()) * tan(p.GetTheta());
  TVector3 intersection(0., 0., std::numeric_limits<double>::quiet_NaN());
  if (r > mInnerRadius && r < mOuterRadius) {
    // Don't care about phi angle, so set x and y arbitrarily.
    intersection.SetXYZ(r, 0., z);
  }  // if
  return intersection;
}

TVector3 RadialTracker::ComputePath(const erhic::VirtualParticle& p) const {
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
  if (2 == xyz.size()) {
    // Arrange the points so that the path vector is going outward
    // from the origin
    if (xyz.front().Mag() > xyz.back().Mag()) {
      path = xyz.front() - xyz.back();
    } else {
      path = xyz.back() - xyz.front();
    }  // if
  }  // if
  return path;
}

double RadialTracker::L(const erhic::VirtualParticle& p) const {
  double Length = ComputePath(p).Mag();
  return Length;
}

double RadialTracker::LPrime(const erhic::VirtualParticle& p) const {
  return ComputePath(p).Perp();
}

int RadialTracker::NPoints(const erhic::VirtualParticle& p) const {
  int n;
  float n_float;
  if (LPrime(p) == (mOuterRadius - mInnerRadius)) {
    n = mNFitPoints;
  } else {
    n_float = mNFitPoints * ComputePath(p).Perp() /
              (mOuterRadius - mInnerRadius);
    n = floor(n_float + 0.5);
  }  // if
  return n;
}

bool RadialTracker::Accepts(const erhic::VirtualParticle& p) const {
  // Require the transverse path length to exceed half of the
  // radial width per fit point, otherwise the detector essentially
  // doesn't "see" the particle.
  if (NPoints(p) > 2) {
    return true;
  }  // if
  return false;
}

double RadialTracker::GetThetaMin() const {
  if (mZMax > 0.) {
    return atan2(mInnerRadius, mZMax);
  } else {
    return atan2(mOuterRadius, mZMax);
  }  // if
}

double RadialTracker::GetThetaMax() const {
  if (mZMin > 0.) {
    return atan2(mOuterRadius, mZMin);
  } else {
    return atan2(mInnerRadius, mZMin);
  }  // if
}

}  // namespace Smear
