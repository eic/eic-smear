/**
 \file
 Implementation of class Smear::PlanarTracker.
 
 \author    Will Foreman
 \date      2011-08-19
 \copyright 2011 Brookhaven National Lab
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

void PlanarTracker::Print(Option_t* /* option */) const {
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

double PlanarTracker::GetThetaMin() const {
  if (mZMax > 0.) {
    return atan2(mInnerRadius, mZMax);
  } else {
    return atan2(mOuterRadius, mZMax);
  }  // if
}

double PlanarTracker::GetThetaMax() const {
  if (mZMin > 0.) {
    return atan2(mOuterRadius, mZMin);
  } else {
    return atan2(mInnerRadius, mZMin);
  }  // if
}

TVector3 PlanarTracker::ComputeIntersectionWithRadius(
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

TVector3 PlanarTracker::ComputeIntersectionWithPlane(
      const erhic::VirtualParticle& p, double z) const {
  // Compute the radius of intersection, adusting for any z ooffset
  // in the particle vertex.
  const double r = (z - p.GetVertex().z()) * tan(p.GetTheta());
  TVector3 intersection(0., 0., std::numeric_limits<double>::quiet_NaN());
  if (r > mInnerRadius && r < mOuterRadius) {
    // Don't care about phi angle, so set x and y arbitrarily.
    intersection.SetXYZ(r, 0., z);
  }  // if
  return intersection;
}

double PlanarTracker::L(const erhic::VirtualParticle& p) const {
  float r1 = 0., r2 = 0., pi = 3.1415926535;
  double zPlane1 = 0., zPlane2 = 0., Length = 0.;
  // Find radial/z points of first and last plane intersection.
  // (While there are certainly eaiser ways to carry out this
  // calculation, I used a manual check at each z-plane since
  // it is easier to generalize to arbitrarily spaced planes.)
  for (int i = 0; i < mNPlanes; i++) {
    zPlane1 = mZMin + i * (mZMax - mZMin) / (mNPlanes - 1);
    r1 = fabs(tan(p.GetTheta()) * zPlane1);
    if ((r1 > mInnerRadius) && (r1 < mOuterRadius)) {
      if (((p.GetTheta() < pi / 2.) && (zPlane1 > 0.)) ||
         ((p.GetTheta() > pi / 2.) && (zPlane1 < 0.))) {
        zPlane2 = zPlane1 + (NPoints(p) - 1) * (mZMax - mZMin) / (mNPlanes - 1);
        r2 = fabs(tan(p.GetTheta()) * zPlane2);
        break;
      }  // if
    }  // if
  }  // for
  Length = sqrt((zPlane2 - zPlane1) * (zPlane2 - zPlane1) +
                (r2 - r1) * (r2 - r1));
  return Length;
}

double PlanarTracker::LPrime(const erhic::VirtualParticle& p) const {
  float r1 = 0., r2 = 0., pi = 3.1415926535;
  double zPlane1 = 0., zPlane2 = 0.;
  // Find radial points of first/ and last plane intersection.
  for (int i = 0; i < mNPlanes; i++) {
    zPlane1 = mZMin + i * (mZMax - mZMin) / (mNPlanes-1);
    r1 = fabs(tan(p.GetTheta()) * zPlane1);
    if ((r1 > mInnerRadius) && (r1 < mOuterRadius)) {
      if (((p.GetTheta() < pi / 2.) && (zPlane1 > 0.)) ||
         ((p.GetTheta() > pi / 2.) && (zPlane1 < 0.))) {
        zPlane2 = zPlane1 + (NPoints(p) - 1) * (mZMax - mZMin) / (mNPlanes - 1);
        r2 = fabs(tan(p.GetTheta()) * zPlane2);
        break;
      }  // if
    }  // if
  }  // for
  return fabs(r2 - r1);
}

int PlanarTracker::NPoints(const erhic::VirtualParticle& p) const {
  double zPlane = 0., pi = 3.1415926535;
  float r = 0.;
  int i, n = 0;
  for (i = 0; i < mNPlanes; i++) {
    zPlane = mZMin + i * (mZMax - mZMin) / (mNPlanes - 1);
    // radial intersect of particle at plane position zPlane
    r = fabs(tan(p.GetTheta()) * zPlane);
    // At each plane, check if particle passes through
    if ((r > mInnerRadius) && (r < mOuterRadius)) {
      if ((p.GetTheta() < pi / 2.) && (zPlane > 0)) n++;
      if ((p.GetTheta() > pi / 2.) && (zPlane < 0)) n++;
    }  // if
  }  // for
  return n;
}

bool PlanarTracker::Accepts(const erhic::VirtualParticle& p) const {
  if (NPoints(p) >= 3) {
    return true;
  } else {
    return false;
  }  // if
}

}  // namespace Smear
