/**
 \file
 Implementation of class Smear::Tracker.
 
 \author    Will Foreman
 \date      2011-08-19
 \copyright 2011 Brookhaven National Lab
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

}  // anonymous namespace

namespace Smear {

Tracker::Tracker(double magneticField, double nRadiationLengths,
                 double resolution)
: mFactor(720)  // Assume no vertex constraint.
, mMagField(magneticField)
, mNRadLengths(nRadiationLengths)
, mSigmaRPhi(resolution) {
}

Tracker::~Tracker() {
}

double Tracker::MultipleScatteringContribution(
    const erhic::VirtualParticle& p) const {
  // Technically should be a factor of particle charge in the numerator
  // but this is effectively always one.
  double val = 0.016 / 0.3 * p.GetP() * sqrt(mNRadLengths) / LPrime(p) /
  p.Get4Vector().Beta() / mMagField;
  if (TMath::IsNaN(val)) {
    std::cerr << "MS nan!" << std::endl;
  }  // if
  return val;
}

double Tracker::IntrinsicContribution(
    const erhic::VirtualParticle& p) const {
  // The factor
  //    sqrt(720 * N^3 / ((N-1)(N+1)(N+2)(N+3)))
  // is a more exact version of the factor
  //    sqrt(720 / (N+5))
  // The constant factor in the sqrt depends on whether there is
  // a vertex constraint assumed.
  double val = sqrt(mFactor * pow(NPoints(p), 3.)) / 0.3 * pow(p.GetP(), 2.)
  * mSigmaRPhi / mMagField / pow(LPrime(p), 2.)
  / sqrt((NPoints(p)-1)*(NPoints(p)+1)*(NPoints(p)+2)*(NPoints(p)+3));
  if (TMath::IsNaN(val)) {
    std::cerr << "Intrinsic nan!" << std::endl;
  }  // if
  return val;
}

double Tracker::Resolution(const erhic::VirtualParticle& p) const {
  // Don't compute resolution for very small path length
  // as L in the denominator causes it to blow up.
  if (!Accepts(p)) {
    return 0.;
  }  // if
  return sqrt(pow(MultipleScatteringContribution(p), 2.) +
              pow(IntrinsicContribution(p), 2.));
}

void Tracker::Smear(const erhic::VirtualParticle& pIn,
                    ParticleMCS& pOut) {
  if (Accepts(pIn) && Accept.Is(pIn)) {
    double y = GetVariable(pIn, kP);
    // Randomly generate a smeared value from the resolution
    // and set it in the smeared particle.
    SetVariable(pOut, Distribution.Generate(y, Resolution(pIn)), kP);
    // Ensure E, p are positive definite
    HandleBogusValues(pOut, kP);
    if (pOut.GetP() < 0.) {
      std::cerr << "p " << pOut.GetP() << std::endl;
    }  // if
    if (TMath::IsNaN(pOut.GetP())) {
      std::cerr << "p nan" << std::endl;
    }  // if
  }  // if
}

void Tracker::SetVertexConstraint(bool constrain) {
  if (constrain) {
    mFactor = 320;
  } else {
    mFactor = 720;
  }  // if
}

}  // namespace Smear
