/**
 \file
 Declaration of various Smear namespace functions.
 
 \author    Michael Savastio
 \date      2011-08-19
 \copyright 2011 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_SMEAR_SMEAR_H_
#define INCLUDE_EICSMEAR_SMEAR_SMEAR_H_

#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>

#include <TF1.h>
#include <TF2.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TVector2.h>

#include "eicsmear/erhic/Kinematics.h"
#include "eicsmear/erhic/Particle.h"
#include "eicsmear/erhic/ParticleIdentifier.h"
#include "eicsmear/erhic/VirtualParticle.h"
#include "eicsmear/smear/ParticleMCS.h"

namespace Smear {

/**
 Enumerator listing particle wise kinematic variables.
 Naming is self explanitory.
 These will be used when specifying arguments and outputs
 for device parametrizations.
 */
enum KinType {
  kE, kP, kTheta, kPhi, kPz, kPt, kInvalidKinType
};

/** Classes of particles */
enum EGenre {
  kAll = 0, kElectromagnetic = 1, kHadronic = 2
};

/** Particle charged **/
enum ECharge {
  kNeutral, kCharged, kAllCharges
};

/**
 Determine particle "genre".
 
 This refers to whether the particle in the argument is
 "electromagnetic" or "hadronic" from the perspective of calorimetry.
 
 Returns one of the following (see enum EGenre):
 kElectromagnetic: photon or electron/positron.
 kHadronic:        stable hadron.
 kAll:             neither of the above.
 */
inline int PGenre(const erhic::VirtualParticle& prt) {
  int genre(kAll);
  const int id = abs(prt.Id());  // Sign doesn't matter
  if (1 == prt.GetStatus()) {  // Only check stable particles.
    if (id == 11 || id == 22) {
      genre = kElectromagnetic;
    } else if (id >110) {
      genre = kHadronic;
    }  // if
  }  // if
  return genre;
}

/**
 Fix a polar angle so that it lies within [0,pi].
 TODO Nothing Smear-specific here - move to general functions file.
 */
inline double FixTheta(double theta) {
  while (theta < 0. || theta > TMath::Pi()) {
    if (theta < 0.) {
      theta = -theta;
    }  // if
    if (theta > TMath::Pi()) {
      theta = TMath::TwoPi() - theta;
    }  // if
  }
  return theta;
}

/**
 Fix an azimuthal angle so that it lies within [0,2*pi).
 TODO Nothing Smear-specific here - move to general functions file.
 */
inline double FixPhi(double phi) {
  return TVector2::Phi_0_2pi(phi);
}

/**
 Returns the kinematic variable associated with kin from the input particle.
 */
inline double GetVariable(const erhic::VirtualParticle& prt, KinType kin) {
  double z(0.);
  switch (kin) {
    case kE:
      z = prt.GetE(); break;
    case kP:
      z = prt.GetP(); break;
    case kTheta:
      z = prt.GetTheta(); break;
    case kPhi:
      z = prt.GetPhi(); break;
    case kPz:
      z = prt.GetPz(); break;
    case kPt:
      z = prt.GetPt(); break;
    default:
      break;
  }  // switch
  return z;
}

/**
 Stores z in the ParticleS.K where K is the kinematic variable associated with kin.
 */
inline void SetVariable(ParticleMCS &prt, double z, KinType kin) {
  switch (kin) {
    case kE:
      prt.SetE(z); break;
    case kP:
      prt.SetP(z); break;
    case kTheta:
      prt.SetTheta(z); break;
    case kPhi:
      prt.SetPhi(z); break;
    case kPz:
      prt.SetPz(z); break;
    case kPt:
      prt.SetPt(z); break;
    default:
      break;
  }  // switch
}

/**
 This dictates how the namespace deals with positive definite variables
 which have been smeared to negative values.
 */
inline void HandleBogusValues(ParticleMCS &prt, KinType kin) {
  double fault(0.);
  if (kE == kin && prt.GetE() < 0.) {
    prt.SetE(fault);
  } else if (kP == kin && prt.GetP() < 0.) {
    prt.SetP(fault);
  } else if (kPt == kin && prt.GetPt() < 0.) {
    prt.SetPt(fault);
  } else if (kPz == kin && prt.GetPz() < 0.) {
    prt.SetPz(fault);
  }  // if
}

inline void HandleBogusValues(ParticleMCS& prt) {
  double fault = NAN;  // -999.;
  if (prt.GetE() < 0.) {
    prt.SetE(fault);
  }  // if
  if (prt.GetP() < 0.) {
    prt.SetP(fault);
  }  // if
  if (prt.GetPt() < 0.) {
    prt.SetPt(fault);
  }  // if
  if (prt.GetPz() < 0.) {
    prt.SetPz(fault);
  }  // if
}

inline bool IsCoreType(KinType kin) {
  if (kin == kE || kin == kP || kin == kTheta || kin == kPhi) return true;
  return false;
}

int ParseInputFunction(TString &s, KinType &kin1, KinType &kin2);

}  // namespace Smear

#endif  // INCLUDE_EICSMEAR_SMEAR_SMEAR_H_
