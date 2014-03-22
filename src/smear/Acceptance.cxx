/**
 \file
 Implementation of class Smear::Acceptance.
 
 \author    Michael Savastio
 \date      2011-08-19
 \copyright 2011 Brookhaven National Lab
 */

#include "eicsmear/smear/Acceptance.h"

#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TString.h>

namespace Smear {

Acceptance::~Acceptance() {
}

Acceptance::Acceptance(int genre)
: mGenre(genre)
, mCharge(kAllCharges) {
}

void Acceptance::AddZone(const Zone& z) {
  mZones.push_back(z);
}

void Acceptance::SetGenre(int n) {
  if (n > 0 && n < 3) {
    mGenre = n;
  } else {
    mGenre = 0;
  }  // if
}

void Acceptance::SetCharge(ECharge charge) {
  mCharge = charge;
}

void Acceptance::AddParticle(int n) {
  mParticles.insert(n);
}

bool Acceptance::Is(const erhic::VirtualParticle& prt) const {
  // Check for genre first (em, hadronic, any)
  if (PGenre(prt) == 0 || (mGenre != 0 && PGenre(prt) != mGenre)) {
    return false;
  }  // if
  // Check if the particle charge matches the required value.
  if (mCharge != kAllCharges) { // Don't need to check if accepting all
    // Try to find the particle's charge via its TParticlePDG object.
    TParticlePDG* pdg = prt.Id().Info();
    if (pdg) {
      bool charged = fabs(pdg->Charge()) > 0.;
      // Check the charge against the requested value and return false
      // if it is incorrect.
      if ((kNeutral == mCharge && charged) ||
         (kCharged == mCharge && !charged)) {
        return false;
      }  // if
    } else {
      // The particle is unknown. We can't guarantee it's charge matches
      // the requested value, so return false.
      return false;
    }  // if
  }  // if
  // Check against exclusive particle list
  if (!mParticles.empty() && mParticles.count(prt.Id()) == 0) {
    return false;
  }  // if
  // If there are no Zones, accept everything that passed genre check
  if (mZones.empty()) {
    return true;
  }  // if
  for (unsigned i(0); i < mZones.size(); i++) {
    if (mZones.at(i).Contains(prt)) {
      return true;
    }  // if
  }  // for
  return false;
}

//
// class Acceptance::CustomCut
//

Acceptance::CustomCut::~CustomCut() {
}

Acceptance::CustomCut::CustomCut()
: mFormula("CustomCutFormula", "0")
, dim(0)
, Kin1(kP)
, Kin2(kTheta)
, Min(-TMath::Infinity())
, Max(TMath::Infinity()) {
}

Acceptance::CustomCut::CustomCut(const TString& formula,
                                 double min, double max)
: mFormula("CustomCutFormula", "0")
, dim(0)
, Kin1(kP)
, Kin2(kTheta)
, Min(min)
, Max(max) {
  TString s(formula);
  dim = ParseInputFunction(s, Kin1, Kin2);
  if (!IsCoreType(Kin1) || !IsCoreType(Kin2)) {
    std::cerr <<
    "ERROR! Custom acceptance is not a function of E, p, theta, phi"
    << std::endl;
  }  // if
  if (1 == dim || 2 == dim) {
    mFormula = TFormula("CustomCutFormula", s);
  } else {
    std::cerr <<
    "ERROR! Provided custom acceptance is not of dimension 1 or 2."
    << std::endl;
    return;
  }  // if
  std::cout << "Added custom cut " << formula << std::endl;
}

bool Acceptance::CustomCut::Contains(
                                     const erhic::VirtualParticle& prt) const {
  double x = GetVariable(prt, Kin1);
  double y(0.);
  if (2 == dim) {
    y = GetVariable(prt, Kin2);
  }  // if
  double z = mFormula.Eval(x, y);
  return z >= Min && z < Max;
}

//
// class Acceptance::Zone
//

Acceptance::Zone::~Zone() {
}

Acceptance::Zone::Zone(double thMin, double thMax,
                       double phMin, double phMax,
                       double eMin, double eMax,
                       double pMin, double pMax,
                       double ptmin, double ptmax,
                       double pzmin, double pzmax)
: thetaMin(thMin)
, thetaMax(thMax)
, phiMin(phMin)
, phiMax(phMax)
, EMin(eMin)
, EMax(eMax)
, PMin(pMin)
, PMax(pMax)
, pTMin(ptmin)
, pTMax(ptmax)
, pZMin(pzmin)
, pZMax(pzmax) {
}

void Acceptance::Zone::Add(const CustomCut& cut) {
  CustomCuts.push_back(cut);
}

Bool_t Acceptance::Zone::Contains(const erhic::VirtualParticle& prt) const {
  bool accept(true);
  const double theta = FixTheta(prt.GetTheta());
  const double phi = FixPhi(prt.GetPhi());
  if (theta < thetaMin || theta > thetaMax) {
    accept = false;
  } else if (phi < phiMin || phi > phiMax) {
    accept = false;
  } else if (prt.GetE() < EMin || prt.GetE() > EMax) {
    accept = false;
  } else if (prt.GetP() < PMin || prt.GetP() > PMax) {
    accept = false;
  } else if (prt.GetPz() < pZMin || prt.GetPz() > pZMax) {
    accept = false;
  } else if (prt.GetPt() < pTMin || prt.GetPt() > pTMax) {
    accept = false;
  }  // if
  // If it made it this far, test the custom cut(s)
  if (accept) {
    for (unsigned j(0); j < CustomCuts.size(); ++j) {
      if (!CustomCuts.at(j).Contains(prt)) {
        accept = false;
        break;
      }  // if
    }  // for
  }  // if
  return accept;
}

}  // namespace Smear
