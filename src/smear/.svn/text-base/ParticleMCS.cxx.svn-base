/**
 \file
 Implementation of class Smear::ParticleMCS.
 
 \author    Michael Savastio
 \date      2011-08-19
 \copyright 2011 Brookhaven National Lab
 */

#include "eicsmear/smear/ParticleMCS.h"

#include <iomanip>
#include <iostream>

#include <TMath.h>

namespace Smear {

ParticleMCS::ParticleMCS()
: status(0)
, id(0)
, px(0.)
, py(0.)
, pz(0.)
, E(0.)
, pt(0.)
, p(0.)
, theta(0.)
, phi(0.) {
}

ParticleMCS::ParticleMCS(const TLorentzVector& ep, int pdg, int stat)
: status(stat)
, id(pdg)
, px(ep.Px())
, py(ep.Py())
, pz(ep.Pz())
, E(ep.E())
, pt(ep.Pt())
, p(ep.P())
, theta(ep.Theta())
, phi(ep.Phi()) {
}

ParticleMCS::~ParticleMCS() {
}

TLorentzVector ParticleMCS::Get4Vector() const {
  return TLorentzVector(px, py, pz, E);
}

void ParticleMCS::Print(Option_t* /* unused */) const {
  // Store initial flags for cout so we can restore them later.
  std::ios_base::fmtflags ff = std::cout.flags();
  // Configure flags for output
  std::cout.setf(std::ios_base::fixed, std::ios_base::basefield);
  std::cout.precision(4);
  // Output values.
  std::cout <<
  std::setw(3) << status <<
  std::setw(12) << id <<
  std::setw(10) << px <<
  std::setw(10) << py <<
  std::setw(10) << pz <<
  std::setw(10) << E << std::endl;
  // Restore initial cout flags.
  std::cout.flags(ff);
}

Double_t ParticleMCS::GetEta() const {
  // The default value of -19 is used in the main eRHIC code,
  // so use that for consistency.
  double eta(-19.);
  const double theta = GetTheta();
  if (theta > 0. && theta < TMath::Pi() && !TMath::IsNaN(theta)) {
    eta = -log(tan(theta / 2.));
  }  // if
  return eta;
}

Double_t ParticleMCS::GetRapidity() const {
  // The default value of -19 is used in the main eRHIC code,
  // so use that for consistency.
  double y(-19.);
  // Be careful here, as the energy and momentum can become
  // smeared such that the stored E < pz, giving NaN when
  // taking a log of a negative.
  // In that case, return the default value.
  // E > pz already takes care of the case if E is NaN, so
  // no need to check that again.
  if (E > pz && !TMath::IsNaN(pz)) {
    y = 0.5 * log((E + pz) / (E - pz));
  }  // if
  return y;
}

}  // namespace Smear
