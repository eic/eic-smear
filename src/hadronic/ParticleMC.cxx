/**
 \file
 Implementation of class erhic::hadronic::ParticleMC.
 
 \author    Thomas Burton 
 \date      2012-05-02
 \copyright 2012 Brookhaven National Lab
 */

#include <eicsmear/hadronic/ParticleMC.h>

#include <cmath>
#include <limits>

#include <TMCParticle.h>
#include <TVector2.h>

namespace erhic {
namespace hadronic {

ParticleMC::ParticleMC()
: KS(-1)
, orig(-1)
, id(std::numeric_limits<Int_t>::min())
, px(NAN)
, py(NAN)
, pz(NAN)
, E(NAN)
, p(NAN)
, m(NAN)
, pt(NAN)
, theta(NAN)
, phi(NAN)
, rapidity(NAN)
, eta(NAN)
, xFeynman(NAN)
, xv(NAN)
, yv(NAN)
, zv(NAN) {
}

ParticleMC::ParticleMC(const TMCParticle& mc)
: KS(mc.GetKS())
, orig(-1)
, id(mc.GetKF())
, px(mc.GetPx())
, py(mc.GetPy())
, pz(mc.GetPz())
, E(mc.GetEnergy())
, m(mc.GetMass())
, pt(NAN)
, theta(NAN)
, rapidity(NAN)
, eta(NAN)
, xFeynman(NAN)
, xv(mc.GetVx())
, yv(mc.GetVy())
, zv(mc.GetVz()) {
  TLorentzVector v(mc.GetPx(), mc.GetPy(), mc.GetPz(), mc.GetEnergy());
  p = v.P();
  pt = v.Pt();
  theta = v.Theta();
  phi = TVector2::Phi_0_2pi(v.Phi());
  rapidity = v.Rapidity();
  // Protect against attempting pseudorapidity calculation for
  // particles with angle close to the beam.
  if (v.Pt() > 1.e-5) {
    eta = v.PseudoRapidity();
  } else {
    eta = std::numeric_limits<Double32_t>::infinity();
  }  // if
}

ParticleMC::ParticleMC(const TLorentzVector& ep, const TVector3& v, int pdg,
                       int status, int parent)
: KS(status)
, orig(parent)
, id(pdg)
, px(ep.Px())
, py(ep.Py())
, pz(ep.Pz())
, E(ep.E())
, p(ep.P())
, m(ep.M())
, pt(ep.Pt())
, theta(ep.Theta())
, phi(TVector2::Phi_0_2pi(ep.Phi()))
, rapidity(ep.Rapidity())
, eta(ep.PseudoRapidity())
, xFeynman(NAN)
, xv(v.X())
, yv(v.Y())
, zv(v.Z()) {
}


void ParticleMC::SetParentIndex(UShort_t i) {
  orig = i;
}

void ParticleMC::Set4Vector(const TLorentzVector& v) {
  E = v.Energy();
  px = v.Px();
  py = v.Py();
  pz = v.Pz();
  p = v.P();
  m = v.M();
  pt = v.Pt();
  theta = v.Theta();
  phi = TVector2::Phi_0_2pi(v.Phi());
  rapidity = v.Rapidity();
  eta = v.PseudoRapidity();
}

}  // namespace hadronic
}  // namespace erhic
