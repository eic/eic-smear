/**
 \file
 Implementation of class erhic::Pythia6ParticleBuilder.
 
 \author    Thomas Burton
 \date      2012-01-17
 \copyright 2012 Brookhaven National Lab
 */

#include <iostream>

#include <TLorentzVector.h>
#include <TMCParticle.h>
#include <TVector3.h>

#include "eicsmear/erhic/ParticleMC.h"
#include "eicsmear/erhic/Pythia6ParticleBuilder.h"

namespace erhic {

std::auto_ptr<ParticleMC>
Pythia6ParticleBuilder::Create(const TMCParticle& mc) const {
  std::auto_ptr<ParticleMC> particle(new ParticleMC);
  TLorentzVector p(mc.GetPx(),
                   mc.GetPy(),
                   mc.GetPz(),
                   mc.GetEnergy());
  TVector3 v(mc.GetVx(),
             mc.GetVy(),
             mc.GetVz());
  particle->SetStatus(mc.GetKS());
  particle->SetId(mc.GetKF());
  particle->SetParentIndex(mc.GetParent());
  particle->SetChild1Index(mc.GetFirstChild());
  particle->SetChildNIndex(mc.GetLastChild());
  std::cout << "ChildN index " << mc.GetLastChild() << " " <<
    particle->GetChildNIndex() << std::endl;
  particle->Set4Vector(p);
  particle->SetVertex(v);
  return particle;
}

}  // namespace erhic
