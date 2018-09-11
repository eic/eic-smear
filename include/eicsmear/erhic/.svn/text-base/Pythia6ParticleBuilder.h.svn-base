/**
 \file
 Declaration of class erhic::Pythia6ParticleBuilder.
 
 \author    Thomas Burton
 \date      2012-01-17
 \copyright 2012 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_ERHIC_PYTHIA6PARTICLEBUILDER_H_
#define INCLUDE_EICSMEAR_ERHIC_PYTHIA6PARTICLEBUILDER_H_

#include <memory>

class TMCParticle;

namespace erhic {

class ParticleMC;

/**
 \brief Factory class for Monte Carlo particles.
 */
class Pythia6ParticleBuilder {
 public:
  /**
   Default constructor.
   */
  Pythia6ParticleBuilder() { }

  /**
   Generate a ParticleMC from a ROOT TMCParticle.
   */
  std::auto_ptr<ParticleMC> Create(const TMCParticle&) const;
};

}  // namespace erhic

#endif  // INCLUDE_EICSMEAR_ERHIC_PYTHIA6PARTICLEBUILDER_H_
