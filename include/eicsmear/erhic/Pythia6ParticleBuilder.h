/**
 Pythia6ParticleBuilder.h

 \file
 Declaration of class Pythia6ParticleBuilder.

 Created by TB on 1/17/12.
 Copyright 2012 BNL. All rights reserved.
*/

#ifndef _Pythia6ParticleBuilder_H_
#define _Pythia6ParticleBuilder_H_

#include <memory>

class TMCParticle;

namespace erhic {
   
   class ParticleMC;
   
   class Pythia6ParticleBuilder {
      
   public:
      
      Pythia6ParticleBuilder() { }
      
      std::auto_ptr<ParticleMC> Create(const TMCParticle&) const;
      
   };
   
} // namespace erhic

#endif
