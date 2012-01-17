//
// ParticleMCSmeared.cpp
//
// Created by TB on 10/10/11.
// Copyright 2011 BNL. All rights reserved.
//

#include <limits>

#include "ParticleMCSmeared.h"

namespace Smear {
   
   ParticleMCS::ParticleMCS()
   : id(std::numeric_limits<Int_t>::max())
   , pz(NAN)
   , E(NAN)
   , pt(NAN)
   , p(NAN)
   , theta(NAN)
   , phi(NAN)
   {
   }
   
   ParticleMCS::~ParticleMCS() {
   }
   
   TLorentzVector ParticleMCS::Get4Vector() const {
      return TLorentzVector(px, py, pz, E);
   }
   
} // namespace Smear
