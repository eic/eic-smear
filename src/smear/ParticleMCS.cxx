//
// ParticleMCSmeared.cpp
//
// Created by TB on 10/10/11.
// Copyright 2011 BNL. All rights reserved.
//

#include <limits>

#include "eicsmear/smear/ParticleMCS.h"

namespace Smear {
   ParticleMCS::ParticleMCS()
   : status(0)
   , id(std::numeric_limits<Int_t>::max())
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
} // namespace Smear
