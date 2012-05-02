//
// ParticleMCSmeared.cpp
//
// Created by TB on 10/10/11.
// Copyright 2011 BNL. All rights reserved.
//

#include <limits>

#include "ParticleMCS.h"

namespace Smear {
   ParticleMCS::ParticleMCS()
   : status(0)
   , id(std::numeric_limits<Int_t>::max())
   , px(NAN)
   , py(NAN)
   , pz(NAN)
   , E(NAN)
   , pt(NAN)
   , p(NAN)
   , theta(NAN)
   , phi(NAN) {
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
