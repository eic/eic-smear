/**
 \file
 Implementation of class BeamParticles.
 
 \author    Thomas Burton
 \date      2011-06-24
 \copyright 2011 Brookhaven National Lab
 */

#include "eicsmear/erhic/BeamParticles.h"

#include <cmath>

BeamParticles::BeamParticles() {
  Reset();
}

BeamParticles::BeamParticles(const TLorentzVector& hadronBeam,
                             const TLorentzVector& leptonBeam,
                             const TLorentzVector& scatteredHadron,
                             const TLorentzVector& scatteredLepton,
                             const TLorentzVector& exchangedBoson)
: mBeamHadron(hadronBeam)
, mBeamLepton(leptonBeam)
, mScatteredHadron(scatteredHadron)
, mScatteredLepton(scatteredLepton)
, mBoson(exchangedBoson) {
}

BeamParticles::~BeamParticles() {
}

void BeamParticles::Reset() {
  mBeamHadron.SetXYZT(NAN, NAN, NAN, NAN);
  mBeamLepton.SetXYZT(NAN, NAN, NAN, NAN);
  mScatteredHadron.SetXYZT(NAN, NAN, NAN, NAN);
  mScatteredLepton.SetXYZT(NAN, NAN, NAN, NAN);
  mBoson.SetXYZT(NAN, NAN, NAN, NAN);
}
