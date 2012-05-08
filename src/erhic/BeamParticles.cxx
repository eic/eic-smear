//
// BeamParticles.cxx
//
// Created by TB on 6/24/11.
// Copyright 2011 BNL. All rights reserved.
//

#include <cmath>

#include "eicsmear/erhic/BeamParticles.h"

//ClassImp(BeamParticles )

BeamParticles::BeamParticles() {
   Reset();
}

BeamParticles::BeamParticles(const TLorentzVector& hadronBeam,
                              const TLorentzVector& leptonBeam,
                              const TLorentzVector& scatteredHadron,
                              const TLorentzVector& scatteredLepton,
                              const TLorentzVector& exchangedBoson )
: mBeamHadron(hadronBeam )
, mBeamLepton(leptonBeam )
, mScatteredHadron(scatteredHadron )
, mScatteredLepton(scatteredLepton )
, mBoson(exchangedBoson ) {
}

BeamParticles::~BeamParticles() {
}

void
BeamParticles::Reset() {
   mBeamHadron.SetXYZT(NAN, NAN, NAN, NAN );
   mBeamLepton.SetXYZT(NAN, NAN, NAN, NAN );
   mScatteredHadron.SetXYZT(NAN, NAN, NAN, NAN );
   mScatteredLepton.SetXYZT(NAN, NAN, NAN, NAN );
   mBoson.SetXYZT(NAN, NAN, NAN, NAN );
}
