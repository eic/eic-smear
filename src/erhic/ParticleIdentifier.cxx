/**
 \file
 Implementation of class ParticleIdentifier.
 
 \author    Thomas Burton
 \date      2011-10-10
 \copyright 2011 Brookhaven National Lab
 */

#include "eicsmear/erhic/ParticleIdentifier.h"

#include <algorithm>
#include <iostream>
#include <vector>
#include <string>

#include <TParticlePDG.h>
#include <TDatabasePDG.h>

#include "eicsmear/erhic/EventMC.h"

using std::string;
using std::cout;
using std::cerr;
using std::endl;


// =============================================================================
// Constructor
// =============================================================================
ParticleIdentifier::ParticleIdentifier(const int leptonBeamPdgCode)
: mChargedCurrent(false)
, mLeptonBeamPdgCode(leptonBeamPdgCode)
, mScatteredPdgCode(leptonBeamPdgCode)
{ }

// =============================================================================
// Identify the lepton beam particle
// =============================================================================
bool ParticleIdentifier::isBeamLepton(
         const erhic::VirtualParticle& particle) const {
  // Beam?
  if ( particle.GetParentIndex() != 0 ) return false;
  // Lepton?
  if ( particle.Id() != GetLeptonBeamPdgCode() ) return false;
  // Beam status?
  switch (particle.GetStatus()){
  case 21 : return true;  // pythia6 et al
  case 201 : return true; // SOPHIA
  case 4: return true;    // HepMC
  }
  return false;
}

// =============================================================================
// Identify the scattered lepton
// =============================================================================
bool ParticleIdentifier::isScatteredLepton(const erhic::VirtualParticle& particle) const {
  auto pdgl = TDatabasePDG::Instance()->GetParticle( particle.Id() );
  
  return ( particle.GetStatus() == 1  && string(pdgl->ParticleClass()) == "Lepton" );
  // the old version here ignores flavor change, such as charged current dis
  // return ( particle.GetStatus() == 1  && particle.Id()==mScatteredPdgCode);
}

// =============================================================================
// Return true if the particle matches any of the following:
//  - Outgoing electron with KS 21 - we instead pick up the
//                                   repeat entry, when KS < 21.
//  - Repeat of incoming nucleon.
//  - Intermediate gluon.
//  - Intermediate quark.
// =============================================================================
bool ParticleIdentifier::SkipParticle(
         const erhic::VirtualParticle& particle) const {
  int kI1 = particle.GetStatus();
  int pdgCode = particle.Id();
  int parent = particle.GetParentIndex();
  // Remove duplicate outgoing electron, info is picked up later when KS < 21
  if (21 == kI1 && pdgCode == GetLeptonBeamPdgCode() && parent == 1) {
    return true;
  }  // if
    // Remove repeat of incoming nucelon
  if (21 == kI1 && (pdgCode == 2112 || pdgCode == 2212) && parent == 2) {
    return true;
  }  // if
    // Remove intermediate gluons
  if (21 == kI1 && pdgCode == 21) {
    return true;
  }  // if
    // Remove intermediate (anti)quarks
  if (21 == kI1 && ::abs(pdgCode) < 10) {
    return true;
  }  // if
  return false;
}

// =============================================================================
// Identify the virtual photon based on its PDG code and PYTHIA K(I,1)
// =============================================================================
bool ParticleIdentifier::IsVirtualPhoton(
         const erhic::VirtualParticle& particle) const {
  const int pdg = abs(particle.Id());
  auto status = particle.GetStatus();
  // special case for placeholder particle in ff2ff events
  if ( particle.GetStatus() == 0 && status ==0 ) return true;
  
  if (pdg<22) return false;
  if (pdg>24) return false;
  if ( status == 21 ) return true; // pythia6 et al
  if ( status == 13 ) return true; // pythia8
  return false;
}

// =============================================================================
// Identify the nucleon beam particle
// =============================================================================
bool ParticleIdentifier::isBeamNucleon(
         const erhic::VirtualParticle& particle) const {
  // Beam?
  if ( particle.GetParentIndex() != 0 ) return false;
  // Hadron?
  if ( particle.Id() != 2112 && particle.Id() != 2212 // e+P
       && abs(particle.Id()) < 1000000000 ) return false; // allow ions in eA
  // Beam status?
  switch (particle.GetStatus()){
  case 21 : return true;  // pythia6 et al
  case 201 : return true; // SOPHIA
  case 4: return true;    // HepMC
  }

  return false;
}

// =============================================================================
// beam electron (11) --> scattered nu_e (12)
// beam positron (-11) --> scattered nu_e_bar (-12)
// =============================================================================
Int_t ParticleIdentifier::DetermineScatteredType(Int_t beamType) {
  if (!mChargedCurrent) {
    return beamType;
  }  // if
  int sign = beamType / abs(beamType);
  return sign * (abs(beamType) + 1);
}

// =============================================================================
// =============================================================================
void ParticleIdentifier::SetLeptonBeamPdgCode(const int pdgCode) {
  mLeptonBeamPdgCode = pdgCode;
  if (mChargedCurrent) {
    mScatteredPdgCode = DetermineScatteredType(pdgCode);
  } else {
    mScatteredPdgCode = pdgCode;
  }  // if
}

bool ParticleIdentifier::SetChargedCurrent(bool cc) {
  // If charged current status is changed we need to update
  // the expected scattered lepton type. But we need to call
  // DetermineScatteredType() *after* updating mChargedCurrent.
  bool changed = !mChargedCurrent && cc;
  mChargedCurrent = cc;
  if (changed) {
    mScatteredPdgCode = DetermineScatteredType(mLeptonBeamPdgCode);
  }  // if
  return cc;
}

// =============================================================================
// Identify the beams from an event and store their properties in a
// BeamParticles object.
// See BeamParticles.h for the quantities stored.
// Returns true if all beams are found, false if not.
// =============================================================================
bool ParticleIdentifier::IdentifyBeams(const erhic::VirtualEvent& event,
                                       BeamParticles& beams) {
  beams.Reset();
  std::vector<const erhic::VirtualParticle*> particles;
  bool foundAll = IdentifyBeams(event, particles);
  if (particles.at(0)) {
    beams.SetBeamLepton(particles.at(0)->Get4Vector());
  }  // if
  if (particles.at(1)) {
    beams.SetBeamHadron(particles.at(1)->Get4Vector());
  }  // if
  if (particles.at(2)) {
    beams.SetBoson(particles.at(2)->Get4Vector());
  }  // if
  if (particles.at(3)) {
    beams.SetScatteredLepton(particles.at(3)->Get4Vector());
  }  // if
  return foundAll;
}

// =============================================================================
// =============================================================================
bool ParticleIdentifier::IdentifyBeams(const erhic::VirtualEvent& event,
         std::vector<const erhic::VirtualParticle*>& beams) {
  // Initialise a vector with four NULL pointers,
  // one beam particle of interest.
  const erhic::VirtualParticle* const null(NULL);
  beams.assign(4, null);
  // Configure the particle finder.
  ParticleIdentifier finder;
  if (event.GetTrack(0)) {
    finder.SetLeptonBeamPdgCode(event.GetTrack(0)->Id());
  }  // if
  // Count leptons so we don't overwrite the beam and scattered lepton with
  // subsequent leptons of the same type.
  int leptonCount(0);
  // Set to true once we find the first virtual photon, so we can skip
  // subsequent virtual photons.
  bool foundExchangeBoson(false);
  bool foundAll(false);
  for (unsigned n(0); n < event.GetNTracks(); ++n) {
    const erhic::VirtualParticle* particle = event.GetTrack(n);
    if (!particle) {
      continue;
    }  // if
    // Test for beam lepton/hadron, exchange boson and scattered lepton.
    if (finder.isBeamNucleon(*particle)) {
      // cout << "Found the beam nucleon" << endl;
      beams.at(1) = particle;
    } else if (finder.isBeamLepton(*particle) && 0 == leptonCount) {
      // cout << "Found the beam lepton" << endl;
      beams.at(0) = particle;
      ++leptonCount;
    } else if (finder.isScatteredLepton(*particle) && 1 == leptonCount) {
      // cout << "Found the scattered lepton" << endl;
      beams.at(3) = particle;
      // Protect against additional KS == 1 leptons following this
      ++leptonCount;
    } else if (finder.IsVirtualPhoton(*particle) && !foundExchangeBoson) {
      // cout << "Found the boson" << endl;
      beams.at(2) = particle;
      foundExchangeBoson = true;
      // Check for charged current events, in which the scattered lepton
      // ID will not be the same as the incident lepton ID.
      finder.SetChargedCurrent(abs(particle->Id()) == 24);
    }  // if
    // Break out if we've found all four beams (should happen at/near the
    // start of the event record, so checking the other particles is a
    // waste of time).
    foundAll = std::find(beams.begin(), beams.end(), null) == beams.end();
    if (foundAll) {
      break;
    }  // if
  }  // for
  return foundAll;
}
