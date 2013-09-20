/**
 \file
 Implementation of class erhic::EventMC.
 
 \author    Thomas Burton
 \date      2011-10-31
 \copyright 2011 Brookhaven National Lab
 */

#include "eicsmear/erhic/EventMC.h"

#include <iostream>
#include <list>
#include <string>
#include <vector>

#include <TCollection.h>
#include <TDatabasePDG.h>
#include <TDirectory.h>
#include <TParticlePDG.h>
#include <TTree.h>

#include "eicsmear/erhic/ParticleIdentifier.h"
#include "eicsmear/erhic/ParticleMC.h"

namespace erhic {

// We use these vectors/iterators a lot, so define
// some typedefs for convenience.
typedef std::vector<const VirtualParticle*> TrackVector;
typedef TrackVector::iterator TrackVectorIter;
typedef TrackVector::const_iterator TrackVectorCIter;

EventMC::EventMC()
: number(-1)
, process(-1)
, nTracks(-1)
, ELeptonInNucl(NAN)
, ELeptonOutNucl(NAN)
, particles("erhic::ParticleMC", 100) {
}

EventMC::~EventMC() {
  // No memory to clear. The TClonesArray takes care of the
  // particles when it is destroyed.
}

TrackVector EventMC::GetTracks() const {
  TrackVector tracks;
  TObject* object(NULL);
  TIter next(&particles);
  while ((object = next())) {
    tracks.push_back(static_cast<ParticleMC*>(object));
  }  // while
  return tracks;
}

void EventMC::HadronicFinalState(TrackVector& final_) const {
  // This is simple but not very efficient.
  // Copy vector to a list and use list::remove to get rid
  // of the scattered lepton. Then copy this list back to
  // the vector.
  FinalState(final_);
  std::list<const VirtualParticle*> plist(final_.begin(),
                                          final_.end());
  plist.remove(ScatteredLepton());
  final_ = TrackVector(plist.begin(), plist.end());
}

// Get the particles that belong to the hadronic final state.
// The stored Particle* are pointers to the original particles in the event
// so don't delete them!
void EventMC::FinalState(TrackVector& final_) const {
  final_.clear();
  TIter next(&particles);
  ParticleMC* p(NULL);
  while ((p = static_cast<ParticleMC*>(next()))) {
    if (1 == p->GetStatus()) {
      final_.push_back(p);
    }  // if
  }  // while
}

TLorentzVector EventMC::FinalStateMomentum() const {
  TrackVector final_;
  FinalState(final_);
  TLorentzVector mom;
  for (TrackVectorCIter i = final_.begin(); i != final_.end(); ++i) {
    mom += (*i)->Get4Vector();
  }  // for

  return mom;
}

TLorentzVector EventMC::HadronicFinalStateMomentum() const {
  TrackVector final_;
  HadronicFinalState(final_);
  TLorentzVector mom;
  for (TrackVectorCIter i = final_.begin(); i != final_.end(); ++i) {
    mom += (*i)->Get4Vector();
  }  // for

  return mom;
}

Double_t EventMC::FinalStateCharge() const {
  TrackVector final_;
  FinalState(final_);
  Double_t charge(0);
  TDatabasePDG* pdg = TDatabasePDG::Instance();
  for (TrackVectorCIter i = final_.begin(); i != final_.end(); ++i) {
    TParticlePDG* part = pdg->GetParticle((*i)->Id());
    if (part) {
      charge += part->Charge() / 3.;
    } else {
      std::cout << "Unknown particle: " << (*i)->Id() << std::endl;
    }  // if
  }  // for
  return charge;
}

const ParticleMC* EventMC::BeamLepton() const {
  return GetTrack(0);
}

const ParticleMC* EventMC::BeamHadron() const {
  return GetTrack(1);
}

const ParticleMC* EventMC::ExchangeBoson() const {
  return GetTrack(3);
}

const ParticleMC* EventMC::ScatteredLepton() const {
  return GetTrack(2);
}

void EventMC::Reset() {
  number = -1;
  process = -1;
  nTracks = -1;
  x = QSquared = y = WSquared = nu = ELeptonInNucl = ELeptonOutNucl = NAN;
}

void EventMC::Clear(Option_t* /*option*/) {
  Reset();
  particles.Clear();
}

void EventMC::AddLast(ParticleMC* track) {
  new(particles[particles.GetEntries()]) ParticleMC(*track);
  nTracks = particles.GetEntries();
}

//
// class Reader
//

Reader::Reader(const std::string& treeName)
: mEvent(NULL)
, mTree(NULL) {
  if (gDirectory) gDirectory->GetObject(treeName.c_str(), mTree);
  if (mTree) mTree->SetBranchAddress("event", &mEvent);
}

EventMC* Reader::Read(Long64_t i) {
  EventMC* event(NULL);
  if (mTree) {
    mTree->GetEntry(i);
    event = mEvent;
  }  // if
  return event;
}

}  // namespace erhic
