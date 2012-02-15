//
// EventMC.cxx
//
// Created by TB on 10/10/11.
// Copyright 2011 BNL. All rights reserved.
//

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>

#include <TASImage.h>
#include <TDatabasePDG.h>
#include <TDirectory.h>
#include <TParticlePDG.h>
#include <TSystem.h>
#include <TTree.h>

#include "EventMC.h"
#include "Kinematics.h"
#include "ParticleIdentifier.h"
#include "ParticleMC.h"

namespace erhic {
   
   EventMC::EventMC()
   : number(-1)
   , process(-1)
   , nTracks(-1)
   , ELeptonInNucl(NAN)
   , ELeptonOutNucl(NAN)
   , x(NAN)
   , QSquared(NAN)
   , y(NAN)
   , WSquared(NAN)
   , nu(NAN)
   , yJB(NAN)
   , QSquaredJB(NAN)
   , xJB(NAN)
   , WSquaredJB(NAN)
   , yDA(NAN)
   , QSquaredDA(NAN)
   , xDA(NAN)
   , WSquaredDA(NAN)
   {
   }
   
   
   EventMC::~EventMC() {
      for(unsigned i(0); i < particles.size(); ++i) {
         if(GetTrack(i)) {
            delete GetTrack(i);
         } // if
      } // for
      particles.clear();
   }
   
   
   Bool_t
   EventMC::Compute() {
      
      // Find the beams, exchange boson, scattered lepton
      std::vector<const TrackType*> beams;
      if(not ParticleIdentifier::IdentifyBeams(*this, beams)) {
         std::cerr << "EventMC::Compute(): failed to find beams" << std::endl;
         return kFALSE;
      } // if
      
      const TLorentzVector& lepton = beams.at(0)->Get4Vector();
      const TLorentzVector& hadron = beams.at(1)->Get4Vector();
      const TLorentzVector& scat = beams.at(3)->Get4Vector();
      
      QSquared = 2.0 * lepton.E() * scat.E() * (1.0 + scat.CosTheta() );
      
      double gamma = hadron.Gamma();
      double beta = hadron.Beta();
      
      ELeptonInNucl = gamma * (lepton.E() - beta * lepton.Pz());
      ELeptonOutNucl = gamma * (scat.E() - beta * scat.Pz());
      
      nu = ELeptonInNucl - ELeptonOutNucl;
      
      const double mass = hadron.M();
      
      x = QSquared / (2. * mass * nu );
      
      y = nu / ELeptonInNucl;
      
      WSquared = mass * mass + (1. - x ) * QSquared / x;
      
      // Calculate event-dependent particle quantities for each particle.
      // Do this before calculating event kinematic variables via hadronic
      // methods, as those rely on the particles.
      for(unsigned n(0); n < GetNTracks(); ++n ) {
         GetTrack(n)->ComputeEventDependentQuantities(*this);
      } // for(particles )
      
      ComputeJaquetBlondel(beams);
      ComputeDoubleAngle(beams);
      
      return kTRUE;
   }
   
   
   void
   EventMC::ComputeJaquetBlondel(const std::vector<const TrackType*>& beams ) {
      
      ::JacquetBlondel jacquetBlondel;
      jacquetBlondel.setBeamLepton(beams.at(0)->Get4Vector());
      jacquetBlondel.setBeamHadron(beams.at(1)->Get4Vector());
      
      std::vector<const TrackType*> final;
      HadronicFinalState(final);

      for(unsigned i(0); i < final.size(); ++i ) {
         jacquetBlondel.addParticle(final.at(i)->Get4Vector() );
      } // for
      
      yJB = jacquetBlondel.computeY();
      QSquaredJB = jacquetBlondel.computeQSquared();
      xJB = jacquetBlondel.computeX();
   }
   
   
   void
   EventMC::ComputeDoubleAngle(const std::vector<const TrackType*>& beams ) {
      
      ::DoubleAngle doubleAngle;
      
      doubleAngle.setBeamLepton(beams.at(0)->Get4Vector());
      doubleAngle.setBeamHadron(beams.at(1)->Get4Vector());
      doubleAngle.setLeptonAngle(beams.at(3)->GetTheta() );
      
      std::vector<const TrackType*> final;
      HadronicFinalState(final);
      for(unsigned i(0); i < final.size(); ++i ) {
         doubleAngle.addParticle(final.at(i)->Get4Vector() );
      } // for
      
      yDA = doubleAngle.computeY();
      QSquaredDA = doubleAngle.computeQSquared();
      xDA = doubleAngle.computeX();
   }
   
   
   // Get the particles that belong to the hadronic final state.
   // The stored Particle* are pointers to the original particles in the event
   // so don't delete them!
   void
   EventMC::HadronicFinalState(std::vector<const TrackType*>& final) const {
      
      std::vector<const TrackType*> vp;
      ParticleIdentifier pi;
      if(not pi.IdentifyBeams(*this, vp)) {
         std::cerr <<
         "EventMC::HadronicFinalState(): failed to find all beams"
         << std::endl;
      } // if
      
      FinalState(final);
      
      typedef std::vector<const TrackType*>::iterator Iter;
      
      // Find the scattered lepton, which we define as the first particle with
      // the same species as the lepton beam
      // TODO figure out a better way to specify the lepton pdg code
      for(Iter i = final.begin(); i not_eq final.end(); ++i ) {
         if(not *i) {
            std::cerr << "NULL particle in array?!" << std::endl;
            continue;
         } // if
         if(*i == vp.at(3)) {
            final.erase(i);
            break;
         } // if
      } // for
   }
   
   
   // Get the particles that belong to the hadronic final state.
   // The stored Particle* are pointers to the original particles in the event
   // so don't delete them!
   void
   EventMC::FinalState(std::vector<const TrackType*>& final_ ) const {
      
      final_.clear();
      
      typedef std::vector<TrackType*>::const_iterator Iter;
      
      for(Iter i_ = particles.begin(); i_ not_eq particles.end(); ++i_ ) {
         if(1 == (*i_)->GetStatus() ) {
            final_.push_back(*i_);
         } // if
      } // for
   }
   
   
   TLorentzVector
   EventMC::FinalStateMomentum() const {
      
      std::vector<const TrackType*> final_;
      FinalState(final_ );
      
      typedef std::vector<const TrackType*>::const_iterator Iter;
      
      TLorentzVector mom;
      
      for(Iter i = final_.begin(); i not_eq final_.end(); ++i ) {
         mom += (*i)->Get4Vector();
      } // for
      
      return mom;
   }
   
   
   TLorentzVector
   EventMC::HadronicFinalStateMomentum() const {
      
      std::vector<const TrackType*> final_;
      HadronicFinalState(final_ );
      
      typedef std::vector<const TrackType*>::const_iterator Iter;
      
      TLorentzVector mom;
      
      for(Iter i = final_.begin(); i not_eq final_.end(); ++i ) {
         mom += (*i)->Get4Vector();
      } // for
      
      return mom;
   }
   
   
   Double_t EventMC::FinalStateCharge() const {
      std::vector<const TrackType*> final_;
      FinalState(final_);
      typedef std::vector<const TrackType*>::const_iterator Iter;
      Double_t charge(0);
      TDatabasePDG* pdg = TDatabasePDG::Instance();
      for(Iter i = final_.begin(); i not_eq final_.end(); ++i) {
         TParticlePDG* part = pdg->GetParticle((*i)->Id());
         if(part) {
            charge += part->Charge() / 3.;
         } // if
         else {
            std::cout << "Unknown particle: " << (*i)->Id() << std::endl;
         } // else
      } // for
      return charge;
   }
   
   
   const EventMC::TrackType* EventMC::BeamLepton() const {
      return (particles.empty() ? NULL : GetTrack(0));
   }
   
   
   const EventMC::TrackType* EventMC::BeamHadron() const {
      return (particles.size() > 1 ? GetTrack(1) : NULL);
   }
   
   
   const EventMC::TrackType* EventMC::ExchangeBoson() const {
      return (particles.size() > 3 ? GetTrack(3) : NULL);
   }
   
   
   const EventMC::TrackType* EventMC::ScatteredLepton() const {
      return (particles.size() > 2 ? GetTrack(2) : NULL);
   }
   
   
   void EventMC::Reset() {
      number = -1;
      process = -1;
      nTracks = -1;
      x = QSquared = y = WSquared = nu = ELeptonInNucl = ELeptonOutNucl = NAN;
      particles.clear();
   }
   
   
   //
   // class Reader
   //
   
   
   Reader::Reader(const std::string& treeName)
   : mEvent(NULL)
   , mTree(NULL) {
      if(gDirectory) gDirectory->GetObject(treeName.c_str(), mTree);
      if(mTree) mTree->SetBranchAddress("event", &mEvent);
   }
   
   
   EventMC* Reader::Read(Long64_t i) {
      EventMC* event(NULL);
      if(mTree) {
         mTree->GetEntry(i);
         event = mEvent;
      } // if
      return event;
   }
   
   
   void EventMC::AddLast(TrackType* track) {
      particles.push_back(track);
      nTracks = particles.size();
   }
   
} // namespace erhic
