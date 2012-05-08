//
// EventMC.cxx
//
// Created by TB on 10/10/11.
// Copyright 2011 BNL. All rights reserved.
//

#include "eicsmear/erhic/EventMC.h"

#include <iostream>

#include <TDatabasePDG.h>
#include <TDirectory.h>
#include <TParticlePDG.h>
#include <TTree.h>

#include "eicsmear/erhic/ParticleIdentifier.h"
#include "eicsmear/erhic/ParticleMC.h"

namespace erhic {

   EventMC::EventMC()
   : number(-1)
   , process(-1)
   , nTracks(-1)
   , ELeptonInNucl(NAN)
   , ELeptonOutNucl(NAN) {
   }

   EventMC::~EventMC() {
      for(unsigned i(0); i < particles.size(); ++i) {
         if(GetTrack(i)) {
            delete GetTrack(i);
         } // if
      } // for
      particles.clear();
   }

   // Get the particles that belong to the hadronic final state.
   // The stored Particle* are pointers to the original particles in the event
   // so don't delete them!
   void
   EventMC::HadronicFinalState(std::vector<const VirtualParticle*>& final) const {
      
      std::vector<const VirtualParticle*> vp;
      ParticleIdentifier pi;
      if(not pi.IdentifyBeams(*this, vp)) {
         std::cerr <<
         "EventMC::HadronicFinalState(): failed to find all beams"
         << std::endl;
      } // if
      
      FinalState(final);
//      std::cout << "\t" << final.size() << std::endl;
      
      typedef std::vector<const VirtualParticle*>::iterator Iter;
      
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
//            std::cout << "erased" << std::endl;
            break;
         } // if
      } // for
//      std::cout << "\t" << final.size() << std::endl;
   }
   
   
   // Get the particles that belong to the hadronic final state.
   // The stored Particle* are pointers to the original particles in the event
   // so don't delete them!
   void
   EventMC::FinalState(std::vector<const VirtualParticle*>& final_ ) const {
      
      final_.clear();
      
      typedef std::vector<ParticleMC*>::const_iterator Iter;
      
      for(Iter i_ = particles.begin(); i_ not_eq particles.end(); ++i_ ) {
         if(1 == (*i_)->GetStatus() ) {
            final_.push_back(*i_);
         } // if
      } // for
   }
   
   
   TLorentzVector
   EventMC::FinalStateMomentum() const {
      
      std::vector<const VirtualParticle*> final_;
      FinalState(final_ );
      
      typedef std::vector<const VirtualParticle*>::const_iterator Iter;
      
      TLorentzVector mom;
      
      for(Iter i = final_.begin(); i not_eq final_.end(); ++i ) {
         mom += (*i)->Get4Vector();
      } // for
      
      return mom;
   }
   
   
   TLorentzVector
   EventMC::HadronicFinalStateMomentum() const {
      
      std::vector<const VirtualParticle*> final_;
      HadronicFinalState(final_ );
      
      typedef std::vector<const VirtualParticle*>::const_iterator Iter;
      
      TLorentzVector mom;
      
      for(Iter i = final_.begin(); i not_eq final_.end(); ++i ) {
         mom += (*i)->Get4Vector();
      } // for
      
      return mom;
   }
   
   
   Double_t EventMC::FinalStateCharge() const {
      std::vector<const VirtualParticle*> final_;
      FinalState(final_);
      typedef std::vector<const VirtualParticle*>::const_iterator Iter;
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
   
   
   const ParticleMC* EventMC::BeamLepton() const {
      return (particles.empty() ? NULL : GetTrack(0));
   }
   
   
   const ParticleMC* EventMC::BeamHadron() const {
      return (particles.size() > 1 ? GetTrack(1) : NULL);
   }
   
   
   const ParticleMC* EventMC::ExchangeBoson() const {
      return (particles.size() > 3 ? GetTrack(3) : NULL);
   }
   
   
   const ParticleMC* EventMC::ScatteredLepton() const {
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
   
   
   void EventMC::AddLast(ParticleMC* track) {
      particles.push_back(track);
      nTracks = particles.size();
   }
   
} // namespace erhic
