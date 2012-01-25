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

//#ifdef USE_NAMESPACE_ERHIC   
namespace erhic {
//#endif
   
   int EventMC::smCount(0);
   
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
   ++smCount;
   }
   
   
   EventMC::~EventMC() {
      --smCount;
      for(unsigned i(0); i < particles.size(); ++i) {
         if(GetTrack(i)) {
            delete GetTrack(i);
         } // if
      } // for
      particles.clear();
   }
   
   Bool_t
   EventMC::Compute() {
      
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
      
      // Calculate event-dependent particle quantities.
      // Do this before calculating event kinematic variables via hadronic
      // methods as those rely on the particles.
      
      for(unsigned n(0); n < GetNTracks(); ++n ) {
         GetTrack(n)->ComputeEventDependentQuantities(*this);
      } // for(particles )
      
      ComputeJaquetBlondel(beams);
      ComputeDoubleAngle(beams);
      
      return kTRUE;
   }
   
   
   void
   EventMC::ComputeJaquetBlondel(const std::vector<const TrackType*>& beams ) {
      
//      std::cout << "EventMC::ComputeJaquetBlondel()" << std::endl;
      
      ::JacquetBlondel jacquetBlondel;
      jacquetBlondel.setBeamLepton(beams.at(0)->Get4Vector());
      jacquetBlondel.setBeamHadron(beams.at(1)->Get4Vector());
      
      std::vector<const TrackType*> final;
      HadronicFinalState(final);

      for(unsigned i(0); i < final.size(); ++i ) {
         jacquetBlondel.addParticle(final.at(i)->Get4Vector() );
      } // for
//      std::cout << "Added " << final.size() << " final state hadrons to the JB computer" << std::endl;
      
      yJB = jacquetBlondel.computeY();
      QSquaredJB = jacquetBlondel.computeQSquared();
      xJB = jacquetBlondel.computeX();
      /** \todo WSquaredJB, WSquaredDA calculations */
//      std::cout << QSquaredJB << " " << QSquared << std::endl;
//      assert((QSquaredJB - QSquared) / QSquared < 1.e-4);
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
      
//      std::cout << "EventMC::HadronicFinalState()" << std::endl;
      
      std::vector<const TrackType*> vp;
      ParticleIdentifier pi;
      if(not pi.IdentifyBeams(*this, vp)) {
         std::cerr <<
         "EventMC::HadronicFinalState(): failed to find all beams"
         << std::endl;
      } // if
/*
      std::cout << "Found beams:" << std::endl;
      std::cout << "Beam lepton:" << std::endl;
      vp.at(0)->Dump();
      std::cout << "Beam hadron:" << std::endl;
      vp.at(1)->Dump();
      std::cout << "Exchange boson:" << std::endl;
      vp.at(2)->Dump();
      std::cout << "Scattered lepton:" << std::endl;
      vp.at(3)->Dump();
      */
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
            //            std::cout << "Excluded scattered lepton for hadronic final state"
            //            << std::endl;
            final.erase(i);
            break;
         } // if
      } // for
      
      /*
      for(Iter i = final.begin(); i not_eq final.end(); ++i ) {
         // Try to skip the proton beam remnant
         if(final.at(1) and
            (*i)->Id() == final.at(1)->Id() and
            (*i)->GetP() > 0.7 * final.at(1)->GetP()) {
            final.erase(i);
            break;
         } // if
         
      } // for
       */
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
   
   
   struct Pair {
      int a;
      int b;
      Pair(int x, int y) : a(x), b(y) { }
      /**
       Less-than operator required for entry in sorted container.
       */
      bool operator<(const Pair& other) const {
         if(other.a == a) return b < other.b;
         return a < other.a;
      }
   };
   
#if 0
   void
   EventMC::Dot(const std::string& name) const {
      
      std::ofstream file(name.c_str());
      std::ostringstream oss;
      
      file << "digraph G {" << std::endl;
      oss << "   label = \"Event " << number << "\"';";
      file << oss.str() << std::endl;
      
      std::set<Pair> pairs;
      
      // Keep track of which particles we use in the diagram
      // so we can set node attributes at the end.
      // Use a set to avoid worrying about adding duplicates.
      std::set<const ParticleMC*> used;
      
      // Loop over all tracks in the event and accumulate indices giving
      // parent/child relationships.
      // We look from parent->child and child->parent because they
      // aren't always fully indexed both ways.
      for(unsigned i(0); i < GetNTracks(); ++i) {
         // Check parent of particle
         const ParticleMC* parent = GetTrack(i)->GetParent();
         if(parent) {
            pairs.insert(Pair(parent->GetIndex(), GetTrack(i)->GetIndex()));
            used.insert(parent);
            used.insert(GetTrack(i));
         } // if
         // Check children of particle
         for(unsigned j(0); j < GetTrack(i)->GetNChildren(); ++j) {
            pairs.insert(Pair(GetTrack(i)->GetIndex(),
                              GetTrack(i)->GetChild(j)->GetIndex()));
            used.insert(GetTrack(i));
            used.insert(GetTrack(i)->GetChild(j));
         } // for
      } // for
      
      // Insert relationships between particles.
      for(std::set<Pair>::const_iterator i = pairs.begin();
          i not_eq pairs.end();
          ++i) {
         const ParticleMC* a = GetTrack(i->a - 1);
         const ParticleMC* b = GetTrack(i->b - 1);
         oss.str("");
         oss << "   " << a->GetIndex() << " -> " <<
         b->GetIndex();
         file << oss.str() << std::endl;
      } // for
      
      file << "#   Node attributes" << std::endl;
      
      // Apply labels, formatting to used particles.
      for(std::set<const ParticleMC*>::const_iterator i = used.begin();
          i not_eq used.end();
          ++i) {
         
         std::string shape("ellipse");
         if((*i)->GetStatus() == 1) {
            shape = "box";
         } // if
         if((*i)->GetIndex() < 3) {
            shape = "diamond";
         } // if
         
         oss.str("");
         oss << "   " <<
         (*i)->GetIndex() << " [label=\""
         << (*i)->GetIndex() << " "
         << (*i)->Id().Info()->GetName()
         << "\", shape=" << shape << "];";
         file << oss.str() << std::endl;
      } // for
      
      file << "}";
   }
#endif
#if 0
   TASImage* EventMC::GenerateImage(const std::string& name) const {
      Dot("tmp.gv");
      gSystem->Exec("dot ");
   }
#endif
   
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
   
   
#if 0
   template<typename T>
   typename EventMCFinalState<T>::ParticlePtrList
   EventMCFinalState<T>::HadronicFinalState(const VirtualEvent<Type>& event) const {
      ParticlePtrList particles;
      for(unsigned i(0); i < event.GetNTracks(); ++i) {
         const Type* p = event.GetTrack(i);
         if(1 == p->GetStatus() and p not_eq event.ScatteredLepton()) {
            particles.push_back(p);
         } // if
      } // for
      return particles;
   }
   
   
   template<typename T>
   typename EventMCFinalState<T>::ParticlePtrList
   EventMCFinalState<T>::FinalState(const VirtualEvent<Type>& event) const {
      ParticlePtrList particles;
      for(unsigned i(0); i < event.GetNTracks(); ++i) {
         const Type* p = event.GetTrack(i);
         if(1 == p->GetStatus()) {
            particles.push_back(p);
         } // if
      } // for
      return particles;
   }
   
   
   template<typename T>
   TLorentzVector
   EventMCFinalState<T>::FinalStateMomentum(const VirtualEvent<Type>& event) const {
      ParticlePtrList particles = FinalState(event);
      TLorentzVector p;
      typedef typename ParticlePtrList::const_iterator Iter;
      for(Iter i = particles.begin(); i not_eq particles.end(); ++i) {
         p += (*i)->Get4Vector();
      } // for
      return p;
   }
   
   
   template<typename T>
   TLorentzVector
   EventMCFinalState<T>::HadronicFinalStateMomentum(const VirtualEvent<Type>& event) const {
      ParticlePtrList particles = HadronicFinalState(event);
      TLorentzVector p;
      typedef typename ParticlePtrList::const_iterator Iter;
      for(Iter i = particles.begin(); i not_eq particles.end(); ++i) {
         p += (*i)->Get4Vector();
      } // for
      return p;
   }
   
   
   template<typename T>
   Double_t
   EventMCFinalState<T>::FinalStateCharge(const VirtualEvent<Type>& event) const {
      ParticlePtrList particles = FinalState(event);
      double charge(0.);
      typedef typename ParticlePtrList::const_iterator Iter;
      for(Iter i = particles.begin(); i not_eq particles.end(); ++i) {
         TParticlePDG* pdg = (*i)->Id().Info();
         if(pdg) {
            charge += pdg->Charge() / 3.;
         } // if
      } // for
      return charge;
   }
#endif
   //#ifdef USE_NAMESPACE_ERHIC   
} // namespace erhic
#if 0
namespace {
   erhic::EventMCFinalState<erhic::ParticleMC> pm;
   erhic::EventMCFinalState<Smear::ParticleMCS> ps;
}
#endif
  //#endif
