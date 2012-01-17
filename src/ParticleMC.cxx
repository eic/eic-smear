//
// ParticleMC.cxx
//
// Created by TB on 10/10/11.
// Copyright 2011 BNL. All rights reserved.
//


#include <iostream>
#include <sstream>
//#include <string>

#include <TLorentzRotation.h>
#include <TRotation.h>

//#include <TRefArray.h>

#include "EventBase.h"
//#include "VirtualParticle.h"
#include "functions.h"

namespace erhic {
   
   
   Long64_t ParticleMC::smNInstances(0);

   
   ParticleMC::ParticleMC(const std::string& line)
   : I(-1)
   , KS(-1)
   , id(std::numeric_limits<Int_t>::min())
   , orig(-1)
   , daughter(-1)
   , ldaughter(-1)
   , px(NAN)
   , py(NAN)
   , pz(NAN)
   , E(NAN)
   , m(NAN)
   , pt(NAN)
   , xv(NAN)
   , yv(NAN)
   , zv(NAN)
   , parentId(std::numeric_limits<Int_t>::min())
   , p(NAN)
   , theta(NAN)
   , phi(NAN)
   , rapidity(NAN)
   , eta(NAN)
   , z(NAN)
   , xFeynman(NAN)
   , thetaGamma(NAN)
   , ptVsGamma(NAN)
   , phiPrf(NAN) {
      // Initialise to nonsense values to make input errors easy to spot
      
      static std::stringstream ss;
      ss.str("");
      ss.clear();
      
      ss << line;
      
      ss >>
      I >> KS >> id >> orig >> daughter >> ldaughter >>
      px >> py >> pz >> E >> m >> xv >> yv >> zv;
      
      if(not line.empty()) {
         ComputeDerivedQuantities();
      } // if
      
      ++smNInstances;
   }
   
   
   ParticleMC::~ParticleMC() {
      --smNInstances;
   }
   
   
   void ParticleMC::Print(Option_t*) const {
      std::cout << I << '\t' << KS << '\t' << id << '\t' << orig << '\t' <<
      daughter << '\t' << ldaughter << '\t' << px << '\t' << py << '\t' << pz
      << '\t' << E << '\t' << m << '\t' << xv << '\t' << yv << '\t' << zv <<
      std::endl;
   }
   
   
   void
   ParticleMC::ComputeDerivedQuantities() {
      
      // We calculate quantities that depend only on the properties already read.
      pt = sqrt(pow(px, 2.) + pow(py, 2.));
      p = sqrt(pow(pt, 2.) + pow(pz, 2.));
      
      // Rapidity and pseudorapidity
      Double_t Epluspz = E + pz;
      Double_t Eminuspz = E - pz;
      Double_t Ppluspz = p + pz;
      Double_t Pminuspz = p - pz;
      if(Eminuspz <= 0.0 or Pminuspz == 0.0 or
         Ppluspz == 0.0 or Epluspz <= 0.0) {
         // Dummy values to avoid zero or infinite arguments in calculations
         rapidity = -19.;
         eta = -19.;
      }	//	if...
      else {
         rapidity = 0.5 * log(Epluspz / Eminuspz);
         eta = 0.5 * log(Ppluspz / Pminuspz); 
      }	//	else
      
      theta = atan2(pt, pz);
      phi = TVector2::Phi_0_2pi(atan2(py, px));
   }
   
   
   void
   ParticleMC::ComputeEventDependentQuantities(
                                               EventMC& event
                                               ) {
      try {
         const TLorentzVector& hadron = event.GetTrack(1)->Get4Vector();
         const TLorentzVector& lepton = event.GetTrack(2)->Get4Vector();
         const TLorentzVector& boson = event.GetTrack(3)->Get4Vector();
         
         TLorentzVector pHadronBoson = Get4Vector();
         
         z = hadron.Dot(pHadronBoson) / hadron.Dot(boson);
         
         double sqrts = ::sqrt(4. * lepton.Energy() * hadron.Energy());
         xFeynman = 2. * pz / sqrts;
         
         // Not the most efficient thing recalculating this boost for
         // each particle...
         TLorentzRotation boostToHadronRest(-(hadron.BoostVector()));
         // Boost the particle to the proton rest frame and calculate its
         // pT with respect to the virtual photon:
         TLorentzVector bosonInRf =
            (TLorentzVector(boson) *= boostToHadronRest);
         pHadronBoson *= boostToHadronRest;
         
         thetaGamma =      pHadronBoson.Angle(bosonInRf.Vect());
         ptVsGamma =      pHadronBoson.Pt(bosonInRf.Vect());
         
         TLorentzVector leptonInPrf(lepton);
         leptonInPrf *= boostToHadronRest;
         
         phiPrf = computeHermesPhiH(pHadronBoson, leptonInPrf, bosonInRf);
         
         // Rotate the frame of reference so that the virtual photon momentum
         // defines the positive z direction.
         // gamma*-p coordinates are:
         // ( (q x e`) x q, q x e`, q )
         TVector3 gammaProtonZ = bosonInRf.Vect();
         TVector3 gammaProtonY = gammaProtonZ.Cross(leptonInPrf.Vect());
         TRotation rotation;
         // First argument defines the y direction, the second defines the
         // yz plane:
         rotation.SetYAxis(gammaProtonY,gammaProtonZ );
         // Inverse --> rotate coordinate system to give new object
         // coordinates, not rotate object in fixed coordinates
         rotation = rotation.Inverse();
         pHadronBoson = Get4Vector();
         TLorentzRotation boost = boostToHadronRest;
         boost.Transform(rotation);
         pHadronBoson *= boost;
         
         // Determine the PDG code of the parent particle, if the particle
         // has a parent and the parent is present in the particle array.
         // The index of the particles from the Monte Carlo runs from [1,N]
         // while the index in the array runs from [0,N-1], so subtract 1
         // from the parent index to find its position.
         if(event.GetNTracks() > unsigned(orig - 1)) {
            parentId = event.GetTrack(orig - 1)->Id();
         } // if
      } // try
      catch(std::exception& e) {
         std::cerr <<
         "Exception in Particle::ComputeEventDependentQuantities: " <<
         e.what() << std::endl;
      } // catch
   }
   
   
   TLorentzVector ParticleMC::Get4Vector() const {
      return TLorentzVector(px, py, pz, E);
   }
   
   
   const EventMC* ParticleMC::GetEvent() const {
      return dynamic_cast<EventMC*>(event.GetObject());
   }
   
   
   const ParticleMC* ParticleMC::GetChild(UShort_t u) const {
      
      // Look up this particle's child via the event
      // containing it and the index of the child in that event.
      
      if(not GetEvent()) return NULL;
      
      // index is in the range [1,N]
      unsigned idx = daughter + u;
      if(daughter < 1 or // If first daughter index = 0, it has no children
         u >= NChildren()) { // Insufficient children
         return NULL;
      } // if
      
      --idx; // Convert [1,N] --> [0,N)
      const ParticleMC* p(NULL);
      // Check this index is within the # of particles in the event
      if(idx < GetEvent()->GetNTracks()) {
         p = GetEvent()->GetTrack(idx);
      } // if
      return p;
   }
   
   
   const ParticleMC* ParticleMC::GetParent() const {
      
      // Look up this particle's parent via the event
      // containing it and the index of the parent in that event.
      
      const ParticleMC* p(NULL);
      
      if(GetEvent()) {
         if(GetEvent()->GetNTracks() >= GetParentIndex()) {
            p = GetEvent()->GetTrack(GetParentIndex() - 1);
         } // if
      } // if
      
      return p;
   }
   
   
   Bool_t ParticleMC::HasChild(Int_t pdg) const {
      for(UInt_t i(0); i < NChildren(); ++i) {
         if(not GetChild(i)) continue;
         if(pdg == GetChild(i)->Id()) return true;
      } // for
      return false;
   }
   
   
   TLorentzVector
   ParticleMC::Get4VectorInHadronBosonFrame() const {
      double p_ = ptVsGamma / ::sin(ptVsGamma);
      double e_ = ::sqrt(::pow(p_, 2.) + ::pow(m, 2.));
      double px_ = ptVsGamma * ::cos(phiPrf);
      double py_ = ptVsGamma * ::sin(phiPrf);
      double pz_ = ptVsGamma / ::tan(thetaGamma);
      return TLorentzVector(px_, py_, pz_, e_);
   }
   
   
   void ParticleMC::SetEvent(EventMC* e) {
      event = e;
   }

} // namespace erhic
