//
// Pythia6EventBuilder.cxx
//
// Created by TB on 1/17/12.
// Copyright 2012 BNL. All rights reserved.
//

#include <cmath>
#include <iostream>
#include <memory>

#include <TMCParticle.h>
#include <TObjArray.h>
#include <TProcessID.h>
#include <TPythia6.h>

#include "Kinematics.h"
#include "ParticleMC.h"
#include "Pythia6EventBuilder.h"
#include "Pythia6ParticleBuilder.h"

namespace erhic {
   
   
   Pythia6EventBuilder::Pythia6EventBuilder()
   { }
   
   
   Pythia6EventBuilder::~Pythia6EventBuilder()
   { }
   
   
   EventPythia* Pythia6EventBuilder::Create() {
      
      // Save current object count
      int objectNumber = TProcessID::GetObjectCount();
      
      TPythia6* pythia = TPythia6::Instance();
      
      // Import all particles (not just final-state)
      TObjArray* particles = pythia->ImportParticles("All");
      
      // Construct the EventPythia object from the current
      // state of TPythia6.
      std::auto_ptr<EventPythia> event(new EventPythia);
      
      // Extract the event-wise quantities:
      event->SetNucleon(pythia->GetMSTI(12));
      event->SetTargetParton(pythia->GetMSTI(16));
      event->SetBeamParton(pythia->GetMSTI(15));
      event->SetGenEvent(1);
      event->SetTargetPartonX(pythia->GetPARI(34));
      event->SetBeamPartonX(pythia->GetPARI(33));
      event->SetBeamPartonTheta(pythia->GetPARI(53));
      event->SetLeptonPhi(pythia->GetVINT(313));
      event->SetHardS(pythia->GetPARI(14));
      event->SetHardT(pythia->GetPARI(15));
      event->SetHardU(pythia->GetPARI(16));
      event->SetHardQ2(pythia->GetPARI(22));
      event->SetHardPt2(pythia->GetPARI(18));
      event->SetPhotonFlux(pythia->GetVINT(319));
      event->SetProcess(pythia->GetMSTI(1));
      
      // We need the beam energies to compute the true x, W2 and nu.
      // The beam lepton should be the first particle
      // and the beam hadron the second particle.
      const double eLepton =
         dynamic_cast<TMCParticle*>(particles->At(0))->GetEnergy();
      const double eHadron =
         dynamic_cast<TMCParticle*>(particles->At(1))->GetEnergy();
      const double mHadron =
         dynamic_cast<TMCParticle*>(particles->At(1))->GetMass();
      
      // x, W2 and nu are calculated from y, Q2 and the beam energies.
      // y and Q2 come from PYTHIA.
      // Use (approximate expression) Q2 = sxy, where s = 4.E_e.E_h
      
      double y = pythia->GetVINT(309);
      double Q2 = pythia->GetVINT(307);
      double x = Q2 / y / 4. / eLepton / eHadron;
      double W2 = pow(mHadron, 2.) + Q2 * (1. / x - 1.);
      double nu = (W2 + Q2 - pow(mHadron, 2.)) / 2. / mHadron;
      
      event->SetTrueY(y);
      event->SetTrueQ2(Q2);
      event->SetTrueX(x);
      event->SetTrueW2(W2);
      event->SetTrueNu(nu);
      
      // Now populate the particle list.

      Pythia6ParticleBuilder builder;
      for(int i(0); i < particles->GetEntries(); ++i) {
         TMCParticle* p =
            dynamic_cast<TMCParticle*>(particles->At(i));
         std::auto_ptr<ParticleMC> particle = builder.Create(*p);
         particle->SetIndex(i + 1);
         particle->SetEvent(event.get());
         event->AddLast(particle.release());
      } // for

      // Compute derived event kinematics
      DisKinematics* nm = LeptonKinematicsComputer(*event).Calculate();
      DisKinematics* jb = JacquetBlondelComputer(*event, NULL).Calculate();
      DisKinematics* da = DoubleAngleComputer(*event, NULL).Calculate();
      if(nm) {
         event->SetLeptonKinematics(*nm);
      } // if
      if(jb) {
         event->SetJacquetBlondelKinematics(*jb);
      } // if
      if(da) {
         event->SetDoubleAngleKinematics(*da);
      } // if
      // We also have to set the remaining variables not taken care of
      // by the general DIS event kinematic computations.
      // Find the beams, exchange boson, scattered lepton.
      BeamParticles beams;
      if(ParticleIdentifier::IdentifyBeams(*event, beams)) {
         const TLorentzVector h = beams.BeamHadron();
         TLorentzVector l = beams.BeamLepton();
         TLorentzVector s = beams.ScatteredLepton();
         TVector3 boost = -h.BoostVector();
         l.Boost(boost);
         s.Boost(boost);
         event->SetELeptonInNuclearFrame(l.E());
         event->SetEScatteredInNuclearFrame(s.E());
      } // if

      //Restore Object count 
      // See example in $ROOTSYS/test/Event.cxx
      //To save space in the table keeping track of all referenced objects
      //we assume that our events do not address each other. We reset the 
      //object count to what it was at the beginning of the event.
      TProcessID::SetObjectCount(objectNumber);
      
      
      return event.release();
   }
   
} // namespace erhic
