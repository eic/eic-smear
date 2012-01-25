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
//      std::cout << "owned? " << particles->IsOwner() << std::endl;
      
      // Construct the EventPythia object from the current
      // state of TPythia6.
      /**
       \todo Fill genevent in Pythia6EventBuilder::Create()
       \todo Compute & store F1, F2, R, sigma_rad, SigRadCor, EBrems
       */
      std::auto_ptr<EventPythia> event(new EventPythia);
//      EventPythia* event = new EventPythia;
      
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
//         particle->SetEvent(event);
         event->AddLast(particle.release());
      } // for
      // We need to calculate any remaining quantities that depend
      // on the particle list.
      
      event->Compute();
      
      
      // FORTRAN NGEN(0, 3) corresponds to NGEN[2][0] in the C version
      // i.e. reversal of indices and 3 --> 2 because of differing
      // C and FORTRAN array indexing.
      //      std::cout << pythia->GetPyint5()->NGEN[2][0] << " trials" << std::endl;
      
      //Restore Object count 
      // See example in $ROOTSYS/test/Event.cxx
      //To save space in the table keeping track of all referenced objects
      //we assume that our events do not address each other. We reset the 
      //object count to what it was at the beginning of the event.
      TProcessID::SetObjectCount(objectNumber);
      
      
      return event.release();
   }
   
} // namespace erhic
