//
// ParticleIdentifier.cxx
//
// Created by TB on 6/10/11.
// Copyright 2011 BNL. All rights reserved.
//

#include <iostream>

#include "eicsmear/erhic/EventMC.h"
#include "eicsmear/erhic/ParticleIdentifier.h"

//ClassImp(ParticleIdentifier<erhic::EventMC>)


//	================================================================================================
// Constructor
//	================================================================================================

ParticleIdentifier::ParticleIdentifier(
                                       const int leptonBeamPdgCode
                                       )
: mLeptonBeamPdgCode(leptonBeamPdgCode )
{ }

//	================================================================================================
// Identify the lepton beam particle
//	================================================================================================

bool
ParticleIdentifier::isBeamLepton(const erhic::VirtualParticle& particle )
const
{
   // Test against status 201 for SOPHIA events
   return
   (21 == particle.GetStatus() or 201 == particle.GetStatus()) and
   GetLeptonBeamPdgCode() == particle.Id() and
   //   particle.I < 3;
   particle.GetParentIndex() == 0;
}

//	================================================================================================
// Identify the scattered lepton
//	================================================================================================

bool
ParticleIdentifier::isScatteredLepton(const erhic::VirtualParticle& particle )
const
{
   return 1 == particle.GetStatus() and GetLeptonBeamPdgCode() == particle.Id();
}

//	=================================================================================================
// Return true if the particle matches any of the following:
//    Outgoing electron with KS 21 - we instead pick up the repeat entry, when KS < 21.
//    Repeat of incoming nucleon.
//    Intermediate gluon.
//    Intermediate quark.
//	=================================================================================================

bool
ParticleIdentifier::SkipParticle(const erhic::VirtualParticle& particle )
const
{
   int kI1 = particle.GetStatus();
   int pdgCode = particle.Id();
   int parent = particle.GetParentIndex();
   // Remove duplicate outgoing electron, info is picked up later when KS < 21
   if(21 == kI1 and pdgCode == GetLeptonBeamPdgCode() and parent == 1 ) {
      return true;;
   }	//	if
   
   // Remove repeat of incoming nucelon
   if(21 == kI1 and (pdgCode == 2112 or pdgCode == 2212) and parent == 2 ) {
      return true;;
   }	//	if
   
   // Remove intermediate gluons
   if(21 == kI1 and pdgCode == 21 ) {
      return true;
   }	//	if
   
   // Remove intermediate (anti)quarks
   if(21 == kI1 and ::abs(pdgCode) < 10 ) {
      return true;;
   }	//	if
   
   return false;
}

//	================================================================================================
// Identify the virtual photon based on its PDG code and PYTHIA K(I,1)
//	================================================================================================

bool
ParticleIdentifier::IsVirtualPhoton(const erhic::VirtualParticle& particle )
const
{
   return (22 == particle.Id() or 23 == particle.Id() ) and
   21 == particle.GetStatus();
}

//	================================================================================================
// Identify the nucleon beam particle
//	================================================================================================

bool
ParticleIdentifier::isBeamNucleon(const erhic::VirtualParticle& particle )
const
{
   // Test against status 201 for SOPHIA events
   return (21  == particle.GetStatus() or 201 == particle.GetStatus()) and
   (2112 == particle.Id() or 2212 == particle.Id()) and
   particle.GetParentIndex() == 0;
}


/**
 Identify the beams from an event and store their properties in a
 BeamParticles object.
 See BeamParticles.h for the quantities stored.
 Returns true if all beams are found, false if not.
 */

bool
ParticleIdentifier::IdentifyBeams(
                                  const erhic::VirtualEvent& event,
                                  BeamParticles& beams
                                  ) {
   beams.Reset();
   
   const unsigned nParticles = event.GetNTracks();
   
   // Count leptons so we don't overwrite the beam and scattered lepton with
   // subsequent leptons of the same type.
   int leptonCount(0);
   
   // Set to true once we find the first virtual photon, so we can skip
   // subsequent virtual photons.
   bool haveFoundVirtualPhoton(false );
   
   ParticleIdentifier finder;
   
   for(unsigned n(0); n < nParticles; ++n ) {
      
      const erhic::VirtualParticle* ptr = event.GetTrack(n);
      if(not ptr) {
         continue;
      } // if
      
      const erhic::VirtualParticle& particle = *ptr;
      
      // The first particle is the beam lepton, so set this value in the
      // ParticleIdentifier
      if(0 == n ) {
         finder.SetLeptonBeamPdgCode(particle.Id() );
      } // if
      
      if(finder.isBeamNucleon(particle ) ) {
         beams.SetBeamHadron(particle.Get4Vector() );
      }	//	if...
      else if(finder.isBeamLepton(particle ) and 0 == leptonCount ) {
         beams.SetBeamLepton(particle.Get4Vector() );
         ++leptonCount;
      }	//	else if...
      else if(finder.isScatteredLepton(particle ) and 1 == leptonCount ) {
         // Protect against additional KS == 1 electrons following this
         beams.SetScatteredLepton(particle.Get4Vector() );
         // Protect against additional KS == 1 electrons following this
         ++leptonCount;
      }	//	if
      else if(finder.IsVirtualPhoton(particle ) and not haveFoundVirtualPhoton ) {
         beams.SetBoson(particle.Get4Vector() );
         haveFoundVirtualPhoton = true;
      }	//	if
   } // for
   
   if(isnan(beams.BeamHadron().E() ) or
      isnan(beams.BeamLepton().E() ) or
      isnan(beams.ScatteredLepton().E())) {
//      isnan(beams.GetBoson().E() ) ) {
      return false;
   } // if
   
   return true;
}


bool
ParticleIdentifier::IdentifyBeams(
                                  const erhic::VirtualEvent& event,
                                  std::vector<const erhic::VirtualParticle*>& beams
                                  ) {
   //   std::cout << "ParticleIdentifier::IdentifyBeams"<<std::endl;
   const Particle* const null(NULL);
   
   beams.assign(4, null);
   
   const unsigned nParticles = event.GetNTracks();
   //   std::cout<<nParticles<<" particles"<<std::endl;
   
   // Count leptons so we don't overwrite the beam and scattered lepton with
   // subsequent leptons of the same type.
   int leptonCount(0);
   
   // Set to true once we find the first virtual photon, so we can skip
   // subsequent virtual photons.
   bool haveFoundVirtualPhoton(false );
   
   ParticleIdentifier finder;
   
   for(unsigned n(0); n < nParticles; ++n ) {
      //      std::cout<<"particle "<<n<<std::endl;
      const erhic::VirtualParticle* ptr = event.GetTrack(n);
      if(not ptr) {
         continue;
      } // if
      
      const erhic::VirtualParticle& particle = *ptr;
      //      particle.Dump();
      // The first particle is the beam lepton, so set this value in the
      // ParticleIdentifier
      if(0 == n ) {
         finder.SetLeptonBeamPdgCode(particle.Id() );
         //         std::cout<<"lepton pdg "<<particle.Id()<<std::endl;
      } // if
      
      if(finder.isBeamNucleon(particle) ) {
         beams.at(1) = &particle;
         //         std::cout<<"hadron beam"<<std::endl;
      }	//	if...
      else if(finder.isBeamLepton(particle) and 0 == leptonCount ) {
         beams.at(0) = &particle;
         ++leptonCount;
         //         std::cout <<"lepton beam"<<std::endl;
      }	//	else if...
      else if(finder.isScatteredLepton(particle) and 1 == leptonCount ) {
         // Protect against additional KS == 1 electrons following this
         beams.at(3) = &particle;
         //         std::cout<<"scattered"<<std::endl;
         // Protect against additional KS == 1 electrons following this
         ++leptonCount;
      }	//	if
      else if(finder.IsVirtualPhoton(particle ) and not haveFoundVirtualPhoton ) {
         beams.at(2) = &particle;
         //         std::cout<<"virtual photon"<<std::endl;
         haveFoundVirtualPhoton = true;
      }	//	if
   } // for
   
   return std::find(beams.begin(), beams.end(), null) == beams.end();
}
