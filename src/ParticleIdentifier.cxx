//
// ParticleIdentifier.cxx
//
// Created by TB on 6/10/11.
// Copyright 2011 BNL. All rights reserved.
//

#include <iostream>

#include "ParticleIdentifier.h"

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
ParticleIdentifier::isBeamLepton(const Particle& particle )
const
{
   // Test against status 201 for SOPHIA events
   return
   (21 == particle.KS or 201 == particle.KS) and
   GetLeptonBeamPdgCode() == particle.id and
   //   particle.I < 3;
   particle.orig == 0;
}

//	================================================================================================
// Identify the scattered lepton
//	================================================================================================
bool
ParticleIdentifier::isScatteredLepton(const Particle& particle )
const
{
   return 1 == particle.Status() and GetLeptonBeamPdgCode() == particle.Id();
}

//	=================================================================================================
// Return true if the particle matches any of the following:
//    Outgoing electron with KS 21 - we instead pick up the repeat entry, when KS < 21.
//    Repeat of incoming nucleon.
//    Intermediate gluon.
//    Intermediate quark.
//	=================================================================================================
bool
ParticleIdentifier::SkipParticle(const Particle& particle )
const
{
   int kI1 = particle.Status();
   int pdgCode = particle.Id();
   int parent = particle.ParentIndex();
   // Remove duplicate outgoing electron, info is picked up later when KS < 21
   if(21 == kI1 and pdgCode == GetLeptonBeamPdgCode() and parent == 1 ) {
      return true;;
   }	//	if
   
   // Remove repeat of incoming nucelon
   if(21 == kI1 and pdgCode == 2212 and parent == 2 ) {
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
ParticleIdentifier::IsVirtualPhoton(const Particle& particle )
const
{
   return (22 == particle.Id() or 23 == particle.Id() ) and
   21 == particle.Status();
}

//	================================================================================================
// Identify the nucleon beam particle
//	================================================================================================
bool
ParticleIdentifier::isBeamNucleon(const Particle& particle )
const
{
   // Test against status 201 for SOPHIA events
   return (21  == particle.KS or 201 == particle.KS) and
   2212 == particle.id and// particle.I < 3;
   particle.orig == 0;
}


/**
 Identify the beams from an event and store their properties in a
 BeamParticles object.
 See BeamParticles.h for the quantities stored.
 Returns true if all beams are found, false if not.
 */
bool
ParticleIdentifier::IdentifyBeams(
                                  const erhic::EventMC& event,
                                  BeamParticles& beams
                                  ) {
   beams.Reset();
   
   const unsigned nParticles = event.NTracks();
   
   // Count leptons so we don't overwrite the beam and scattered lepton with
   // subsequent leptons of the same type.
   int leptonCount(0);
   
   // Set to true once we find the first virtual photon, so we can skip
   // subsequent virtual photons.
   bool haveFoundVirtualPhoton(false );
   
   ParticleIdentifier finder;
   
   for(unsigned n(0); n < nParticles; ++n ) {
      
      const Particle* ptr = event.GetTrack(n);
      if(not ptr) {
         continue;
      } // if
      
      const Particle& particle = *ptr;
      
      // The first particle is the beam lepton, so set this value in the
      // ParticleIdentifier
      if(0 == n ) {
         finder.SetLeptonBeamPdgCode(particle.Id() );
      } // if
      
      if(finder.isBeamNucleon(particle ) ) {
         beams.SetBeamHadron(particle.PxPyPzE() );
      }	//	if...
      else if(finder.isBeamLepton(particle ) and 0 == leptonCount ) {
         beams.SetBeamLepton(particle.PxPyPzE() );
         ++leptonCount;
      }	//	else if...
      else if(finder.isScatteredLepton(particle ) and 1 == leptonCount ) {
         // Protect against additional KS == 1 electrons following this
         beams.SetScatteredLepton(particle.PxPyPzE() );
         // Protect against additional KS == 1 electrons following this
         ++leptonCount;
      }	//	if
      else if(finder.IsVirtualPhoton(particle ) and not haveFoundVirtualPhoton ) {
         beams.SetBoson(particle.PxPyPzE() );
         haveFoundVirtualPhoton = true;
      }	//	if
   } // for
   
   if(isnan(beams.BeamHadron().E() ) or
      isnan(beams.BeamLepton().E() ) or
      isnan(beams.ScatteredLepton().E() ) or
      isnan(beams.GetBoson().E() ) ) {
      return false;
   } // if
   
   return true;
}

bool
ParticleIdentifier::IdentifyBeams(
                                  const erhic::EventMC& event,
                                  std::vector<const Particle*>& beams
                                  ) {
   //   std::cout << "ParticleIdentifier::IdentifyBeams"<<std::endl;
   const Particle* const null(NULL);
   
   beams.assign(4, null);
   
   const unsigned nParticles = event.NTracks();
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
      const Particle* ptr = event.GetTrack(n);
      if(not ptr) {
         continue;
      } // if
      
      const Particle& particle = *ptr;
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
