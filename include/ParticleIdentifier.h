#ifndef _ERHIC_BUILDTREE_PARTICLEIDENTIFIER_
#define _ERHIC_BUILDTREE_PARTICLEIDENTIFIER_

#include <cmath>
#include <limits>

#include "BeamParticles.h"
#include "Particle.h"
#include "VirtualEvent.h"

/**
 Implements methods to identify particles based on their species and status
 codes.
 */
struct ParticleIdentifier {
   
   /**
    Default constructor.
    Initialise with the PDG code of the lepton beam.
    The default is an invalid value.
    */
//   ParticleIdentifier(const int leptonPdg = std::numeric_limits<int>::max() );
   ParticleIdentifier(const int leptonPdg = ~unsigned(0)/2 );
   
   virtual ~ParticleIdentifier() { }
   
   /**
    Returns whether the particle is the beam lepton.
    */
   virtual bool isBeamLepton(const Particle& ) const;
   
   /**
    Returns whether the particle is the beam hadron.
    */
   virtual bool isBeamNucleon(const Particle& ) const;
   
   /**
    Returns whether the particle is the scattered lepton beam particle.
    */
   virtual bool isScatteredLepton(const Particle& ) const;
   
   /**
    Returns whether the particle is a virtual photon.
    */
   virtual bool IsVirtualPhoton(const Particle& ) const;
   
   /**
    Returns whether the particles should be skipped by the tree building code.
    */
   virtual bool SkipParticle(const Particle& ) const;
   
   /**
    Sets the PDG code to use when identifying the lepton beam.
    */
   virtual void SetLeptonBeamPdgCode(int );
   
   /**
    Returns the PDG code to use when identifying the lepton beam.
    */
   virtual int GetLeptonBeamPdgCode() const;
   
   /**
    Identify the beams from an event and store their properties in a
    BeamParticles object.
    See BeamParticles.h for the quantities stored.
    Returns true if all beams are found, false if not.
    Important: finding the scattered hadron beam is not implemented.
    */
//   template<typename T>
//   static bool IdentifyBeams(const T&, BeamParticles& );
//   static bool IdentifyBeams(const EventBase&, BeamParticles& );
   static bool IdentifyBeams(const erhic::EventMC&, BeamParticles&);
   
   /**
    Identify the beams from an event and store their properties in a
    vector of pointers to the particle objects in the event.
    Do not delete the pointers, as they belong to the event.
    The returned vector has four entries, in this order: incident lepton,
    incident hadron, exchanged boson, scattered lepton.
    Any particle not found yields a NULL pointer in the vector.
    Returns true if all beams are found (i.e. no NULL pointers), false if not.
    Important: finding the scattered hadron beam is not implemented.
    */
//   template<typename T>
//   static bool IdentifyBeams(const T&, std::vector<const Particle*>& );
//   static bool IdentifyBeams(const EventBase&, std::vector<const Particle*>& );
   static bool IdentifyBeams(const erhic::EventMC&, std::vector<const Particle*>&);
   
protected:
   
   Int_t mLeptonBeamPdgCode;
};

inline void ParticleIdentifier::SetLeptonBeamPdgCode(const int pdgCode ) {
   mLeptonBeamPdgCode = pdgCode;
}

inline int ParticleIdentifier::GetLeptonBeamPdgCode() const {
   return mLeptonBeamPdgCode;
}

#endif

// 2011.07.29: Fixed bug in SkipParticle() causing particles with
//    id<10 to be skipped instead of abs(id)<10.
