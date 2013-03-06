#ifndef _EICSMEAR_ParticleIdentifier_
#define _EICSMEAR_ParticleIdentifier_

#include <cmath>
#include <limits>

#include "eicsmear/erhic/BeamParticles.h"
#include "eicsmear/erhic/Particle.h"
#include "eicsmear/erhic/VirtualEvent.h"

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
   ParticleIdentifier(const int leptonPdg = ~unsigned(0)/2);
   
   virtual ~ParticleIdentifier() { }
   
   /**
    Returns whether the particle is the beam lepton.
    */
   virtual bool isBeamLepton(const erhic::VirtualParticle& ) const;
   
   /**
    Returns whether the particle is the beam hadron.
    */
   virtual bool isBeamNucleon(const erhic::VirtualParticle& ) const;
   
   /**
    Returns whether the particle is the scattered lepton beam particle.
    */
   virtual bool isScatteredLepton(const erhic::VirtualParticle& ) const;
   
   /**
    Returns whether the particle is a virtual photon.
    */
   virtual bool IsVirtualPhoton(const erhic::VirtualParticle& ) const;
   
   /**
    Returns whether the particles should be skipped by the tree building code.
    */
   virtual bool SkipParticle(const erhic::VirtualParticle& ) const;
   
   /**
    Sets the PDG code to use when identifying the lepton beam.
    */
   virtual void SetLeptonBeamPdgCode(int );
   
   /**
    Returns the PDG code to use when identifying the lepton beam.
    */
   virtual int GetLeptonBeamPdgCode() const;
   
   /**
    Look for charged current events.
    
    In this case, the scattered lepton searched for will be the neutrino
    corresponding to the incident lepton beam type (e.g. W- for electron,
    W+ for proton).
    */
   virtual bool SetChargedCurrent(bool);
   
   /**
    Identify the beams from an event and store their properties in a
    BeamParticles object.
    See BeamParticles.h for the quantities stored.
    Returns true if all beams are found, false if not.
    Important: finding the scattered hadron beam is not implemented.
    */
//   
//   static bool IdentifyBeams(const T&, BeamParticles& );
//   static bool IdentifyBeams(const EventBase&, BeamParticles& );
   static bool IdentifyBeams(const erhic::VirtualEvent&, BeamParticles&);
   
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
//   
//   static bool IdentifyBeams(const T&, std::vector<const Particle*>& );
//   static bool IdentifyBeams(const EventBase&, std::vector<const Particle*>& );
   static bool IdentifyBeams(const erhic::VirtualEvent&,
                             std::vector<const erhic::VirtualParticle*>&);
   
protected:
   
   /**
    Determine the scattered lepton type from an incident lepton type.
    
    For neutral current events these are equal.
    For charged current events the scattered lepton will be a neutrino, the
    type of which depends on the incident lepton type e.g.
      beam electron (11) --> scattered nu_e (12)
      beam positron (-11) --> scattered nu_e_bar (-12)
    Even though we would not envisage mu or tau beams (!) the algorithm is
    still robust for these inputs.
    */
   Int_t DetermineScatteredType(Int_t);
   
   Bool_t mChargedCurrent;
   Int_t mLeptonBeamPdgCode;
   Int_t mScatteredPdgCode;
};

inline int ParticleIdentifier::GetLeptonBeamPdgCode() const {
   return mLeptonBeamPdgCode;
}

#endif

// 2011.07.29: Fixed bug in SkipParticle() causing particles with
//    id<10 to be skipped instead of abs(id)<10.
