/**
 \file
 Declaration of class Smear::EventSmear.
 
 \author    Michael Savastio
 \date      2011-10-10
 \copyright 2011 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_SMEAR_EVENTSMEAR_H_
#define INCLUDE_EICSMEAR_SMEAR_EVENTSMEAR_H_

#include <cmath>
#include <list>
#include <vector>

#include <TObject.h>

#include "eicsmear/erhic/EventDis.h"
#include "eicsmear/erhic/Kinematics.h"
#include "eicsmear/smear/ParticleMCS.h"
#include "eicsmear/erhic/VirtualParticle.h"

namespace Smear {

/*
 A generator-independent DIS event with smeared kinematics and particles.
 */
class Event : public erhic::EventDis {
 public:
  /**
   Default constructor.
   */
  Event();

  /**
   Destructor.
   */
  virtual ~Event();

  /**
   Clear the particle list, sets event properties to default values.
   */
  virtual void Reset();

  /**
   Clears particle array, leaves event variables unchanged.
   */
  virtual void ClearParticles();

  /**
   Returns the number of tracks in the event.
   */
  virtual UInt_t GetNTracks() const;

  /**
   Returns the nth track.
   Returns NULL if the track number is out of the range [0, GetNTracks()).
   @param [in] The track index, in the range [0, GetNTracks()).
   */
  virtual const ParticleMCS* GetTrack(UInt_t) const;

  /**
   Returns the nth track.
   Returns NULL if the track number is out of the range [0, GetNTracks()).
   @param [in] The track index, in the range [0, GetNTracks()).
   */
  virtual ParticleMCS* GetTrack(UInt_t);

  virtual void SetQ2(double Q2) { QSquared = Q2; }

  virtual void SetX(double xB) { x = xB; }

  virtual void SetY(double inelasticity) { y = inelasticity; }

  virtual void SetW2(double W2) { WSquared = W2; }

  virtual void SetNu(double Nu) { nu = Nu; }

  /**
   Returns a pointer to the incident lepton beam particle.
   Returns a NULL pointer if the particle cannot be located in the event.
   IMPORTANT - DO NOT DELETE THE POINTER OR BAD THINGS WILL HAPPEN!
   
   In the standard eRHIC Monte Carlo format, the incident lepton beam
   is assumed to be the first particle in the particle list.
   This is the behaviour implemented here.
   Derived classes can implement other selection mechanisms depending on
   their data format.
   */
  virtual const ParticleMCS* BeamLepton() const;

  /**
   Returns a pointer to the incident hadron beam particle.
   See also notes in BeamLepton().
   
   In the standard eRHIC Monte Carlo format, the incident hadron beam
   is assumed to be the second particle in the particle list.
   */
  virtual const ParticleMCS* BeamHadron() const;

  /**
   Returns a pointer to the exchanged boson.
   See also notes in BeamLepton().
   
   In the standard eRHIC Monte Carlo format, the exchanged boson
   is assumed to be the third particle in the particle list.
   */
  virtual const ParticleMCS* ExchangeBoson() const;

  /**
   Returns a pointer to the lepton beam particle after scattering.
   See also notes in BeamLepton().
   
   In the standard eRHIC Monte Carlo format, the scattered lepton beam
   is assumed to be the first final-state particle in the particle list
   with the same PDG code as the incident lepton beam.
   */
  virtual const ParticleMCS* ScatteredLepton() const;

  /**
   Add a new track to the end of the track list.
   The track must be allocated via new and is subsequently owned
   by the Event.
   */
  virtual void AddLast(ParticleMCS* particle);

  /**
   Yields all particles that belong to the hadronic final state.
   This is the same as the result of FinalState(), minus the scattered
   beam lepton.
   */
  void HadronicFinalState(ParticlePtrList&) const;

  /**
   Returns a vector of pointers to all tracks in the event.
   Note that this includes NULL pointers to tracks that were not detected.
   Do not delete the pointers.
   */
  std::vector<const erhic::VirtualParticle*> GetTracks() const;

  /**
   Set which particle is the scattered lepton.
   */
  virtual void SetScattered(int index);

  /**
   Prints the attributes of this event to standard output.
   Prints event-wise kinematic values, and all tracks in the event.
   */
  virtual void Print(Option_t* = "") const;

 protected:
  Int_t nTracks;  ///< Number of particles (intermediate + final)
  std::vector<ParticleMCS*> particles;  ///< The smeared particle list
  Int_t mScatteredIndex;

  ClassDef(Smear::Event, 1)
};

inline UInt_t Event::GetNTracks() const {
  return particles.size();
}

inline const Smear::ParticleMCS* Event::GetTrack(UInt_t u) const {
  return (u < particles.size() ? particles.at(u) : NULL);
}

inline Smear::ParticleMCS* Event::GetTrack(UInt_t u) {
  return (u < particles.size() ? particles.at(u) : NULL);
}

inline const ParticleMCS* Event::BeamLepton() const {
  return (particles.empty() ? NULL : particles.front());
}

inline const ParticleMCS* Event::BeamHadron() const {
  return (particles.size() > 1 ? particles.at(1) : NULL);
}

inline const ParticleMCS* Event::ExchangeBoson() const {
  return NULL;
}

}  // namespace Smear

typedef Smear::Event EventS;

#endif  // INCLUDE_EICSMEAR_SMEAR_EVENTSMEAR_H_
