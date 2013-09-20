/**
 \file
 Declaration of class BeamParticles.
 
 \author    Thomas Burton
 \date      2011-06-24
 \copyright 2011 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_ERHIC_BEAMPARTICLES_H_
#define INCLUDE_EICSMEAR_ERHIC_BEAMPARTICLES_H_

#include <Rtypes.h>
#include <TLorentzVector.h>

/**
 Wrapper class around energy-momentum 4-vectors defining the incident and
 scattered beams and the exchanged boson.
 */
class BeamParticles {
 public:
  /**
   Default constructor.
   Initialises all vector components to not-a-number.
   */
  BeamParticles();

  /**
   Constructor initialsing all particles.
   */
  BeamParticles(const TLorentzVector& hadronBeam,
                const TLorentzVector& leptonBeam,
                const TLorentzVector& scatteredHadron,
                const TLorentzVector& scatteredLepton,
                const TLorentzVector& exchangedBoson);

  /**
   Destructor.
   */
  virtual ~BeamParticles();

  /**
   Sets all the 4-vectors' components to not-a-number.
   */
  void Reset();

  void SetBeamHadron(const TLorentzVector&);

  void SetBeamLepton(const TLorentzVector&);

  void SetScatteredHadron(const TLorentzVector&);

  void SetScatteredLepton(const TLorentzVector&);

  void SetBoson(const TLorentzVector&);

  const TLorentzVector& BeamHadron() const;

  const TLorentzVector& BeamLepton() const;

  const TLorentzVector& GetScatteredHadron() const;

  const TLorentzVector& ScatteredLepton() const;

  const TLorentzVector& GetBoson() const;

 protected:
  TLorentzVector mBeamHadron;       ///< Incident hadron beam
  TLorentzVector mBeamLepton;       ///< Incident lepton beam
  TLorentzVector mScatteredHadron;  ///< Scattered hadron beam
  TLorentzVector mScatteredLepton;  ///< Scattered lepton beam
  TLorentzVector mBoson;            ///< Exchanged boson

  ClassDef(BeamParticles, 1)
};

inline const TLorentzVector& BeamParticles::BeamHadron() const {
  return mBeamHadron;
}

inline const TLorentzVector& BeamParticles::BeamLepton() const {
  return mBeamLepton;
}

inline const TLorentzVector& BeamParticles::GetScatteredHadron() const {
  return mScatteredHadron;
}

inline const TLorentzVector& BeamParticles::ScatteredLepton() const {
  return mScatteredLepton;
}

inline const TLorentzVector& BeamParticles::GetBoson() const {
  return mBoson;
}

inline void BeamParticles::SetBeamHadron(const TLorentzVector& vec) {
  mBeamHadron = vec;
}

inline void BeamParticles::SetBeamLepton(const TLorentzVector& vec) {
  mBeamLepton = vec;
}

inline void BeamParticles::SetScatteredHadron(const TLorentzVector& vec) {
  mScatteredHadron = vec;
}

inline void BeamParticles::SetScatteredLepton(const TLorentzVector& vec) {
  mScatteredLepton = vec;
}

inline void BeamParticles::SetBoson(const TLorentzVector& vec) {
  mBoson = vec;
}

#endif  // INCLUDE_EICSMEAR_ERHIC_BEAMPARTICLES_H_
