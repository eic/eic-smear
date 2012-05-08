#ifndef _ERHIC_BUILDTREE_BEAMPARTICLES_
#define _ERHIC_BUILDTREE_BEAMPARTICLES_

//
// BeamParticles.h
//
// Created by TB on 6/24/11.
// Copyright 2011 BNL. All rights reserved.
//

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
   
   BeamParticles(const TLorentzVector& hadronBeam,
                  const TLorentzVector& leptonBeam,
                  const TLorentzVector& scatteredHadron,
                  const TLorentzVector& scatteredLepton,
                  const TLorentzVector& exchangedBoson );
   
   virtual ~BeamParticles();
   
   /**
    Sets all the 4-vectors' components to not-a-number.
    */
   void Reset();
   
   void SetBeamHadron(const TLorentzVector& );
   void SetBeamLepton(const TLorentzVector& );
   void SetScatteredHadron(const TLorentzVector& );
   void SetScatteredLepton(const TLorentzVector& );
   void SetBoson(const TLorentzVector& );
   
   const TLorentzVector& BeamHadron() const;
   const TLorentzVector& BeamLepton() const;
   const TLorentzVector& GetScatteredHadron() const;
   const TLorentzVector& ScatteredLepton() const;
   const TLorentzVector& GetBoson() const;
   
protected:
   
   TLorentzVector mBeamHadron;      ///< Incident hadron beam
   TLorentzVector mBeamLepton;      ///< Incident lepton beam
   TLorentzVector mScatteredHadron; ///< Scattered hadron beam
   TLorentzVector mScatteredLepton; ///< Scattered lepton beam
   TLorentzVector mBoson;           ///< Exchanged boson
   
   ClassDef(BeamParticles, 1 )
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

inline void BeamParticles::SetBeamHadron(const TLorentzVector& vec ) {
   mBeamHadron = vec;
}

inline void BeamParticles::SetBeamLepton(const TLorentzVector& vec ) {
   mBeamLepton = vec;
}

inline void BeamParticles::SetScatteredHadron(const TLorentzVector& vec ) {
   mScatteredHadron = vec;
}

inline void BeamParticles::SetScatteredLepton(const TLorentzVector& vec ) {
   mScatteredLepton = vec;
}

inline void BeamParticles::SetBoson(const TLorentzVector& vec ) {
   mBoson = vec;
}

#endif
