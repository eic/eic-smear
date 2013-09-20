/**
 \file
 Declaration of class erhic::hadronic::ParticleMC.
 
 \author    Thomas Burton 
 \date      2012-05-02
 \copyright 2012 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_HADRONIC_PARTICLEMC_H_
#define INCLUDE_EICSMEAR_HADRONIC_PARTICLEMC_H_

#include <iostream>

#include <TLorentzVector.h>
#include <TVector3.h>

#include "eicsmear/erhic/Pid.h"
#include "eicsmear/erhic/VirtualParticle.h"

class TMCParticle;  // ROOT Monte Carlo particle class

namespace erhic {
namespace hadronic {

/**
 A realisation of erhic::VirtualParticle for tracks from a
 hadron-hadron Monte Carlo event.
 */
class ParticleMC : public erhic::VirtualParticle {
 public:
  /**
   Destructor
   */
  virtual ~ParticleMC() { }

  /**
   Default constructor
   */
  ParticleMC();

  /**
   Initialise from a PYTHIA TParticleMC
   */
  explicit ParticleMC(const TMCParticle&);

  /**
   Initialise from energy-momentum 4-vector, vertex 3-vector and PDG,
   status, parent index.
   */
  ParticleMC(const TLorentzVector&, const TVector3&, int, int, int);

  /**
   Returns identity information for the Particle species.
   */
  virtual erhic::Pid Id() const;

  /**
   Returns the momentum-energy four-vector (px, py, pz, E).
   */
  virtual TLorentzVector Get4Vector() const;

  /**
   Returns the x component of 3-momentum.
   */
  virtual Double_t GetPx() const;

  /**
   Returns the y component of 3-momentum.
   */
  virtual Double_t GetPy() const;

  /**
   Returns the z component of 3-momentum.
   */
  virtual Double_t GetPz() const;

  /**
   Returns total energy.
   */
  virtual Double_t GetE() const;

  /**
   Returns the magnitude of 3-momentum (GeV).
   */
  virtual Double_t GetP() const;

  /**
   Returns invariant mass (GeV/c<sup>2</sup>).
   */
  virtual Double_t GetM() const;

  /**
   Returns momentum perpendicular to the beam direction.
   */
  virtual Double_t GetPt() const;

  /**
   Returns the polar angle in the range [0, pi] radians.
   */
  virtual Double_t GetTheta() const;

  /**
   Returns the polar angle in the range [0, 2pi] radians.
   */
  virtual Double_t GetPhi() const;

  /**
   Returns the rapidity.
   */
  virtual Double_t GetRapidity() const;

  /**
   Returns the pseudorapidity.
   */
  virtual Double_t GetEta() const;

  /**
   Returns the origin point of the particle in cm.
   (0,0,0) indicates a particle originating in the collision.
   */
  virtual TVector3 GetVertex() const;

  /**
   A general "status" code for the particle
   (definition depends on implementation).
   */
  virtual UShort_t GetStatus() const;

  /**
   Index of this particle's precursor in the event.
   Returns 0 if the particle has no direct parent.
   */
  virtual UShort_t GetParentIndex() const;

  /**
   Returns Feynman-x.
   x<sub>F</sub> = 2*p<sub>z</sub>/sqrt(s).
   */
  virtual Double_t GetXFeynman() const;

  /**
   Sets the status code.
   */
  virtual void SetStatus(UShort_t);

  /**
   Sets the parent index, in the range [1, N - 1] for particles with
   parents, or 0 for those without.
   */
  virtual void SetParentIndex(UShort_t);

  /**
   Sets the Feynman-x
   */
  virtual void SetXFeynman(double xf);

  /**
   Sets the four-momentum of the particle.
   Changes are propagated to derived quantities.
   */
  virtual void Set4Vector(const TLorentzVector&);

  /**
   Sets the origin coordinates.
   */
  virtual void SetVertex(const TVector3&);

 protected:
  UShort_t KS;  ///< Status code: see PYTHIA manual
  UShort_t orig;  ///< I of parent particle
  Int_t id;  ///< PDG code identifying the particle
  Double32_t px;  ///< x component of momentum (GeV/c)
  Double32_t py;  ///< y component of momentum (GeV/c)
  Double32_t pz;  ///< z component of momentum (GeV/c)
  Double32_t E;  ///< Total energy (GeV)
  Double32_t p;  ///< Magnitude of momentum (GeV/c)
  Double32_t m;  ///< Invariant mass (GeV/c2)
  Double32_t pt;  ///< Momentum transverse to the beam direction (GeV/c)
  Double32_t theta;  ///< Polar angle (radians [0, pi])
  Double32_t phi;  ///< Angle of azimuth (radians [0, 2pi])
  Double32_t rapidity;  ///< Rapidity
  Double32_t eta;  ///< Pseudorapidity
  Double32_t xFeynman;  ///< Feynman x = 2 * pz / centre of mass energy
  Double32_t xv;  ///< x vertex position (cm)
  Double32_t yv;  ///< y vertex position (cm)
  Double32_t zv;  ///< z vertex position (cm)

  ClassDef(erhic::hadronic::ParticleMC, 1)
};

inline Pid ParticleMC::Id() const {
  return Pid(id);
}

inline TLorentzVector ParticleMC::Get4Vector() const {
  return TLorentzVector(GetPx(), GetPy(), GetPz(), GetE());
}

inline Double_t ParticleMC::GetPx() const {
  return px;
}

inline Double_t ParticleMC::GetPy() const {
  return py;
}

inline Double_t ParticleMC::GetPz() const {
  return pz;
}

inline Double_t ParticleMC::GetE() const {
  return E;
}

inline Double_t ParticleMC::GetP() const {
  return p;
}

inline Double_t ParticleMC::GetM() const {
  return m;
}

inline Double_t ParticleMC::GetPt() const {
  return pt;
}

inline Double_t ParticleMC::GetTheta() const {
  return theta;
}

inline Double_t ParticleMC::GetPhi() const {
  return phi;
}

inline Double_t ParticleMC::GetRapidity() const {
  return rapidity;
}

inline Double_t ParticleMC::GetEta() const {
  return eta;
}

inline TVector3 ParticleMC::GetVertex() const {
  return TVector3(xv, yv, zv);
}

inline UShort_t ParticleMC::GetStatus() const {
  return KS;
}

inline UShort_t ParticleMC::GetParentIndex() const {
  return orig;
}

inline Double_t ParticleMC::GetXFeynman() const {
  return xFeynman;
}

inline void ParticleMC::SetStatus(UShort_t i) {
  KS = i;
}

inline void ParticleMC::SetXFeynman(double xf) {
  xFeynman = xf;
}

inline void ParticleMC::SetVertex(const TVector3& v) {
  xv = v.x();
  yv = v.y();
  zv = v.z();
}

}  // namespace hadronic
}  // namespace erhic

#endif  // INCLUDE_EICSMEAR_HADRONIC_PARTICLEMC_H_
