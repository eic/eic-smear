/**
 \file
 Declaration of class erhic::VirtualParticle.
 
 \author    Thomas Burton
 \date      2011-08-19
 \copyright 2011 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_ERHIC_VIRTUALPARTICLE_H_
#define INCLUDE_EICSMEAR_ERHIC_VIRTUALPARTICLE_H_

#include <TLorentzVector.h>
#include <TVector3.h>

#include "eicsmear/erhic/Pid.h"

namespace erhic {

/**
 \brief Abstract base class for a general particle.
 */
class VirtualParticle : public TObject {
 public:
  /**
   Destructor
   */
  virtual ~VirtualParticle() { }

  /**
   Returns identity information for the Particle species.
   */
  virtual Pid Id() const = 0;

  /**
   Returns the momentum-energy four-vector (px, py, pz, E).
   */
  virtual TLorentzVector Get4Vector() const = 0;

  /**
   Returns the x component of 3-momentum.
   */
  virtual Double_t GetPx() const = 0;

  /**
   Returns the y component of 3-momentum.
   */
  virtual Double_t GetPy() const = 0;

  /**
   Returns the z component of 3-momentum.
   */
  virtual Double_t GetPz() const = 0;

  /**
   Returns total energy.
   */
  virtual Double_t GetE() const = 0;

  /**
   Returns the magnitude of 3-momentum (GeV).
   */
  virtual Double_t GetP() const = 0;

  /**
   Returns invariant mass (GeV/c<sup>2</sup>).
   */
  virtual Double_t GetM() const = 0;

  /**
   Returns momentum perpendicular to the beam direction.
   */
  virtual Double_t GetPt() const = 0;

  /**
   Returns the polar angle in the range [0, pi] radians.
   */
  virtual Double_t GetTheta() const = 0;

  /**
   Returns the polar angle in the range [0, 2pi] radians.
   */
  virtual Double_t GetPhi() const = 0;

  /**
   Returns the rapidity.
   */
  virtual Double_t GetRapidity() const = 0;

  /**
   Returns the pseudorapidity.
   */
  virtual Double_t GetEta() const = 0;

  /**
   Returns the origin point of the particle in cm.
   (0,0,0) indicates a particle originating in the collision.
   */
  virtual TVector3 GetVertex() const = 0;

  /**
   A general "status" code for the particle
   (definition depends on implementation).
   */
  virtual UShort_t GetStatus() const = 0;

  /**
   Returns the index of the parent particle
   */
  virtual UShort_t GetParentIndex() const = 0;

  /**
   Sets the origin coordinates
   */
  virtual void SetVertex(const TVector3& /* (x, y, z) */) = 0;

  /**
   Sets the four-momentum of the particle.
   Changes are propagated to derived quantities.
   */
  virtual void Set4Vector(const TLorentzVector& /* (px, py, pz, E) */) = 0;

  ClassDef(erhic::VirtualParticle, 1)
};

}  // namespace erhic

#endif  // INCLUDE_EICSMEAR_ERHIC_VIRTUALPARTICLE_H_
