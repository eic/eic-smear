/**
 \file
 Declarations of kinematic calculation classes.
 
 \author    Thomas Burton
 \date      2011-07-07
 \copyright 2011 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_ERHIC_KINEMATICS_H_
#define INCLUDE_EICSMEAR_ERHIC_KINEMATICS_H_

#include <list>
#include <vector>

#include <Rtypes.h>
#include <TLorentzVector.h>

#include "eicsmear/erhic/ParticleMC.h"
#include "eicsmear/erhic/ParticleIdentifier.h"
#include "eicsmear/erhic/VirtualParticle.h"

namespace erhic {

class EventDis;
class VirtualParticle;

/**
 A collection of DIS kinematic variables.
 */
struct DisKinematics : public TObject {
  /**
   Default constructor.
   */
  DisKinematics();
  /**
   Constructor initialising each variable.
   */
  DisKinematics(double x, double y, double nu, double Q2, double W2);
  Double32_t mX;
  Double32_t mQ2;
  Double32_t mW2;
  Double32_t mNu;
  Double32_t mY;

  ClassDef(erhic::DisKinematics, 1)
};

/**
 Abstract base class for computations of event kinematics.
 */
class KinematicsComputer {
 public:
  virtual ~KinematicsComputer() { }
  virtual TObject* Calculate() = 0;
  ClassDef(erhic::KinematicsComputer, 1)
};

/**
 Computes DIS event kinematics from the scattered lepton.
 Uses lepton momentum in place of energy if momentum is available,
 as this is typically measured more precisely for EIC kinematics.
 */
class LeptonKinematicsComputer : public KinematicsComputer {
 public:
  virtual ~LeptonKinematicsComputer() { }
  /**
   Determine the beam info from the input event
   */
  explicit LeptonKinematicsComputer(const EventDis&);
  virtual DisKinematics* Calculate();

 protected:
  std::vector<const VirtualParticle*> mBeams;

  ClassDef(erhic::LeptonKinematicsComputer, 1)
};

/**
 Computes DIS event kinematics from final-state hadrons using
 the Jacquet-Blondel method.
 \todo Revisit implementation, giving option for using particle energy
 or momentum when computing "energy", and think how to handle mass for
 smeared particles.
 */
class JacquetBlondelComputer : public KinematicsComputer {
 public:
  virtual ~JacquetBlondelComputer();
  /**
   Initialise with the event to compute.
   If the second argument is non-NULL, use the beam information from it
   in the computation.
   If it is NULL, determine the beam information automatically from
   the event.
   This allows the same class to be used with smeared calculations, where
   the beam information isn't associated with the smeared event itself.
   */
  explicit JacquetBlondelComputer(const EventDis&);
  virtual DisKinematics* Calculate();

 protected:
  virtual Double_t ComputeY() const;
  virtual Double_t ComputeQSquared() const;
  virtual Double_t ComputeX() const;
  /// The event for which kinematics are being calculated.
  const EventDis& mEvent;
  /// Array of final-state particles used in computing kinematics.
  std::vector<const VirtualParticle*> mParticles;

  ClassDef(erhic::JacquetBlondelComputer, 1)
};

/**
 Computes DIS event kinematics from a mixture of hadronic and lepton
 variables using the double-angle method.
 */
class DoubleAngleComputer : public KinematicsComputer {
 public:
  virtual ~DoubleAngleComputer();
  /**
   Initialise with the event to compute.
   If the second argument is non-NULL, use the beam information from it
   in the computation.
   If it is NULL, determine the beam information automatically from
   the event.
   This allows the same class to be used with smeared calculations, where
   the beam information isn't associated with the smeared event itself.
   */
  explicit DoubleAngleComputer(const EventDis&);
  virtual DisKinematics* Calculate();

 protected:
  const EventDis& mEvent;
  /**
   Scattering angle of struck quark.
   */
  virtual Double_t ComputeQuarkAngle() const;
  virtual Double_t ComputeY() const;
  virtual Double_t ComputeQSquared() const;
  virtual Double_t ComputeX() const;
  /// Stores whether the particle list has changed since the last
  /// computation of the quark angle.
  mutable Bool_t mHasChanged;
  /// Caches the quark angle
  mutable Double_t mAngle;
  std::vector<const VirtualParticle*> mParticles;

  ClassDef(erhic::DoubleAngleComputer, 1)
};

}  // namespace erhic

#endif  // INCLUDE_EICSMEAR_ERHIC_KINEMATICS_H_
