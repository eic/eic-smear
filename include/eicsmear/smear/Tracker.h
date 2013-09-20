/**
 \file
 Declaration of class Smear::Tracker.
 
 \author    Will Foreman
 \date      2011-08-19
 \copyright 2011 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_SMEAR_TRACKER_H_
#define INCLUDE_EICSMEAR_SMEAR_TRACKER_H_

#include <Rtypes.h>  // For ClassDef

#include "eicsmear/smear/Distributor.h"
#include "eicsmear/smear/Smear.h"  // KinType
#include "eicsmear/smear/Smearer.h"

namespace erhic {

class VirtualParticle;

}  // namespace erhic

namespace Smear {

class ParticleMCS;

/**
 A cylindrical tracking detector.
 Implements both intrinsic and multiple-scattering resolution.
 Abstract base class: inheriting classes must implement L(), LPrime(),
 NPoints() and Accepts().
 */
class Tracker : public Smearer {
 public:
  /**
   Constructor.
   */
  Tracker(double magneticField = 2., double nRadiationLengths = 0.01,
          double resolution = 0.001);

  /**
   Destructor.
   */
  virtual ~Tracker();

  /**
   Returns the resolution at the kinematics of this particle.
   */
  virtual double Resolution(const erhic::VirtualParticle&) const;

  /**
   Smear the properties of the input particle and store the
   smeared values in the ParticleMCS.
   */
  void Smear(const erhic::VirtualParticle&, ParticleMCS&);

  /**
   Returns the path length of the particle through the tracker in metres.
   */
  virtual double L(const erhic::VirtualParticle&) const = 0;

  /**
   Returns the transverse path length of the particle through the
   tracker in metres.
   */
  virtual double LPrime(const erhic::VirtualParticle&) const = 0;

  /**
   Returns the number of measurement points for the particle.
   */
  virtual int NPoints(const erhic::VirtualParticle&) const = 0;

  /**
   Returns true if the particle falls within the angular acceptance
   of the detector.
   */
  virtual bool Accepts(const erhic::VirtualParticle&) const = 0;

  /**
   Returns the minimum theta of particles accepted by the tracker (radians).
   */
  virtual double GetThetaMin() const = 0;

  /**
   Returns the maximum theta of particles accepted by the tracker (radians).
   */
  virtual double GetThetaMax() const = 0;

  /**
   Set whether a vertex constraint should be used when calculating the
   intrinsic resolution. Without it a factor of sqrt(720) is included in
   the resolution; with it the factor is sqrt(320). By defuault no
   constraint is assumed.
   */
  void SetVertexConstraint(bool constrain);

 protected:
  /**
   Multiple scattering contribution, given by
   delta(p)/p = 0.0136 * z * sqrt(NRL) / (0.3 * B * L * beta)
   z = charge, NRL = # radiation lengths, B = mag field, L = track length,
   beta = particle velocity.
   */
  virtual double MultipleScatteringContribution(
    const erhic::VirtualParticle&) const;

  /**
   The intrinsic resolution of the detector, depending on momentum,
   magnetic field, the detector dimensions, the number of fit points
   and the point resolution.
   */
  virtual double IntrinsicContribution(const erhic::VirtualParticle&) const;

  Int_t mFactor;  ///< Factor in intrinsic resolution calculation
                  ///< dependent on vertex constraint.
  double mMagField;  ///< Magnetic field strength in Tesla
  double mNRadLengths;  ///< Number of radiation lengths (dimensionless)
  double mSigmaRPhi;  ///< Point resolution
  Distributor Distribution;  ///< Random distribution

  ClassDef(Smear::Tracker, 1)
};

}  // namespace Smear

#endif  // INCLUDE_EICSMEAR_SMEAR_TRACKER_H_
