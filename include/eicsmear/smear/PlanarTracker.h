/**
 \file
 Declaration of class Smear::PlanarTracker.
 
 \author    Will Foreman
 \date      2011-08-19
 \copyright 2011 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_SMEAR_PLANARTRACKER_H_
#define INCLUDE_EICSMEAR_SMEAR_PLANARTRACKER_H_

#include <Rtypes.h>  // For ClassDef

#include "eicsmear/smear/Smear.h"  // KinType
#include "eicsmear/smear/Tracker.h"

namespace erhic {

class VirtualParticle;

}  // namespace erhic

namespace Smear {

class ParticleMCS;

/**
 A cylindrical tracking detector.
 Implements both intrinsic and multiple-scattering resolution.
 */
class PlanarTracker : public Tracker {
 public:
  /**
   Default constructor.
   B = 2, NRL = 0.03, sigma(r-phi) = 0.001, N = 25, inner radius = 0.25 m,
   outer radius = 1 m, length = 6 m.
   */
  PlanarTracker();

  /**
   Constructor for a tracker with arbitrary positioning along z.
   */
  PlanarTracker(double innerRadius, double outerRadius,
                double zMin, double zMax,
                double magneticField, double nRadiationLengths,
                double sigmaRPhi, double nPlanes);

  /**
   Destructor.
   */
  virtual ~PlanarTracker();

  /**
   Returns a new copy of this Tracker.
   The argument is not used.
   */
  virtual PlanarTracker* Clone(const char* = "") const;

  /**
   Print information about this device to standard output.
   */
  virtual void Print(Option_t* = "") const;

  /**
   Returns the path length of the particle through the tracker in metres.
   */
  double L(const erhic::VirtualParticle&) const;

  /**
   Returns the transverse path length of the particle through the
   tracker in metres.
   */
  double LPrime(const erhic::VirtualParticle&) const;

  /**
   Returns the number of measurement points for the particle.
   This is the number of planes hit in the planar tracker for the particle.
   */
  virtual int NPoints(const erhic::VirtualParticle&) const;

  /**
   Returns true if the particle falls within the angular acceptance
   defined via the input parameters.
   This is defined as the particle's transverse path length being
   greater than half of (outer radius - inner radius) / number of fit points
   i.e. the particle has to pass through enough of the detector to cause
   at least one measureable point.
   */
  virtual bool Accepts(const erhic::VirtualParticle&) const;

  /**
   Returns the minimum theta of particles accepted by the tracker (radians).
   */
  virtual double GetThetaMin() const;

  /**
   Returns the maximum theta of particles accepted by the tracker (radians).
   */
  virtual double GetThetaMax() const;

 protected:
  /** 
   Compute the intersection point of the particle with a radial surface. 
   Returns the 3-vector of the intersection point if the z of the 
   intersection is within (zmin, zmax) of this detector, 
   or (0, 0, NaN) if not. 
   */ 
  TVector3 ComputeIntersectionWithRadius(const erhic::VirtualParticle&,
                                         double radius) const;
  TVector3 ComputeIntersectionWithPlane(const erhic::VirtualParticle&,
                                        double z) const;
  TVector3 ComputePath(const erhic::VirtualParticle&) const;

  double mNPlanes;  ///< Number of planes
  double mInnerRadius;  ///< Inner radius (m)
  double mOuterRadius;  ///< Outer radius (m)
  double mZMin;  ///< Lower (most negative) z face
  double mZMax;  ///< Upper (most positive) z face

  ClassDef(Smear::PlanarTracker, 1)
};

inline PlanarTracker* PlanarTracker::Clone(const char*) const {
  return new PlanarTracker(*this);
}

}  // namespace Smear

#endif  // INCLUDE_EICSMEAR_SMEAR_PLANARTRACKER_H_
