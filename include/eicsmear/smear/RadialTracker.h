/**
 \file
 Declaration of class Smear::RadialTracker.
 
 \author    Will Foreman
 \date      2011-08-19
 \copyright 2011 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_SMEAR_RADIALTRACKER_H_
#define INCLUDE_EICSMEAR_SMEAR_RADIALTRACKER_H_

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
class RadialTracker : public Tracker {
 public:
  /**
   Default constructor.
   B = 2, NRL = 0.03, sigma(r-phi) = 0.001, N = 25, inner radius = 0.25 m,
   outer radius = 1 m, length = 6 m.
   */
  RadialTracker();

  /**
   Constructor for a tracker with arbitrary positioning along z.
   */
  RadialTracker(double innerRadius, double outerRadius,
                double zMin, double zMax,
                double magneticField, double numberOfRadiationLengths,
                double sigmaRPhi, double numberOfPoints);

  /**
   Destructor.
   */
  virtual ~RadialTracker();

  /**
   Returns a new copy of this Tracker.
   The argument is not used.
   */
  virtual RadialTracker* Clone(const char* = "") const;

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

  /**
   Compute the intersection point of the particle with an x-y plane at z.
   Returns the 3-vector of the intersetion point if the radius of the
   intersection is within (inner radius, outer radius), or
   (0, 0, NaN) if not.
   */
  TVector3 ComputeIntersectionWithPlane(const erhic::VirtualParticle&,
                                        double z) const;

  /**
   Computes the path vector, defined as (v2 - v1), where v1 and v2
   are the position vectors of the particle's intersections with
   the cylinder's surface.
   Returns (0, 0, 0) for particles that don't intersect, or if
   something goes wrong.
   */
  TVector3 ComputePath(const erhic::VirtualParticle&) const;

  double mNFitPoints;  ///< Number of fit points
  double mInnerRadius;  ///< Inner radius (m)
  double mOuterRadius;  ///< Outer radius (m)
  double mZMin;  ///< Lower (most negative) z face
  double mZMax;  ///< Upper (most positive) z face

  ClassDef(Smear::RadialTracker, 1)
};

inline RadialTracker* RadialTracker::Clone(const char*) const {
  return new RadialTracker(*this);
}

}  // namespace Smear

#endif  // INCLUDE_EICSMEAR_SMEAR_RADIALTRACKER_H_
