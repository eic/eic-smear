/**
 Tracker.h
 
 \file
 Declaration of class Tracker.
 
 \author Thomas Burton
 \date 5/8/12
 \copyright 2012 BNL. All rights reserved.
 */

#ifndef _EICSMEAR_SMEAR_TRACKER_H_
#define _EICSMEAR_SMEAR_TRACKER_H_

#include <Rtypes.h> // For ClassDef

#include "eicsmear/smear/Distributor.h"
#include "eicsmear/smear/Smear.h" // KinType
#include "eicsmear/smear/Smearer.h"

namespace erhic {
   class VirtualParticle;
} // namespace erhic

namespace Smear {

   class ParticleMCS;

   /**
    A cylindrical tracking detector.
    Implements both intrinsic and multiple-scattering resolution.
    */
   class Tracker : public Smearer {
   public:

      /**
       Default constructor.
       B = 2, NRL = 0.03, sigma(r-phi) = 0.001, N = 25, inner radius = 0.25 m,
       outer radius = 1 m, length = 6 m.
      */
      Tracker();

      /**
       Constructor for a tracker symetrically positioned around z = 0.
       Initialise with all detector characteristics.
       Lengths are in metres.
       */
      Tracker(double innerRadius, double outerRadius, double length,
              double magneticField, double numberOfRadiationLengths,
              double sigmaRPhi, double numberOfPoints);

      /**
       Constructor for a tracker with arbitrary positioning along z.
      */
      Tracker(double innerRadius, double outerRadius,
              double zMin, double zMax,
              double magneticField, double numberOfRadiationLengths,
              double sigmaRPhi, double numberOfPoints);

      /** Destructor */
      virtual ~Tracker();

      /**
       Returns a new copy of this Tracker.
       The argument is not used.
       */
      virtual Tracker* Clone(const char* = "") const;

      /**
       Smear the kinematics of prt and store the result in prtOut.
       */
      void Smear(const erhic::VirtualParticle&, ParticleMCS&);

      /**
       Returns the minimum theta of particles accepted by the tracker (radians).
       */
      double GetThetaMin();

      /**
       Returns the maximum theta of particles accepted by the tracker (radians).
      */
      double GetThetaMax();

      /**
       Print information about this device to standard output.
      */
      virtual void Print(Option_t* = "") const;

      /**
       Returns the resolution at the kinematics of this particle.
      */
      double Resolution(const erhic::VirtualParticle&) const;

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
       Returns true if the particle falls within the angular acceptance
       defined via the input parameters.
       This is defined as the particle's transverse path length being
       greater than half of (outer radius - inner radius) / number of fit points
       i.e. the particle has to pass through enough of the detector to cause
       at least one measureable point.
      */
      bool Accepts(const erhic::VirtualParticle&) const;

   protected:

      /**
       Multiple scattering contribution, given by
       delta(p)/p = 0.0136 * z * sqrt(NRL) / (0.3 * B * L * beta)
       z = charge, NRL = # radiation lengths, B = mag field, L = track length,
       beta = particle velocity.
      */
      double MultipleScatteringContribution(const erhic::VirtualParticle&) const;

      /**
       The intrinsic resolution of the detector, depending on momentum,
       magnetic field, the detector dimensions, the number of fit points
       and the point resolution.
      */
      double IntrinsicContribution(const erhic::VirtualParticle&) const;

      /**
       Compute the intersection point of the particle with a radial surface.
       Returns the 3-vector of the intersection point if the z of the
       intersection is within (zmin, zmax) of this detector,
       or (0, 0, NaN) if not.
      */
      TVector3 ComputeIntersectionWithRadius(
         const erhic::VirtualParticle&, double radius) const;

      /**
       Compute the intersection point of the particle with an x-y plane at z.
       Returns the 3-vector of the intersetion point if the radius of the
       intersection is within (inner radius, outer radius), or
       (0, 0, NaN) if not.
      */
      TVector3 ComputeIntersectionWithPlane(
         const erhic::VirtualParticle&, double z) const;

      /**
       Computes the path vector, defined as (v2 - v1), where v1 and v2
       are the position vectors of the particle's intersections with
       the cylinder's surface.
       Returns (0, 0, 0) for particles that don't intersect, or if
       something goes wrong.
      */
      TVector3 ComputePath(const erhic::VirtualParticle&) const;

      double mMagField; ///< Magnetic field strength in Tesla
      double mNRadLengths; ///< Number of radiation lengths (dimensionless)
      double mSigmaRPhi; ///< Point resolution
      double mNFitPoints; ///< Number of fit points
      double mInnerRadius; ///< Inner radius (m)
      double mOuterRadius; ///< Outer radius (m)
      double mZMin; ///< Lower (most negative) z face
      double mZMax; ///< Upper (most positive) z face
      Distributor Distribution; ///< Random distribution

      ClassDef(Smear::Tracker, 1)
   };

   inline Tracker* Tracker::Clone(const char*) const {
      return new Tracker(*this);
   }
} // namespace Smear

#endif // _EICSMEAR_SMEAR_TRACKER_H_
