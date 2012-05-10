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
    Specialized smearing class for simulating tracking detectors.
    Implements both intrinsic and multiple-scattering resolution
    based on specified detector dimensions and properties.
    \todo Remove SetParticle(), P and bMoreCrit - can pass Particle as function
    argument when required, as have a ComputePathLength(Particle) method
    serving bMoreCrit's function.
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
       Constructor.
       Initialise with all detector characteristics.
       Lengths are in metres.
       */
      Tracker(double innerRadius, double outerRadius, double length,
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
       Polar angle of the end of the outer surface.
      */
      double ThetaCrit();

      /**
       Returns the minimum theta of particle accepted by a tracker with
       the length and radii you specified.
       If you must set a different acceptance in theta, it is recomended
       that you allow it no closer to the beam line than this.
       */
      double GetThetaMin();

   protected:

      double MultipleScatteringContribution(const erhic::VirtualParticle&);
      double IntrinsicContribution(const erhic::VirtualParticle&);
      double EvaluateRes(const erhic::VirtualParticle&);
      double L(const erhic::VirtualParticle&);
      double LPrime(const erhic::VirtualParticle&);
      double tanTheta(const erhic::VirtualParticle&);

      double mMagField; ///< Magnetic field strength in Tesla
      double mNRadLengths; ///< Number of radiation lengths (dimensionless)
      double mSigmaRPhi; ///< Point resolution
      double mNFitPoints; ///< Number of fit points
      double mInnerRadius; ///< Inner radius (m)
      double mOuterRadius; ///< Outer radius (m)
      double mLength; ///< Total tracker length (m)
      Distributor Distribution; ///< Random distribution

      ClassDef(Smear::Tracker, 1)
   };

   inline double Tracker::tanTheta(const erhic::VirtualParticle& p) {
      return fabs(tan(p.GetTheta()));
   }

   inline double Tracker::ThetaCrit() {
      return atan(2. * mOuterRadius / mLength);
   }

   inline double Tracker::GetThetaMin() {
      return atan(2. * mInnerRadius / mLength);
   }

   inline Tracker* Tracker::Clone(const char*) const {
      return new Tracker(*this);
   }
} // namespace Smear

#endif // _EICSMEAR_SMEAR_TRACKER_H_
