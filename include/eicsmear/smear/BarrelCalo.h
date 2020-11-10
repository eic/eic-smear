/**
 \file
 Declaration of class Smear::BarrelCalo.
 
 \author    Kolja Kauder
 \date      2020-11-05
 \copyright 2020 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_SMEAR_BARRELCALO_H_
#define INCLUDE_EICSMEAR_SMEAR_BARRELCALO_H_

#include <Rtypes.h>  // For ClassDef

#include "eicsmear/smear/Distributor.h"
#include "eicsmear/smear/Smear.h"  // KinType
#include "eicsmear/smear/Smearer.h"
#include "eicsmear/smear/Device.h"

namespace erhic {

class VirtualParticle;

}  // namespace erhic

namespace Smear {

class ParticleMCS;

/**
 A cylindrical calorimeter.
 Main motivation for this class is realistic angular resolution,
 but for compactness we derive from Device and allow for an E resolution string.
 Information needed is spatial resolution, such as in slide 2 here:
 https://indico.bnl.gov/event/8231/contributions/37910/attachments/28335/43607/talk_eic_yr-cal_2020_05.pdf
 plus geometric and material properties 
 */
class BarrelCalo : public Device {
 public:
  /**
   Constructor. KinType is always kE
   */
  BarrelCalo(const TString SpatialFormula, const EGenre genre,
	     const double InnerRadius, const double OuterRadius,	   
	     const bool Projective=true, const double ProjPhiAngle=0, const double ProjThetaAngle=0,
	     const TString EResFormula="");

  // delete default device ctor
  BarrelCalo(KinType, const TString& formula, EGenre) = delete;
  
  /**
   Destructor.
   */
  virtual ~BarrelCalo();

  /**
   Returns a dynamically allocated copy of this object.
   The argument is unused and is present for compatibility with
   ROOT::TObject::Clone().
   */
  virtual BarrelCalo* Clone(const char* = "") const;

  /** Copy ctor.
      Important to have, as it's used by Clone(), and Detector uses clones, not originals
   */
  BarrelCalo(const BarrelCalo& that);

  
  /**
   Smear the properties of the input particle and store the
   smeared values in the ParticleMCS.
   */
  void Smear(const erhic::VirtualParticle&, ParticleMCS&);


  // =====================
  // Math:
  // =====================
  // R = z * tan (theta), where R = distance from z-axis (inner cylinder radius)
  // <==> theta = atan ( R / z ) == atan2 ( R ,  z ) if needed
  // sigma := sigma_x = sigma_y = sigma_z

  // BARREL: R=const
  // -- sigma(phi) = sigma / R
  // -- sigma(theta) non-projective:
  // sigma(th) = | del theta / del z | sigma
  //           = R / (R^2 + z^2)  * sigma
  //           = sin^2(theta) / R * sigma

  // -- sigma(theta) projective: tower is perpendicular to trajectory, removing one sine:
  // sigma(th) = sin(theta) / R * sigma
  // 
  // 
  // -- for approximately projective cases, replace sigma with sqrt ( sigma^2 + (X*sin( angle ))^2)
  // while this could be folded over into the constant term, the material choices matter quite a bit
  //
  // --> sigma ( phi )   = 1/R            *  sqrt( sigma^2 + (X*sin( angle ))^2)
  // --> sigma ( theta ) = sin(theta) / R *  sqrt( sigma^2 + (X*sin( angle ))^2)

 protected:

  /** RadLength can stand for lambda or X0
      lambda  - HCal nuclear interaction length in mm
      X0      - EMCal radiation length in mm
  */
  double mRadLength;

  FormulaString* mSpatialFormula; ///< in mm; typically, something like "sqrt ( pow( 3.0, 2 ) / E + pow ( 1.0,2))"
  
  double mInnerRadius; ///< in mm
  double mOuterRadius=0; ///< in mm; not needed for calculations, just added for completeness

  bool mProjective;         ///< projective calorimeter?
  double mProjPhiAngle=0;    ///< deviation from perfect phi projectivity in radians
  double mProjThetaAngle=0;  ///< deviation from perfect theta projectivity in radians
  
  // From Device, we inherit
  /* KinType mSmeared;   ///< Smeared variable */
  /* TF1* mKinematicFunction; */
  /* FormulaString* mFormula;  ///< Expression for resolution standard deviation */
  /* std::vector<Smear::KinType> mDimensions;  ///< Variables on which smearing */
  /*                                           ///< is dependent (up to 4) */
  /* Distributor mDistribution;  ///< Random distribution */
  // not all of which are needed.

 private:
  // Assignment is not supported
  BarrelCalo& operator=(const BarrelCalo&) { return *this; }

  
  
  ClassDef(Smear::BarrelCalo, 1)
};

}  // namespace Smear

#endif  // INCLUDE_EICSMEAR_SMEAR_BARRELCALO_H_
