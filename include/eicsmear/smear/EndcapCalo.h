/**
 \file
 Declaration of class Smear::EndcapCalo.
 
 \author    Kolja Kauder
 \date      2020-11-05
 \copyright 2020 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_SMEAR_ENDCAPCALO_H_
#define INCLUDE_EICSMEAR_SMEAR_ENDCAPCALO_H_

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
class EndcapCalo : public Device {
 public:
  /**
   Constructor. KinType is always kE
   */
  EndcapCalo(const TString SpatialFormula, const EGenre genre,
	     const double Zposition,
	     const bool Projective=false, const double ProjPhiAngle=0, const double ProjThetaAngle=0,
	     const TString EResFormula="");

  // delete standard device ctor
  EndcapCalo(KinType, const TString& formula, EGenre) = delete;

  /** 
      Default constructor
  */
 EndcapCalo() :
  EndcapCalo ( "", kAll, 2000 ){};
  
  
  /**
   Destructor.
   */
  virtual ~EndcapCalo();

  /**
   Returns a dynamically allocated copy of this object.
   The argument is unused and is present for compatibility with
   ROOT::TObject::Clone().
   */
  virtual EndcapCalo* Clone(const char* = "") const;

  /** Copy ctor.
      Important to have, as it's used by Clone(), and Detector uses clones, not originals
   */
  EndcapCalo(const EndcapCalo& that);

  
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
  // ENDCAP: z=const
  // theta = atan ( sqrt ( x^2 + y^2 ) /z )
  // phi = atan ( y/x )
  // with sigma_x = sigma_y =: sigma

  // -- sigma( phi ) = sqrt ( (del phi / del x )^2 + (del phi / del y )^2 ) * sigma = 1/R * sigma 
  //              = sigma / (z * tan (theta ) ) 
    
  // -- sigma(th) =  sqrt ( (del theta / del x )^2 + (del theta / del y )^2 ) * sigma
  //              = z / (R^2 + z^2)  * sigma
  //              = cos^2(theta) / z * sigma
  // https://www.wolframalpha.com/input/?i=diff+%28+atan+%28+R%2Fz%29%2C+R%29
  // https://www.wolframalpha.com/input/?i=simplify+%28+z+%2F+%28R%5E2%2Bz%5E2%29%2C+R+%3D+z+*+tan%28theta%29%29
  // 
  // However, there is uncertainty in the z direction, for non-projective endcap.
  // again, account for by replacing sigma with sqrt( sigma^2 + (X*sin( angle ))^2)
  // --> sigma(th) = cos^2(theta)/z * sqrt( sigma^2 + (X*sin( angle ))^2)
  // where now angle is theta



 protected:

  /** RadLength can stand for lambda or X0
      lambda  - HCal nuclear interaction length in mm
      X0      - EMCal radiation length in mm
  */
  double mRadLength;

  FormulaString* mSpatialFormula; ///< in mm; typically, something like "sqrt ( pow( 3.0, 2 ) / E + pow ( 1.0,2))"
  
  double mZposition; ///< in mm; can be positive or negative
  
  // We're currently not foreseeing projective endcaps
  // and the functionality will not be implemented.
  // Just using this as a placeholder for future safety
  bool mProjective=false;         ///< projective calorimeter?
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
  EndcapCalo& operator=(const EndcapCalo&) { return *this; }

  
  
  ClassDef(Smear::EndcapCalo, 1)
};

}  // namespace Smear

#endif  // INCLUDE_EICSMEAR_SMEAR_ENDCAPCALO_H_
