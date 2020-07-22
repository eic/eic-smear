/**
   \file
   Declaration of various Smear namespace types and constants.
 
   \author    Kolja Kauder
   \date      2020-06-30
   \copyright 2020 Brookhaven National Lab
*/

#ifndef INCLUDE_EICSMEAR_SMEAR_SMEARCONSTANTS_H_
#define INCLUDE_EICSMEAR_SMEAR_SMEARCONSTANTS_H_

namespace Smear {

  /**
     Enumerator listing particle wise kinematic variables.
     Naming is self explanitory.
     These will be used when specifying arguments and outputs
     for device parametrizations.
  */
  enum KinType {
    kE, kP, kTheta, kPhi, kPz, kPt, kInvalidKinType
  };

  /** Classes of particles */
  enum EGenre {
    kAll = 0, kElectromagnetic = 1, kHadronic = 2
  };

  /** Particle charged **/
  enum ECharge {
    kNeutral, kCharged, kAllCharges
  };

} // namespace
#endif // INCLUDE_EICSMEAR_SMEAR_SMEARCONSTANTS_H_
