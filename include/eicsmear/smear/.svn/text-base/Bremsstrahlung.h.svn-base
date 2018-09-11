/**
 \file
 Declaration of class Smear::Bremsstrahlung.
 
 \author    Michael Savastio
 \date      2011-08-19
 \copyright 2011 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_SMEAR_BREMSSTRAHLUNG_H_
#define INCLUDE_EICSMEAR_SMEAR_BREMSSTRAHLUNG_H_

#include <memory>

#include <TF1.h>

#include "eicsmear/smear/Device.h"
#include "eicsmear/erhic/ParticleMC.h"
#include "eicsmear/smear/ParticleMCS.h"

namespace erhic {

class VirtualParticle;

}  // namespace erhic

namespace Smear {

/**
 \brief A specialized Device class for modelling radiative losses.
 */
struct Bremsstrahlung : public Device {
  /**
   Constructor.
   Photon energies are randomly generated in the range
   [epsilon, E - epsilon] for particle energy E.
   traversed is the distance of material through with the particle passes
   and radLength is the radiation length of that material (both in cm).
   */
  Bremsstrahlung(double epsilon = 0.01,
                 double traversed = 10.,
                 double radLength = 47.1);

  /**
   Copy constructor
   */
  Bremsstrahlung(const Bremsstrahlung&);

  /**
   Returns a pointer to a duplicate of this object.
   */
  virtual Bremsstrahlung* Clone(Option_t* option = "not used") const;

  /**
   Smear the properties of a Particle and assign them to a ParticleS.
   */
  virtual void Smear(const erhic::VirtualParticle&, ParticleMCS&);

 protected:

  /**
   Returns dSigmga/dK at k = x[0].
   The arguments have this form to interface with ROOT::TF1.
   The second argument is unused.
   */
  double dSigmadK(double* x, double*);

  /**
   Compute the number of photons emitted.
   */
  int NGamma();

  void FixParticleKinematics(ParticleMCS&);

  /**
   Set the radiating particle type and configure the dSigma/dK
   function.
   */
  void SetParticle(const erhic::VirtualParticle&);

  /**
   Configure the dSigma/dK function, setting the energy range over
   which to generate photons.
   The energy range is computed from the energy of the current mParticle.
   If the resultant energy range is invalid (e.g. max < min) the range
   is not set and the function returns false.
   */
  bool SetupPDF();

  std::auto_ptr<erhic::ParticleMC> mParticle;  //!< Copy of the current particle

  double mKMin;
  double mKMax;
  double mEpsilon;
  double mTraversed;
  double mRadLength;

  TF1* mPdf;  //!< dSigma/dK function

  ClassDef(Smear::Bremsstrahlung, 1)
};

}  // namespace Smear

#endif  // INCLUDE_EICSMEAR_SMEAR_BREMSSTRAHLUNG_H_
