/**
 \file
 Declaration of class Smear::Device.
 
 \author    Michael Savastio
 \date      2011-08-19
 \copyright 2011 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_SMEAR_DEVICE_H_
#define INCLUDE_EICSMEAR_SMEAR_DEVICE_H_

#include <cmath>
#include <vector>

#include <Math/ParamFunctor.h>  // For ROOT::TMath::ParamFunctor
#include <TF1.h>
#include <TF2.h>
#include <TRandom3.h>
#include <TString.h>

#include "eicsmear/smear/Acceptance.h"
#include "eicsmear/smear/Distributor.h"
#include "eicsmear/smear/Smear.h"
#include "eicsmear/smear/Smearer.h"

class TRandom3;

namespace erhic {

class VirtualParticle;

}  // namespace erhic

namespace Smear {

class FormulaString;
class ParticleMCS;

/**
 Performs smearing of a single kinematic variable according to a simple
 expression defined via a string.
 */
class Device : public Smearer {
 public:
  /**
   Constructor.
   The first argument is the type of kinematic variable to smear.
   The second argument is a formula giving the width of the
   resolution in the variable selected with the first argument, i.e.
   sigma(A) = f(B, C...)
   where A, B, C... are selected from: E, P, theta, phi, pZ and pT, A
   is the variable type given for the first argument and B, C... are
   the variables listed in the formula.
   For example, for resolution in pT of 1% pT times sin of polar angle:
   Smear::Device(Smear::kPt, "0.01 * pT * sin(theta)");
   See ROOT::TFormula for the form of valid expressions.
   Formulae can be a function of up to four variables.
   The third argument allows selection of the types of particles that are
   smeared: electromagnetic, hadronic or all.
   */
  Device(KinType = kE, const TString& formula = "0", EGenre = kAll);

  /**
   Constructor for smearing with an arbitrary function of a single variable.
   The first argument is a function of E, P, theta, phi, pT, or pZ.
   See ROOT::TFormula for syntax.
   For example, to smear in 1/pT:
   Smear::Device('1/pT', '<some resolution function>')
   */
  Device(const TString&, const TString& resolution = "0", EGenre = kAll);

  /**
   Copy constructor.
   */
  Device(const Device&);

  /**
   Destructor.
   */
  virtual ~Device();

  /**
   Returns a dynamically allocated copy of this object.
   The argument is unused and is present for compatibility with
   ROOT::TObject::Clone().
   */
  virtual Device* Clone(const char* = "") const;

  /**
   Smear the kinematic value of the input particle and store the
   result in the ParticleMCS.
   Smearing works in the following way. If we smear a variable
   VirtualParticle.X == x with parametrization f(x), 
   then z[x,f(x)] will be stored in ParticleMCS.X, where z is a
   randomly generated number from a distribution of which x is
   the mean and f(xi) is the standard deviation.
   By default a Gaussian distribution is used.
   Use SetDistribution() for other distributions
   (and see Smear::Distributor).
   */
  virtual void Smear(const erhic::VirtualParticle&, ParticleMCS&);

  /**
   Set the random distribution from which to sample smeared kinematics.
   By default a Gaussian distribution is used.
   */
  virtual void SetDistribution(const Distributor&);

  /**
   Print information about this device to standard output.
   */
  virtual void Print(Option_t* = "") const;

 protected:
  bool Init(const TString&, const TString&, int);

  KinType mSmeared;   ///< Smeared variable
  TF1* mKinematicFunction;
  FormulaString* mFormula;  ///< Expression for resolution standard deviation
  std::vector<Smear::KinType> mDimensions;  ///< Variables on which smearing
                                            ///< is dependent (up to 4)
  Distributor mDistribution;  ///< Random distribution

 private:
  // Assignment is not supported
  Device& operator=(const Device&) { return *this; }

  ClassDef(Smear::Device, 1)
};

inline void Device::SetDistribution(const Distributor& d) {
  mDistribution = d;
}

}  // namespace Smear

#endif  // INCLUDE_EICSMEAR_SMEAR_DEVICE_H_
