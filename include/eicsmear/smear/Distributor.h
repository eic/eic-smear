/**
 \file
 Declaration of class Smear::Distributor.
 
 \author    Michael Savastio
 \date      2011-08-19
 \copyright 2011 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_SMEAR_DISTRIBUTOR_H_
#define INCLUDE_EICSMEAR_SMEAR_DISTRIBUTOR_H_

#include <Rtypes.h>

class TF1;
class TString;

namespace Smear {

/**
 Distribution function for random sampling.
 Used by devices to generate smearing on some distribution.
 By default, smearing is generated on a Gaussian
 of which the Monte Carlo value is the mean, and the standard deviation
 is given by the device parametrization.
 It is intended that you access this object via the methods provided
 in the Device class.
 \todo This could do with some reworking, to make it compatible with
 an arbitrary function via a functor (which I think is the best way
 to go). Move to its own files.
 */
class Distributor {
 public:
  /**
   Default constructor.
   Gaussian distribution.
   */
  Distributor();

  /**
   Initialise with a ROOT::TF1-style string.
   The function should depend on two parameters, the first of which
   specifies the function midpoint and the second of which species a width.
   See ROOT::TF1 documentation for details of the string format.
   For example, for a step function, one could use
   Distributor("abs(x-[0])<[1]");
   where [0] stands for the first parameter (midpoint) and [1] is the
   second parameter (width).
   The second and third arguments limit the range of thrown values around
   the midpoint to [midpoint - lower, midpoint + upper].
   Values of <= 0 give no limit.
   The 4th and 5th arguments define the valid range of the function.
   */
  Distributor(const TString& formula, double lower, double upper,
              double minimum = -1.e6, double maximum = 1.e6);

  /**
   Destructor.
   */
  virtual ~Distributor();

  /**
   Generate a random value based on a given midpoint and width.
   */
  virtual double Generate(double midpoint, double width);

 protected:
  double mPlus;
  double mMinus;
  TF1* mDistribution;

  ClassDef(Smear::Distributor, 1)
};

}  // namespace Smear

#endif  // INCLUDE_EICSMEAR_SMEAR_DISTRIBUTOR_H_
