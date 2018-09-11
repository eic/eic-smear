/**
 \file
 Implementation of class Smear::Distributor.
 
 \author    Michael Savastio
 \date      2011-08-19
 \copyright 2011 Brookhaven National Lab
 */

#include "eicsmear/smear/Distributor.h"

#include <TF1.h>
#include <TRandom.h>
#include <TUUID.h>

namespace Smear {

Distributor::Distributor()
: mPlus(0.)
, mMinus(0.)
, mDistribution(NULL) {
}

Distributor::Distributor(const TString& formula, double lower, double upper,
                         double minimum, double maximum)
: mPlus(0.)
, mMinus(0.)
, mDistribution(new TF1(TUUID().AsString(), formula, minimum, maximum)) {
  mDistribution->SetParameters(0., 1.);
  if (lower > 0.) {
    mMinus = lower;
  }  // if
  if (upper > 0.) {
    mPlus = upper;
  }  // if
}

Distributor::~Distributor() {
  if (mDistribution) {
    delete mDistribution;
    mDistribution = NULL;
  }  // if
}

double Distributor::Generate(double mean, double sigma) {
  double random(0.);
  if (!mDistribution) {
    random = gRandom->Gaus(mean, sigma);
  } else {
    mDistribution->SetParameters(mean, sigma);
    if (mMinus > 0. || mPlus > 0.) {
      random = mDistribution->GetRandom(mean - mMinus, mean + mPlus);
    } else {
      random = mDistribution->GetRandom();
    }  // if
  }  // if
  return random;
}

}  // namespace Smear
