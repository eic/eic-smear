/**
 Distributor.cxx

 \file
 Implementation of class Distributor.

 \author Thomas Burton 
 \date 5/9/12
 \copyright 2012 BNL. All rights reserved.
*/

#include "eicsmear/smear/Distributor.h"

#include <TF1.h>
#include <TUUID.h>

namespace Smear {

   Distributor::Distributor()
   : mPlus(0.)
   , mMinus(0.)
   , mDistribution(new TF1(TUUID().AsString(),
                   "exp(-pow(x-[0],2)/pow([1],2))", -1.e6, 1.e6)) {
      mDistribution->SetParameters(0., 1.);
   }

   Distributor::Distributor(const TString& formula, double lower, double upper,
                            double minimum, double maximum)
   : mPlus(0.)
   , mMinus(0.)
   , mDistribution(new TF1(TUUID().AsString(), formula, minimum, maximum)) {
      mDistribution->SetParameters(0., 1.);
      if(lower > 0.) {
         mMinus = lower;
      } // if
      if(upper > 0.) {
         mPlus = upper;
      } // if
   }

   Distributor::~Distributor() {
      if(mDistribution) {
         delete mDistribution;
         mDistribution = NULL;
      } // if
   }

   double Distributor::Generate(double mean, double sigma) {
      mDistribution->SetParameters(mean, sigma);
      return mDistribution->GetRandom(mean - mMinus, mean + mPlus);
   }
} // namespace Smear
