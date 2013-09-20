/**
 \file
 Declaration of smearing global functions.
 
 \author    Michael Savastio
 \date      2011-08-19
 \copyright 2011 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_SMEAR_FUNCTIONS_H_
#define INCLUDE_EICSMEAR_SMEAR_FUNCTIONS_H_

#include <Rtypes.h>  // For Long64_t
#include <TString.h>

namespace Smear {

  class Detector;

}  // namespace Smear

/**
 \fn
 Processes a ROOT Monte Carlo event file to produce a file with 
 information smeared for detector effects.
 */
int SmearTree(const Smear::Detector&, const TString& inFileName,
              const TString& outFileName = "", Long64_t nEvents = -1);

#endif  // INCLUDE_EICSMEAR_SMEAR_FUNCTIONS_H_
