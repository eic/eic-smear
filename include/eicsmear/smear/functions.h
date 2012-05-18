/**
 functions.h

 \file
 Forward declarations of functions in namespace Smear that require CINT
 dictionaries.

 \author Thomas Burton 
 \date 5/8/12
 \copyright 2012 BNL. All rights reserved.
*/

#ifndef _EICSMEAR_SMEAR_FUNCTIONS_H_
#define _EICSMEAR_SMEAR_FUNCTIONS_H_

#include <Rtypes.h> // For Long64_t
#include <TString.h>

namespace Smear {
   class Detector;
} // namespace Smear

/**
 \fn
 Processes a ROOT Monte Carlo event file to produce a file with 
 information smeared for detector effects.
 */
int SmearTree(const Smear::Detector&, const TString& inFileName,
              const TString& outFileName = "", Long64_t nEvents = -1);

#endif // _EICSMEAR_SMEAR_FUNCTIONS_H_
