/**
 \file
 Global function declarations.
 
 \author    Thomas Burton
 \date      2011-07-07
 \copyright 2011 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_FUNCTIONS_H_
#define INCLUDE_EICSMEAR_FUNCTIONS_H_

#include <string>

#include <Rtypes.h>
#include <TString.h>

#include "eicsmear/smear/Smear.h"
#include "eicsmear/smear/Detector.h"

class TLorentzVector;

namespace erhic {
  class EventMC;
  class VirtualEventFactory;
  /** 
      Simple namespace-wide constant to determine the version
  */
  const std::string EicSmearVersionString = "1.1.9";

}

/**
 \fn
 Returns the first non-blank character in a line.
 Returns a NULL termination character if there are no non-blank characters
 in the line.
 */
char getFirstNonBlank(const std::string&);

/**
 \fn
 Calculate the hadron azimuthal angle around the virtual photon direction with
 respect to the lepton scattering plane in the proton rest frame.
 We use the HERMES convention, returning an angle in the range [0,2pi].
 The vectors passed as arguments should already be boosted to the proton rest
 frame. Incident and scattered leptons, incident protons and virtual photons
 all return -999.
 */
double computeHermesPhiH(const TLorentzVector& hadronInPrf,
                         const TLorentzVector& leptonInPrf,
                         const TLorentzVector& photonInPrf);

// Forward declarations of functions for which we wish to build dictionaries
// using rootcint

/**
 \fn
 Function for generating a ROOT TTree file from a plain-text Monte Carlo file.
 */
Long64_t BuildTree(const std::string& inputFileName,
                   const std::string& outputDirName = ".",
                   const Long64_t maxEvent = 0,
                   const std::string& logFileName = "");

/**
 \fn
 Function for generating a HepMC file from en EICTree ROOT file.
 */
Long64_t TreeToHepMC(const std::string& inputFileName,
		     const std::string& outputDirName = ".",
		     Long64_t maxEvent = 0,
		     const bool createHepMC2=false);

/**
 \fn
 Identical to TreeToHepMC with an added tweak to filter on certain events.
 Currently, only accept events with a charmed final state
 */
Long64_t FilterTreeToHepMC(const std::string& inputFileName,
			   const std::string& outputDirName = ".",
			   Long64_t maxEvent = 0);

/**
 Produces a DOT file describing the particle content of the event.
 */
class EventToDot {
 public:
  /**
   Destructor.
   */
  virtual ~EventToDot() { }

  /**
   Write a DOT file describing an event.
   */
  void Generate(const erhic::EventMC&, const std::string& outputName) const;

  ClassDef(EventToDot, 0)
};

#endif  // INCLUDE_EICSMEAR_FUNCTIONS_H_
