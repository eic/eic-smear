// For Doxygen documentation:
/*! \file */ 

#ifndef _ERHIC_BUILDTREE_FUNCTIONS_
#define _ERHIC_BUILDTREE_FUNCTIONS_

#include <string>

#include <Rtypes.h>
#include <TString.h>

#include "Smear.h"
#include "Detector.h"

class TLorentzVector;

namespace erhic {
   class EventMC;
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
 Function for generating a ROOT TTree file from a plain-text Monte Carl file.
 */
Long64_t BuildTree(const TString& inputFileName,
                   const TString& outputDirName = ".",
                   const Long64_t maxEvent = 0,
                   std::string logFileName = "");


/**
 \fn
 Processes a ROOT Monte Carlo event file to produce a file with 
 information smeared for detector effects.
 */
int SmearTree(Smear::Detector det,
              TString inFileName,
              TString outFileName = "",
              Long64_t nEvents = -1);

/**
 Produces a DOT file describing the particle content of the event.
 */
class EventToDot {
   
public:
   
   virtual ~EventToDot() { }
   
   void Generate(const erhic::EventMC&,
                 const std::string& outputName) const;
   
   ClassDef(EventToDot, 0)
};

#endif
