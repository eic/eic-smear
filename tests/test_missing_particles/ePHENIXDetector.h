/**
 \file ePHENIXDetector.cpp
 Example smearing script for the ePHENIX detector

 \author    Thomas Burton
 \date      2014-03-21
 \copyright 2014 Brookhaven National Lab
 */

/*
 Example smearing script for the ePHENIX detector.
 
 It defines:
  - EPhenixMomentum: a special smearing class for the ePHENIX momentum
                     performance
  - BuildDetector: a function to generate the full ePHENIX detector description
 
 This script must be compiled in ROOT before running.
 Therefore, you must first make sure that
 (1) libeicsmear is loaded, via
       gSystem->Load("libeicsmear");
 (2) the eicsmear headers are accessible, by doing e.g.
       gSystem->AddIncludePath(" -I/path/to/eic-smear/include");
 You can then compile it in a ROOT session via
       .L ePHENIXDetector.cpp+g
 */

#include <algorithm>  // For std::max
#include <cmath>
#include <list>

#include <TCollection.h>  // For TIter
#include <TGraph.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <TRandom.h>

#include "eicsmear/erhic/VirtualParticle.h"
#include "eicsmear/smear/Acceptance.h"
#include "eicsmear/smear/Device.h"
#include "eicsmear/smear/Detector.h"
#include "eicsmear/smear/Smearer.h"
#include "eicsmear/smear/ParticleMCS.h"
#include "eicsmear/smear/PerfectID.h"

/**
 Smearing class describing ePHENIX momentum resolution.
 
 The ePHENIX momentum resolution is too complicated to handle with a simple
 parameterisation via the Smear::Device class.
 Therefore we define a custom Smearer class to implement the resolution.
 See the email in comments at the end of the file for details of the
 resolution values we use.
 */
class EPhenixMomentum : public Smear::Smearer {
public:
    /**
     Destructor.
     */
    ~EPhenixMomentum() override;

    /**
     Constructor.

     If multipleScattering is true, apply the multiple scattering resolution in
     the region where it is known, 2 < eta < 4.
     Otherwise apply only the linear resolution term.
     */
    explicit EPhenixMomentum(bool multipleScattering = true);

    /**
     Initialise the object

     This is called automatically by the smearing routine, so the user needn't
     call it themselves before smearing.
     */
    void Initialise();

    /**
     Returns a pointer to the graphs.

     Initialises all graphs if not yet done.
     */
    TMultiGraph *Graphs();

    /**
     Duplicate this object.
     */
    EPhenixMomentum *Clone(const char * /* unused */) const override;

    /**
     Smears the input ParticleMC and stores the results in the ParticleMCS.
     */
    void Smear(const erhic::VirtualParticle &particle, Smear::ParticleMCS &smeared) override;

    /**
     Return sigma(P) due to multiple scattering

     Note this means sigma(P) (GeV) *not* sigma(P)/P
     */
    virtual double computeMultipleScattering(const erhic::VirtualParticle &particle) const;

    /**
     Draw all graphs
     */
    virtual void Draw(Option_t *option = "ac") override;

private:
    TMultiGraph *mGraphDrawer;  ///< Collection of graphs for each eta range
    Bool_t mMultipleScattering;  ///< Flag to include multiple scattering

    void Streamer(TBuffer &R__b) override;
};

/**
 Smearing parameterisations for the ePHENIX detector.
 
 These parameterisations are non-exhaustive: they do not cover elements such
 as particle identification, and they are only for the central elements of
 the detector - essentially, just TPC and B/EEMC.
 
 If multipleScattering == true, apply multiple scattering term to momentum
 resolution (currently only implemented for 2 < eta < 4). Otherwise just use
 the linear resolution term.

 Note: you must gSystem->Load("libeicsmear") BEFORE loading this script,
 as ROOT needs to understand what a Smear::Detector is.
 */
Smear::Detector BuildEphoenix(bool multipleScattering = true);
