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
  virtual ~EPhenixMomentum();
  /**
     Constructor.
   
     If multipleScattering is true, apply the multiple scattering resolution in
     the region where it is known, 2 < eta < 4.
     Otherwise apply only the linear resolution term.
  */
  EPhenixMomentum(bool multipleScattering = true);
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
  TMultiGraph* Graphs();
  /**
     Duplicate this object.
  */
  virtual EPhenixMomentum* Clone(const char* /* unused */) const;
  /**
     Smears the input ParticleMC and stores the results in the ParticleMCS.
  */
  virtual void Smear(const erhic::VirtualParticle& particle,
                     Smear::ParticleMCS& smeared);
  /**
     Return sigma(P) due to multiple scattering
   
     Note this means sigma(P) (GeV) *not* sigma(P)/P
  */
  virtual double computeMultipleScattering(
					   const erhic::VirtualParticle& particle) const;
  /**
     Draw all graphs
  */
  virtual void Draw(Option_t* option = "ac");

  /**
     Helper function to convert eta to theta (radians)
     
     Detector acceptances require theta, not eta
  */
  static double etaToTheta(double eta) {
    return 2. * std::atan(std::exp(-eta));
  }


private:
  TMultiGraph* mGraphDrawer;  ///< Collection of graphs for each eta range
  Bool_t mMultipleScattering;  ///< Flag to include multiple scattering
  // ClassDef(EPhenixMomentum, 1)
};

EPhenixMomentum::EPhenixMomentum(bool multipleScattering)
  : mGraphDrawer(nullptr), mMultipleScattering(multipleScattering) {
  // Only use charged particles, as this is a tracker
  Accept.SetCharge(Smear::kCharged);
}

EPhenixMomentum::~EPhenixMomentum() {
  if (mGraphDrawer) {
    delete mGraphDrawer;
    mGraphDrawer = NULL;
  }  // if
}

void EPhenixMomentum::Initialise() {
  // Pseudorapidity values
  const double eta[75] = {
    -3.0, -2.9, -2.8, -2.7, -2.6, -2.5, -2.4, -2.3, -2.2, -2.1, -2.0, -2.0,
    -1.9, -1.8, -1.7, -1.6, -1.5, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9,
    -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,  0.0,  0.1,  0.2,  0.3,
    0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1.0,  1.0,  1.1,  1.2,  1.3,  1.4,
    1.5,  1.6,  1.7,  1.8,  1.9,  2.0,  2.1,  2.2,  2.3,  2.4,  2.5,  2.5,
    2.6,  2.7,  2.8,  2.9,  3.0,  3.1,  3.2,  3.3,  3.4,  3.5,  3.6,  3.7,
    3.8,  3.9,  4.0
  };
  // sigma(1/p) values at each value of eta
  const double sigma[75] = {
    0.010500, 0.009820, 0.009140, 0.008460, 0.007780, 0.007100, 0.006520,
    0.005940, 0.005360, 0.004780, 0.004200, 0.009811, 0.008759, 0.007479,
    0.006242, 0.005436, 0.004524, 0.006887, 0.005374, 0.004485, 0.003822,
    0.003164, 0.002592, 0.002791, 0.002991, 0.003187, 0.003374, 0.003547,
    0.003700, 0.003827, 0.003921, 0.003980, 0.004000, 0.003980, 0.003921,
    0.003827, 0.003700, 0.003547, 0.003374, 0.003187, 0.002991, 0.002791,
    0.002592, 0.001200, 0.001280, 0.001360, 0.001440, 0.001520, 0.001600,
    0.001840, 0.002080, 0.002320, 0.002560, 0.002800, 0.003160, 0.003520,
    0.003880, 0.004240, 0.004600, 0.002300, 0.002620, 0.002940, 0.003260,
    0.003580, 0.003900, 0.004400, 0.004900, 0.005400, 0.005900, 0.006400,
    0.007220, 0.008040, 0.008860, 0.009680, 0.010500
  };
  mGraphDrawer = new TMultiGraph;
  mGraphDrawer->SetTitle(";#eta;Momentum Resolution, #delta(1/p) (1/GeV)");
  // Use different graphs for each region to avoid problems at discontinuities
  // both with the calculated values and with drawing
  int points[5] = {11, 6, 26, 16, 16};  // Number of points in each graph
  int colours[5] = {kBlue, kMagenta, kRed, kGreen + 3, kGreen - 2};
  int position(0);  // Tracks cumulative array index
  for (int i(0); i < 5; ++i) {
    TGraph* graph = new TGraph(points[i], eta + position, sigma + position);
    TString name("graph");
    name += i;  // Append graph number to name
    graph->SetName(name);
    // Apply formatting
    graph->SetLineWidth(3);
    graph->SetLineColor(colours[i]);
    graph->SetMarkerColor(colours[i]);
    mGraphDrawer->Add(graph);
    position += points[i];  // Update total array index
  }  // for
}

TMultiGraph* EPhenixMomentum::Graphs() {
  if (!mGraphDrawer) {
    Initialise();
  }  // if
  return mGraphDrawer;
}

EPhenixMomentum* EPhenixMomentum::Clone(const char*) const {
  EPhenixMomentum* clone = new EPhenixMomentum(*this);
  // We need to "de-initialise" the clone so it doesn't retain the
  // same graphs as this instance.
  clone->mGraphDrawer = NULL;
  return clone;
}

void EPhenixMomentum::Smear(const erhic::VirtualParticle& particle,
                            Smear::ParticleMCS& smeared) {
  TGraph* graph(NULL);
  TGraph* iterGraph(NULL);  // For iterator position
  TIter iter(Graphs()->GetListOfGraphs());
  // Loop through graphs for each eta range until we find the one for
  // the eta of this particle
  while ((iterGraph = static_cast<TGraph*>(iter()))) {
    // Array of eta points in the graph
    double* eta = iterGraph->GetX();
    if (particle.GetEta() > eta[0] &&
        particle.GetEta() < eta[iterGraph->GetN() - 1] ) {
      graph = iterGraph;
      break;
    }  // if
  }  // while
  if (graph) {
    // Stored values are sigma(1/p).
    // We want to get sigma(p), which is equivalent to sigma(1/p) * p^2.
    double sigmaP = graph->Eval(particle.GetEta()) // sigma(1/p)
      * std::pow(particle.GetP(), 2.);
    // Add multiple scattering in quadrature with the linear term if requested
    if (mMultipleScattering) {
      sigmaP = std::sqrt(std::pow(sigmaP, 2.) +
                         std::pow(computeMultipleScattering(particle), 2.));
    }  // if
    double sigma = gRandom->Gaus(particle.GetP(), sigmaP);
    smeared.SetP(std::max(sigma, 0.));  // Store p, must not be negative
  }  // if
}

double EPhenixMomentum::computeMultipleScattering(
						  const erhic::VirtualParticle& particle) const {
  // Here is a duplicate of the relevant information from Sasha's email:
  //  "As for multiple scattering term, what has been simulated so far is eta
  //   ranges 2-3 and 3-4. Numbers keep changing after inclusion of different
  //   material (associated with detectors and read-out). The current very
  //   conservative advice is to use 3% for eta=2-3; in eta=3-4, log of this
  //   term changes ~linearly from 3% at eta=3 to ~15% at eta=4."
  // Note that the 3% etc is sigma(P) / P.
  const double eta = particle.GetEta();
  if (eta > 2. && eta <= 3.) {
    return 0.03 * particle.GetP();  // Flat sigma(P)/P = 3% in this range
  } else if (eta > 3. && eta <= 4.) {
    // The behaviour Sasha describes is (I think!)
    // log(M.S.) ~= log(3%) + (log(15%) - log(3%)) * (eta - 3) [3 < eta < 4]
    double logSigmaP = std::log(0.03) +
      (std::log(0.15) - std::log(0.03)) * (eta - 3);
    // Remember to multiply by momentum, as the 3-15% is sigma(P)/P and we
    // want sigma(P)
    return std::exp(logSigmaP) * particle.GetP();
  }  // if
  return 0.;  // Outside 2 < eta < 4
}

void EPhenixMomentum::Draw(Option_t* option) {
  Graphs()->Draw(option);
  if (gPad) {
    gPad->SetLogy(1);
  }  // if
}

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
Smear::Detector BuildEphoenix(bool multipleScattering) {
  EPhenixMomentum momentum(multipleScattering);
  // Define acceptance zones for different ePHENIX regions:
  // - electron-going: -4 < eta < -1
  // - barrel: -1 < eta < 1
  // - hadron-going: 1 < eta < 4
  // As the same calorimeter performance is used in the barrel and hadron-going
  // directions we define them as a single zone
  // Note etamin -> thetamax, etamax -> thetamin
  Smear::Acceptance::Zone electronDirection(EPhenixMomentum::etaToTheta(-1.), EPhenixMomentum::etaToTheta(-4.));
  Smear::Acceptance::Zone barrelAndHadronDirection(EPhenixMomentum::etaToTheta(4.),EPhenixMomentum::etaToTheta(-1.));

  // Electron-going electromagnetic calorimeter
  Smear::Device electronEcal("E", "0.015*sqrt(E) + 0.01*E",
                             Smear::kElectromagnetic);
  electronEcal.Accept.AddZone(electronDirection);
  // Barrel and hadron-going electromagnetic calorimeter
  Smear::Device barrelAndHadronEcal("E", "0.12*sqrt(E) + 0.02*E",
                                    Smear::kElectromagnetic);
  barrelAndHadronEcal.Accept.AddZone(barrelAndHadronDirection);
  // Barrel and hadron-going hadronic calorimeter
  Smear::Device barrelAndHadronHcal("E", "sqrt(E)", Smear::kHadronic);
  barrelAndHadronHcal.Accept.AddZone(barrelAndHadronDirection);
  // Assume perfect theta and phi performance i.e. momentum resolution
  // dominates over angular resolution
  Smear::Device theta("theta", "0");
  Smear::Device phi("phi", "0");
  // PID performance is unparameterised as of now
  Smear::PerfectID pid;
  // Combine the devices into a detector.
  Smear::Detector ephenix;
  ephenix.AddDevice(momentum);
  ephenix.AddDevice(electronEcal);
  ephenix.AddDevice(barrelAndHadronEcal);
  ephenix.AddDevice(barrelAndHadronHcal);
  ephenix.AddDevice(theta);
  ephenix.AddDevice(phi);
  ephenix.AddDevice(pid);
  ephenix.SetEventKinematicsCalculator("NM JB DA");
  return ephenix;
}

/*
  Here is the email from Sasha giving the calorimeter performances used above:

  From: 	Alexander Bazilevsky <shura@bnl.gov>
  Subject: 	Re: a favour
  Date: 	26 December, 2013 1:47:19 PM EST
  To: 	elke-caroline aschenauer <elke@bnl.gov>, Thomas P Burton <tpb@bnl.gov>

  Hello Elke, Thomas,

  Of course numbers are still drifting, but for now what you could use is roughly the following:

  EMCal sigma_E/E:
  e-going direction: 1.5%/sqrt(E) \oplus 1%
  h-going direction and barrel: 12%/sqrt(E) \oplus 2%

  HCal sigma_E/E:
  h-going direction and barrel: 100%/sqrt(E)

  Tracking sigma_p/p:
  is quite complicated: the linear term is in attached plot; we also calculated constant (multiple scattering) term from Geant, which appeared to vary from under 1% at eta~1 to 3% at eta~3 and to larger values (~10%) towards eta~4. I'll give you numerical parametrization one of the following days, as soon as I get them from the author of these calculations (Jin Huang).

  Regards,

  Sasha.

  --------------------------------------------------------------------------------

  Here is the email from Sasha giving the momentum resolution points used above:

  From: 	Alexander Bazilevsky <shura@bnl.gov>
  Subject: 	Re: a favour
  Date: 	2 January, 2014 11:01:39 AM EST
  To: 	elke-caroline aschenauer <elke@bnl.gov>, Thomas P Burton <tpb@bnl.gov>

  Hello Elke, Thomas, 

  Ok, the linear term in momentum resolution (as in the plot I sent you before): 

  eta        d(1/p) in (1/GeV)
  -3.0       0.010500
  -2.9       0.009820
  -2.8       0.009140
  -2.7       0.008460
  -2.6       0.007780
  -2.5       0.007100
  -2.4       0.006520
  -2.3       0.005940
  -2.2       0.005360
  -2.1       0.004780
  -2.0       0.004200

  -2.0       0.009811
  -1.9       0.008759
  -1.8       0.007479
  -1.7       0.006242
  -1.6       0.005436
  -1.5       0.004524

  -1.5       0.006887
  -1.4       0.005374
  -1.3       0.004485
  -1.2       0.003822
  -1.1       0.003164
  -1.0       0.002592
  -0.9       0.002791
  -0.8       0.002991
  -0.7       0.003187
  -0.6       0.003374
  -0.5       0.003547
  -0.4       0.003700
  -0.3       0.003827
  -0.2       0.003921
  -0.1       0.003980
  0.0       0.004000
  0.1       0.003980
  0.2       0.003921
  0.3       0.003827
  0.4       0.003700
  0.5       0.003547
  0.6       0.003374
  0.7       0.003187
  0.8       0.002991
  0.9       0.002791
  1.0       0.002592

  1.0       0.001200
  1.1       0.001280
  1.2       0.001360
  1.3       0.001440
  1.4       0.001520
  1.5       0.001600
  1.6       0.001840
  1.7       0.002080
  1.8       0.002320
  1.9       0.002560
  2.0       0.002800
  2.1       0.003160
  2.2       0.003520
  2.3       0.003880
  2.4       0.004240
  2.5       0.004600

  2.5       0.002300
  2.6       0.002620
  2.7       0.002940
  2.8       0.003260
  2.9       0.003580
  3.0       0.003900
  3.1       0.004400
  3.2       0.004900
  3.3       0.005400
  3.4       0.005900
  3.5       0.006400
  3.6       0.007220
  3.7       0.008040
  3.8       0.008860
  3.9       0.009680
  4.0       0.010500

  As for multiple scattering term, what has been simulated so far is eta ranges 2-3 and 3-4. Numbers keep changing after inclusion of different material (associated with detectors and read-out). The current very conservative advice is to use 3% for eta=2-3; in eta=3-4, log of this term changes ~linearly from 3% at eta=3 to ~15% at eta=4. 

  That's what we have for now. 

  Regards, 

  Sasha.
*/
