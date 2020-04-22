// Based on
// Electron-Ion Collider Detector Requirements and R&D Handbook
// Version 1.2 February 25, 2020
// Page 26
// Retrieved from
// http://www.eicug.org/web/content/detector-rd
// http://www.eicug.org/web/sites/default/files/EIC_HANDBOOK_v1.2.pdf

// Additional details from discussions within the calorimetry working group:
// 
// >    These are two configuration for hadron endcap we talked today to consider.
// > Base configuration:
// > EM - W/ScFi, granularity 2.5 cm x 2.5 cm, 12% stochastic, 2% constant terms
// > HCal- Fe/Sc, granularity 10 x 10 cm, 50% stochastic, 10% constant

// > Configuration 2:
// > EM is the same as above.
// > HCal -Pb/Sc, 10 x 10 cm, 45% stochastic, 6% constant terms.

// > Alternative ECAL
// > Shashlyk, Pb/Sc 5.5 x5.5, 8% stochastic, 2% constant.

// > Thanks,
// > Oleg [Tsai]


#include "eicsmear/erhic/VirtualParticle.h"
#include "eicsmear/smear/Acceptance.h"
#include "eicsmear/smear/Device.h"
#include "eicsmear/smear/Detector.h"
#include "eicsmear/smear/Smearer.h"
#include "eicsmear/smear/ParticleMCS.h"
#include "eicsmear/smear/PerfectID.h"
#include <eicsmear/smear/Smear.h>
#include <eicsmear/erhic/ParticleMC.h>
#include "Math/Vector4D.h"

// declare static --> local to this file, won't clash with others
static double ThetaFromEta( const double eta );
static double EtaFromTheta( const double theta );

Smear::Detector BuildHandBookDetector() {

  gSystem->Load("libeicsmear");

  // Create the detector object to hold all devices
  Smear::Detector det;

  // The framework provides implementations of three kinematic calculation methods
  // from smeared values
  // NM - "Null method" uses the scattered lepton.
  // DA - Double Angle method
  // JB - Jacquet-Blondel method
  // Important Notes:
  // - All methods rely on measured energy and momentum components (lepton for NM, hadrons for DA, JB).
  //   In order to rely on these methods, you therefore _need_ to cover the relevant
  //   regions with smearers for momentum, angles and energy.
  // - This is exacerbated by the fact that in the current implementation a smearing of e.g. P
  //   may fool the framework into assuming theta is measured to be the initialization value 0,
  //   leading to nan or inf or otherwise wrong calculations.
  //   NM catches P=0 and replaces it in calculations with E,
  //   but JB and DA either don't work at all or report nonsenical values.
  // - It may be advantageous to use a P measurement in place of E in a more general case.
  //   In the future, we plan to change to a weighted mean approach.
  // The summary of the above is that the user should be mindful of the limitations and assumptions in
  // the calculation of smeared kinematics. Sophisticated calculations are best redone in user code
  // from the smeared particles directly.
  det.SetEventKinematicsCalculator("NM DA JB");

    
  // IMPORTANT: There are two traps (will be addressed in future releases):
  //            1) If you smear E but don't provide phi, theta smearers, those values will be
  //               set to 0, not to a fault value and not to the truth level
  //            2) If you do provide such a smearer, pt and pz will be changed
  //               by a consistency enforcer in Detector::Smear()


  // As detailed above, we always need phi and theta smearers.
  // In the absence of handbook values, assume 0.001 radians throughout the acceptance
  // That's a reasonably conservative guesstimate for tracking,
  // and probably a bit too good for most calorimeters.
  // In reality, you would have granularity and clustering assumptions
  // which are best handled in an afterburner (like tracking efficiency).

  // Phi and theta smearing for all particles
  // ----------------------------------------
  // total coverage of the handbook for tracker and hcal is -3.5 < eta < 3.5
  Smear::Acceptance::Zone AngleZoneHadronic(ThetaFromEta ( 3.5 ),ThetaFromEta ( -3.5 ));
  Smear::Device SmearThetaHadronic(Smear::kTheta, "0.001");
  SmearThetaHadronic.Accept.AddZone(AngleZoneHadronic);
  SmearThetaHadronic.Accept.SetGenre(Smear::kHadronic);    
  det.AddDevice(SmearThetaHadronic);
  
  Smear::Device SmearPhiHadronic(Smear::kPhi, "0.001");
  SmearPhiHadronic.Accept.AddZone(AngleZoneHadronic);
  SmearPhiHadronic.Accept.SetGenre(Smear::kHadronic);
  det.AddDevice(SmearPhiHadronic);
  
  // emcal stretches to -4.5 < eta < 4.5
  Smear::Acceptance::Zone AngleZoneEmcal(ThetaFromEta ( 4.5 ),ThetaFromEta ( -4.5 ));
  Smear::Device SmearThetaEmcal(Smear::kTheta, "0.001");
  SmearThetaEmcal.Accept.AddZone(AngleZoneEmcal);
  SmearThetaEmcal.Accept.SetGenre(Smear::kElectromagnetic);    
  det.AddDevice(SmearThetaEmcal);
  
  Smear::Device SmearPhiEmcal(Smear::kPhi, "0.001");
  SmearPhiEmcal.Accept.AddZone(AngleZoneEmcal);
  SmearPhiEmcal.Accept.SetGenre(Smear::kElectromagnetic);
  det.AddDevice(SmearPhiEmcal);

  // Tracking
  // --------
  // Note: Smear::kCharged checks pdg charge, so includes muons (good)
  // eta = -3.5 --  -2.5
  // sigma_p/p ~ 0.1% p+2.0%
  Smear::Acceptance::Zone TrackBack1Zone(ThetaFromEta ( -2.5 ),ThetaFromEta ( -3.5 ));
  Smear::Device TrackBack1P(Smear::kP, "sqrt( pow ( 0.001*P*P, 2) + pow ( 0.02*P, 2) )");
  TrackBack1P.Accept.AddZone(TrackBack1Zone);
  TrackBack1P.Accept.SetCharge(Smear::kCharged);
  // TrackBack1P.Accept.SetGenre(Smear::kHadronic);
  det.AddDevice(TrackBack1P);

  // eta = -2.5 --  -1
  // sigma_p/p ~ 0.05% p+1.0%
  Smear::Acceptance::Zone TrackBack2Zone(ThetaFromEta ( -1 ),ThetaFromEta ( -2.5 ));
  Smear::Device TrackBack2P(Smear::kP, "sqrt( pow ( 0.0005*P*P, 2) + pow ( 0.01*P, 2) )");
  TrackBack2P.Accept.AddZone(TrackBack2Zone);
  TrackBack2P.Accept.SetCharge(Smear::kCharged);
  // TrackBack2P.Accept.SetGenre(Smear::kHadronic);
  det.AddDevice(TrackBack2P);

  // eta = -1 -- +1
  // sigma_p/p ~ 0.05% p+0.5%
  Smear::Acceptance::Zone TrackBarrelZone(ThetaFromEta ( 1 ),ThetaFromEta ( -1 ));
  Smear::Device TrackBarrelP(Smear::kP, "sqrt( pow ( 0.0005*P*P, 2) + pow ( 0.005*P, 2) )");
  TrackBarrelP.Accept.AddZone(TrackBarrelZone);
  TrackBarrelP.Accept.SetCharge(Smear::kCharged);
  // TrackBarrelP.Accept.SetGenre(Smear::kHadronic);
  det.AddDevice(TrackBarrelP);

  // eta = 1 -- 2.5
  // sigma_p/p ~ 0.05% p+1.0%
  Smear::Acceptance::Zone TrackFwd2Zone(ThetaFromEta ( 2.5 ),ThetaFromEta ( 1 ));
  Smear::Device TrackFwd2P(Smear::kP, "sqrt( pow ( 0.0005*P*P, 2) + pow ( 0.01*P, 2) )");
  TrackFwd2P.Accept.AddZone(TrackFwd2Zone);
  TrackFwd2P.Accept.SetCharge(Smear::kCharged);
  // TrackFwd2P.Accept.SetGenre(Smear::kHadronic);
  det.AddDevice(TrackFwd2P);

  // eta = 2.5 -- 3.5
  // sigma_p/p ~ 0.1% p+2.0%
  Smear::Acceptance::Zone TrackFwd1Zone(ThetaFromEta ( 3.5 ),ThetaFromEta ( 2.5 ));
  Smear::Device TrackFwd1P(Smear::kP, "sqrt( pow ( 0.001*P*P, 2) + pow ( 0.02*P, 2) )");
  TrackFwd1P.Accept.AddZone(TrackFwd1Zone);
  TrackFwd1P.Accept.SetCharge(Smear::kCharged);
  // TrackFwd1P.Accept.SetGenre(Smear::kHadronic);
  det.AddDevice(TrackFwd1P);

  // TODO: Add low-Q^2 tagger: -6.9<eta<-5.8: Delta_theta/theta < 1.5%; 10^-6 < Q2 < 10^-2 GeV2
  // TODO: Add proton spectrometer:  eta>6.2: sigma_intrinsic(|t|)/|t| < 1%; Acceptance: 0.2 < pT < 1.2 GeV/c
  // TODO: Add material budget X/X0 <~5%
  // TODO: Add vertexing: sigma_xyz ~ 20 microns, d0(z) ~ d0(r phi) ~ (20 microns)/(pT [GeV])  + 5 microns
  
  // EM Calorimeters
  // ---------------
  // Note: Smear::kElectromagnetic == gamma + e. Does not include muons (good)
  
  // Calorimeter resolution usually given as sigma_E/E = const% + stocastic%/Sqrt{E}
  // EIC Smear needs absolute sigma: sigma_E = Sqrt{const*const*E*E + stoc*stoc*E}

  // IMPORTANT: The Handbook table 2 does not provide constant terms,
  //            which may be important for some measurements.
  //            We adapted the handbook here to values provided by the calo DWG where available
  
  // Back
  // eta = -4.5 -- -2
  // stoch. = 2%
  // Using const. = 1% (PANDA-ish)
  // Note: p. 41f calls for stoch = 1--1.5%, const <~ 0.5%
  Smear::Acceptance::Zone EmcalBackZone(ThetaFromEta ( -2 ),ThetaFromEta ( -4.5 ));
  Smear::Device EmcalBack(Smear::kE, "sqrt( pow ( 0.01*E,2 ) + pow( 0.01,2)*E)");
  EmcalBack.Accept.AddZone(EmcalBackZone);
  EmcalBack.Accept.SetGenre(Smear::kElectromagnetic);
  det.AddDevice(EmcalBack);

  // MidBack
  // eta = -2 -- -1
  // stoch. = 7%
  // Could use const. = 1% (bit worse than KLOE)
  // but using the Calo group's alternative instead, Shashlyk, Pb/Sc 8% stochastic, 2% constant.
  Smear::Acceptance::Zone EmcalMidBackZone(ThetaFromEta ( -1 ),ThetaFromEta ( -2 ));
  // Smear::Device EmcalMidBack(Smear::kE, "sqrt( pow ( 0.01*E,2 ) + pow( 0.07,2)*E)");
  Smear::Device EmcalMidBack(Smear::kE, "sqrt( pow ( 0.02*E,2 ) + pow( 0.08,2)*E)");
  EmcalMidBack.Accept.AddZone(EmcalMidBackZone);
  EmcalMidBack.Accept.SetGenre(Smear::kElectromagnetic);
  det.AddDevice(EmcalMidBack);
  
  // Forward
  // eta = -1 -- 4.5
  // stoch. = 10-12%.
  // From calo group: using 12% and const.=2% (W/SciFi)
  Smear::Acceptance::Zone EmcalFwdZone(ThetaFromEta ( 4.5 ),ThetaFromEta ( -1 ));
  Smear::Device EmcalFwd(Smear::kE, "sqrt( pow ( 0.02*E,2 ) + pow( 0.12,2)*E)");
  EmcalFwd.Accept.AddZone(EmcalFwdZone);
  EmcalFwd.Accept.SetGenre(Smear::kElectromagnetic);
  det.AddDevice(EmcalFwd);

  // TODO: Add PID
  // Could turn on perfect PID
  // Smear::Acceptance::Zone acceptall(etaToTheta(15.),etaToTheta(-15.));
  // Smear::PerfectID pid;
  // pid.Accept.AddZone(acceptall);
  // det.AddDevice( pid );
  
  // Hadronic  Calorimeters
  // ----------------------
  // IMPORTANT: The Handbook table 2 does not provide constant terms,
  //            which may be important for some measurements.
  //            The barrel options modeled after existing instruments will have one,
  //            leaving it at 0 otherwise
  // Note: kHadronic == |pdg|>110. 
    
  // Back
  // eta = -3.5 -- -1
  // stoch. = 50%
  // Using configuration 2 from calo group, Pb/Sc 45% stochastic, 6% constant terms.
  Smear::Acceptance::Zone HcalBackZone(ThetaFromEta ( -1 ),ThetaFromEta ( -3.5 ));
  Smear::Device HcalBack(Smear::kE, "sqrt(pow( 0.06*E, 2) + pow ( 0.45,2) *E)");
  HcalBack.Accept.AddZone(HcalBackZone);
  HcalBack.Accept.SetGenre(Smear::kHadronic);
  det.AddDevice(HcalBack);
  
  // Barrel
  // Handbook: TBD
  // Providing two options
  // ~CMS
  Smear::Device HcalBarrel(Smear::kE, "sqrt( pow( 0.07*E, 2)  + pow( 0.85, 2)*E)");
  // Zeus
  // Smear::Device HcalBarrel(Smear::kE, "sqrt( pow( 0.02*E, 2) + pow( 0.35,2) *E)");

  // eta = -1 -- 1
  Smear::Acceptance::Zone HcalBarrelZone(ThetaFromEta ( 1 ),ThetaFromEta ( -1 ));
  HcalBarrel.Accept.AddZone(HcalBarrelZone);
  HcalBarrel.Accept.SetGenre(Smear::kHadronic);
  det.AddDevice(HcalBarrel);

  // Forward
  // eta = 1 -- 3.5
  // stoch. = 50%
  // Using configuration 2 from calo group, Pb/Sc 45% stochastic, 6% constant terms.
  Smear::Acceptance::Zone HcalFwdZone(ThetaFromEta ( 3.5 ),ThetaFromEta ( 1 ));
  Smear::Device HcalFwd(Smear::kE, "sqrt(pow( 0.06*E, 2) + pow ( 0.45,2) *E)");
  HcalFwd.Accept.AddZone(HcalFwdZone);
  HcalFwd.Accept.SetGenre(Smear::kHadronic);
  det.AddDevice(HcalFwd);

  return det;
}

// -------------------------------------------------------------------
double ThetaFromEta( const double eta ) {
  if ( !std::isnan(eta) && !std::isinf(eta)   ) {
    return 2.0 * atan( exp( -eta ));
  }
  throw std::runtime_error("ThetaFromEta called with NaN or Inf");
  return -1;
}

// -------------------------------------------------------------------
double EtaFromTheta( const double theta ) {
  // The default value of -19 is used in the main eRHIC code,
  // so use that for consistency.
  double eta = -19.;
  if (theta > 0. && theta < TMath::Pi() && !std::isnan(theta) && !std::isinf(theta)) {
    eta = -log(tan(theta / 2.));
  }
  return eta;
}
// -------------------------------------------------------------------

