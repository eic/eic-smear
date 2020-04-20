// Based on
// Electron-Ion Collider Detector Requirements and R&D Handbook
// Version 1.2 February 25, 2020
// Page 26
// Retrieved from
// http://www.eicug.org/web/content/detector-rd
// http://www.eicug.org/web/sites/default/files/EIC_HANDBOOK_v1.2.pdf


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

Smear::Detector BuildNaiveHandBookDetector() {

  gSystem->Load("libeicsmear");

  // Create the detector object to hold all devices
  Smear::Detector det;
  // The detector will calculate event kinematics from smeared values
  det.SetEventKinematicsCalculator("NM JB DA");
    
  // Tracking
  // --------
  // Note: Smear::kCharged checks pdg charge, so includes muons (good)
  // Note: For kinematic calculations, both E and P of an electron
  //       will enter. Without intervention, the smearer will faithfully
  //       smear both and use wildly smeared values where a human analyzer
  //       would adapt the method, e.g. to a weighted one.
  //       In this implementation, we will ignore that for educational purposes.
  //       smearHandBook.cxx will attempt a solution.  

  // eta = -3.5 --  -2.5
  // sigma_p/p ~ 0.1% p+2.0%
  Smear::Acceptance::Zone TrackBack1Zone(ThetaFromEta ( -2.5 ),ThetaFromEta ( -3.5 ));
  Smear::Device TrackBack1P(Smear::kP, "sqrt( pow ( 0.001*P*P, 2) + pow ( 0.02*P, 2) )");
  TrackBack1P.Accept.AddZone(TrackBack1Zone);
  TrackBack1P.Accept.SetCharge(Smear::kCharged);
  det.AddDevice(TrackBack1P);

  // eta = -2.5 --  -1
  // sigma_p/p ~ 0.05% p+1.0%
  Smear::Acceptance::Zone TrackBack2Zone(ThetaFromEta ( -1 ),ThetaFromEta ( -2.5 ));
  Smear::Device TrackBack2P(Smear::kP, "sqrt( pow ( 0.0005*P*P, 2) + pow ( 0.01*P, 2) )");
  TrackBack2P.Accept.AddZone(TrackBack2Zone);
  TrackBack2P.Accept.SetCharge(Smear::kCharged);
  det.AddDevice(TrackBack2P);

  // eta = -1 -- +1
  // sigma_p/p ~ 0.05% p+0.5%
  Smear::Acceptance::Zone TrackBarrelZone(ThetaFromEta ( 1 ),ThetaFromEta ( -1 ));
  Smear::Device TrackBarrelP(Smear::kP, "sqrt( pow ( 0.0005*P*P, 2) + pow ( 0.005*P, 2) )");
  TrackBarrelP.Accept.AddZone(TrackBarrelZone);
  TrackBarrelP.Accept.SetCharge(Smear::kCharged);
  det.AddDevice(TrackBarrelP);

  // eta = 1 -- 2.5
  // sigma_p/p ~ 0.05% p+1.0%
  Smear::Acceptance::Zone TrackFwd2Zone(ThetaFromEta ( 2.5 ),ThetaFromEta ( 1 ));
  Smear::Device TrackFwd2P(Smear::kP, "sqrt( pow ( 0.0005*P*P, 2) + pow ( 0.01*P, 2) )");
  TrackFwd2P.Accept.AddZone(TrackFwd2Zone);
  TrackFwd2P.Accept.SetCharge(Smear::kCharged);
  det.AddDevice(TrackFwd2P);

  // eta = 2.5 -- 3.5
  // sigma_p/p ~ 0.1% p+2.0%
  Smear::Acceptance::Zone TrackFwd1Zone(ThetaFromEta ( 3.5 ),ThetaFromEta ( 2.5 ));
  Smear::Device TrackFwd1P(Smear::kP, "sqrt( pow ( 0.001*P*P, 2) + pow ( 0.02*P, 2) )");
  TrackFwd1P.Accept.AddZone(TrackFwd1Zone);
  TrackFwd1P.Accept.SetCharge(Smear::kCharged);
  det.AddDevice(TrackFwd1P);

  // Assume sigma = 0.001 in phi and theta throughout tracking
  Smear::Acceptance::Zone TrackAngleZone(ThetaFromEta ( 3.5 ),ThetaFromEta ( -3.5 ));
  Smear::Device TrackTheta(Smear::kTheta, "0.001");
  TrackTheta.Accept.AddZone(TrackAngleZone);
  TrackTheta.Accept.SetCharge(Smear::kCharged);
  det.AddDevice(TrackTheta);
  
  Smear::Device TrackPhi(Smear::kPhi, "0.001");
  TrackPhi.Accept.AddZone(TrackAngleZone);
  TrackPhi.Accept.SetCharge(Smear::kCharged);
  det.AddDevice(TrackPhi);

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
  //            I added a "reasonable" constant terms from tables 4 and 5
  
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
  // Using const. = 1% (bit worse than KLOE)
  Smear::Acceptance::Zone EmcalMidBackZone(ThetaFromEta ( -1 ),ThetaFromEta ( -2 ));
  Smear::Device EmcalMidBack(Smear::kE, "sqrt( pow ( 0.01*E,2 ) + pow( 0.07,2)*E)");
  EmcalMidBack.Accept.AddZone(EmcalMidBackZone);
  EmcalMidBack.Accept.SetGenre(Smear::kElectromagnetic);
  det.AddDevice(EmcalMidBack);
  
  // Forward
  // eta = -1 -- 4.5
  // stoch. = 10-12%. Using 12% and letting const.=0
  Smear::Acceptance::Zone EmcalFwdZone(ThetaFromEta ( 4.5 ),ThetaFromEta ( -1 ));
  Smear::Device EmcalFwd(Smear::kE, "sqrt( pow( 0.12,2)*E)");
  EmcalFwd.Accept.AddZone(EmcalFwdZone);
  EmcalFwd.Accept.SetGenre(Smear::kElectromagnetic);
  det.AddDevice(EmcalFwd);

  // TODO: Add PID
  
  // Hadronic  Calorimeters
  // ----------------------
  // IMPORTANT: The Handbook table 2 does not provide constant terms,
  //            which may be important for some measurements.
  //            The barrel options modeled after existing instruments will have one,
  //            leaving it at 0 otherwise
  // Note: kHadronic = |pdg|>110. 
    
  // Back
  // eta = -3.5 -- -1
  // stoch. = 50%
  Smear::Acceptance::Zone HcalBackZone(ThetaFromEta ( -1 ),ThetaFromEta ( -3.5 ));
  Smear::Device HcalBack(Smear::kE, "sqrt(0.50*0.50*E)");
  HcalBack.Accept.AddZone(HcalBackZone);
  HcalBack.Accept.SetGenre(Smear::kHadronic);
  det.AddDevice(HcalBack);
  
  // Barrel
  // Handbook: TBD
  // Providing two options
  // ~CMS
  Smear::Device HcalBarrel(Smear::kE, "sqrt(0.07*0.07*E*E + 0.85*0.85*E)");
  // Zeus
  // Smear::Device HcalBarrel(Smear::kE, "sqrt(0.02*0.02*E*E + 0.35*0.35*E)");

  // eta = -1 -- 1
  Smear::Acceptance::Zone HcalBarrelZone(ThetaFromEta ( 1 ),ThetaFromEta ( -1 ));
  HcalBarrel.Accept.AddZone(HcalBarrelZone);
  HcalBarrel.Accept.SetGenre(Smear::kHadronic);
  det.AddDevice(HcalBarrel);

  // Forward
  // eta = 1 -- 3.5
  // stoch. = 50%
  Smear::Acceptance::Zone HcalFwdZone(ThetaFromEta ( 3.5 ),ThetaFromEta ( 1 ));
  Smear::Device HcalFwd(Smear::kE, "sqrt(0.50*0.50*E)");
  HcalFwd.Accept.AddZone(HcalFwdZone);
  HcalFwd.Accept.SetGenre(Smear::kHadronic);
  det.AddDevice(HcalFwd);

  return det;
}

// -------------------------------------------------------------------
double ThetaFromEta( const double eta ) {
  if ( !isnan(eta) && !isinf(eta)   ) {
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
  if (theta > 0. && theta < TMath::Pi() && !isnan(theta) && !isinf(theta)) {
    eta = -log(tan(theta / 2.));
  }
  return eta;
}


