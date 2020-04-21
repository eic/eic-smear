// NOT intended tobe representative of the most recent
// BeAST developments (coming soon)
// Rather an example, similar to what was used for the design study,
// using more complicated parameterizations than the Handbook

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


Smear::Detector BuildBeAST() {

  gSystem->Load("libeicsmear");

  // Create the detector object to hold all devices
  Smear::Detector det;
  det.SetEventKinematicsCalculator("NM JB DA"); // The detector will calculate event kinematics from smeared values

  // For proper kinematics calculations,
  // phi and theta need to be smeared (at least with "0") for
  // all particles that see a kE or kP smearer

  // Phi and theta smearing for all particles
  // ----------------------------------------
  // eta = -3.5 -- 3.5
  // charged particles in the tracker
  // includes electrons
  Smear::Acceptance::Zone TrackZone(ThetaFromEta ( 3.5 ),ThetaFromEta ( -3.5 ));
  Smear::Device trackTheta(Smear::kTheta, "((1.0/(1.0*P))*(0.752935 + 0.280370*pow((-log(tan(theta/2.0))), 2) - 0.0359713*pow((-log(tan(theta/2.0))), 4) + 0.00200623*pow((-log(tan(theta/2.0))), 6)) + 0.0282315 - 0.00998623*pow((-log(tan(theta/2.0))), 2) + 0.00117487*pow((-log(tan(theta/2.0))), 4) - 0.0000443918*pow((-log(tan(theta/2.0))), 6))*0.001");
  trackTheta.Accept.AddZone(TrackZone);
  trackTheta.Accept.SetCharge(Smear::kCharged);
  det.AddDevice(trackTheta);
  Smear::Device trackPhi(Smear::kPhi, "((1.0/(1.0*P))*(0.743977 + 0.753393*pow((-log(tan(theta/2.0))), 2) + 0.0634184*pow((-log(tan(theta/2.0))), 4) + 0.0128001*pow((-log(tan(theta/2.0))), 6)) + 0.0308753 + 0.0480770*pow((-log(tan(theta/2.0))), 2) - 0.0129859*pow((-log(tan(theta/2.0))), 4) + 0.00109374*pow((-log(tan(theta/2.0))), 6))*0.001");
  trackPhi.Accept.AddZone(TrackZone);
  trackPhi.Accept.SetCharge(Smear::kCharged);
  det.AddDevice(trackPhi);

  // emcal stretches to -4.5 < eta < 4.5
  // includes gammas
  // Assuming 1 mrad
  Smear::Acceptance::Zone AngleZoneEmcalB(ThetaFromEta ( -3.5 ),ThetaFromEta ( -4.5 ));
  Smear::Device SmearThetaEmcalB(Smear::kTheta, "0.001");
  SmearThetaEmcalB.Accept.AddZone(AngleZoneEmcalB);
  SmearThetaEmcalB.Accept.SetGenre(Smear::kElectromagnetic);
  det.AddDevice(SmearThetaEmcalB);

  Smear::Device SmearPhiEmcalB(Smear::kPhi, "0.001");
  SmearPhiEmcalB.Accept.AddZone(AngleZoneEmcalB);
  SmearPhiEmcalB.Accept.SetGenre(Smear::kElectromagnetic);
  det.AddDevice(SmearPhiEmcalB);

  Smear::Acceptance::Zone AngleZoneEmcalF(ThetaFromEta ( 4.5 ),ThetaFromEta ( 3.5 ));
  Smear::Device SmearThetaEmcalF(Smear::kTheta, "0.001");
  SmearThetaEmcalF.Accept.AddZone(AngleZoneEmcalF);
  SmearThetaEmcalF.Accept.SetGenre(Smear::kElectromagnetic);
  det.AddDevice(SmearThetaEmcalF);

  Smear::Device SmearPhiEmcalF(Smear::kPhi, "0.001");
  SmearPhiEmcalF.Accept.AddZone(AngleZoneEmcalF);
  SmearPhiEmcalF.Accept.SetGenre(Smear::kElectromagnetic);
  det.AddDevice(SmearPhiEmcalF);

  // Haven't covered gammas in the barrel yet
  Smear::Acceptance::Zone AngleZoneGamma(ThetaFromEta ( 3.5 ),ThetaFromEta ( -3.5 ));
  Smear::Device SmearThetaGamma(Smear::kTheta, "0.001");
  SmearThetaGamma.Accept.AddZone(AngleZoneGamma);
  SmearThetaGamma.Accept.AddParticle(22);
  det.AddDevice(SmearThetaGamma);

  Smear::Device SmearPhiGamma(Smear::kPhi, "0.001");
  SmearPhiGamma.Accept.AddZone(AngleZoneGamma);
  SmearPhiGamma.Accept.AddParticle(22);
  det.AddDevice(SmearPhiGamma);

  // Wha'ts left is neutral hadrons.
  // HCal in this example covers +/-4.5
  // Using the same, very optimistic, 1 mrad
  Smear::Acceptance::Zone AngleZoneHCal(ThetaFromEta ( 4.5 ),ThetaFromEta ( -4.5 ));
  Smear::Device SmearThetaHCal(Smear::kTheta, "0.001");
  SmearThetaHCal.Accept.AddZone(AngleZoneHCal);
  SmearThetaHCal.Accept.SetGenre(Smear::kHadronic);
  SmearThetaHCal.Accept.SetCharge(Smear::kNeutral);
  det.AddDevice(SmearThetaHCal);

  Smear::Device SmearPhiHCal(Smear::kPhi, "0.001");
  SmearPhiHCal.Accept.AddZone(AngleZoneHCal);
  SmearPhiHCal.Accept.SetGenre(Smear::kHadronic);
  SmearPhiHCal.Accept.SetCharge(Smear::kNeutral);
  det.AddDevice(SmearPhiHCal);

  // Tracking
  // --------
  // Note: Smear::kCharged checks pdg charge, so includes muons (good)

  // The momentum parametrization (a*p + b) gives sigma_P/P in percent.
  // So Multiply through by P and divide by 100 to get absolute sigma_P
  // Theta and Phi parametrizations give absolute sigma in miliradians

  // Track Momentum
  // eta = -3.5 -- 3.5
  // sigma_p/p ~ 0.1% p+2.0%
  Smear::Device momentum(Smear::kP, "(P*P*(0.0182031 + 0.00921047*pow((-log(tan(theta/2.0))), 2) - 0.00291243*pow((-log(tan(theta/2.0))), 4) + 0.000264353*pow((-log(tan(theta/2.0))), 6)) + P*(0.209681 + 0.275144*pow((-log(tan(theta/2.0))), 2) - 0.0436536*pow((-log(tan(theta/2.0))), 4) + 0.00367412*pow((-log(tan(theta/2.0))), 6)))*0.01");
  momentum.Accept.SetCharge(Smear::kCharged);
  momentum.Accept.AddZone(TrackZone);
  det.AddDevice(momentum);

  // EMCal
  // -----
  // Calorimeter resolution usually given as sigma_E/E = const% + stocastic%/Sqrt{E}
  // EIC Smear needs absolute sigma: sigma_E = Sqrt{const*const*E*E + stoc*stoc*E}
  // eta = -4.5 -- -2
   Smear::Acceptance::Zone emBck(ThetaFromEta ( -2 ),ThetaFromEta ( -4.5 ));
   Smear::Device emcalBck(Smear::kE, "sqrt(0.01*0.01*E*E + 0.015*0.015*E)");
   emcalBck.Accept.AddZone(emBck);
   emcalBck.Accept.SetGenre(Smear::kElectromagnetic);
   det.AddDevice(emcalBck);
   
   // eta = -2 -- -1
   Smear::Acceptance::Zone emMidBck(ThetaFromEta ( -1 ),ThetaFromEta ( -2 ));
   Smear::Device emcalMidBck(Smear::kE, "sqrt(0.01*0.01*E*E + 0.07*0.07*E)");
   emcalMidBck.Accept.AddZone(emMidBck);
   emcalMidBck.Accept.SetGenre(Smear::kElectromagnetic);
   det.AddDevice(emcalMidBck);

   // eta = -1 -- 1
   Smear::Acceptance::Zone emMid(ThetaFromEta ( 1 ),ThetaFromEta ( -1 ));
   Smear::Device emcalMid(Smear::kE, "sqrt(0.01*0.01*E*E + 0.10*0.10*E)");
   emcalMid.Accept.AddZone(emMid);
   emcalMid.Accept.SetGenre(Smear::kElectromagnetic);
   det.AddDevice(emcalMid);

   // eta = 1 -- 4.5
   Smear::Acceptance::Zone emFwd(ThetaFromEta ( 4.5 ),ThetaFromEta ( 1 ));
   Smear::Device emcalFwd(Smear::kE, "sqrt(0.01*0.01*E*E + 0.07*0.07*E)");
   emcalFwd.Accept.AddZone(emFwd);
   emcalFwd.Accept.SetGenre(Smear::kElectromagnetic);
   det.AddDevice(emcalFwd);

   

  // HCal
  // -----
   // Create the Forward/Backward Hadron Calorimeter
   // eta = -4.5 -- -1
   Smear::Acceptance::Zone hBck(ThetaFromEta ( -1 ),ThetaFromEta ( -4.5 ));
   Smear::Device hcalBck(Smear::kE, "sqrt(0.015*0.015*E*E + 0.50*0.50*E)");
   hcalBck.Accept.AddZone(hBck);
   hcalBck.Accept.SetGenre(Smear::kHadronic);
   det.AddDevice(hcalBck);
		    
   // eta = 1 -- 4.5
   Smear::Acceptance::Zone hFwd(ThetaFromEta ( 4.5 ),ThetaFromEta ( 1 ));
   Smear::Device hcalFwd(Smear::kE, "sqrt(0.015*0.015*E*E + 0.50*0.50*E)");
   hcalFwd.Accept.AddZone(hFwd);
   hcalFwd.Accept.SetGenre(Smear::kHadronic);
   det.AddDevice(hcalFwd);
   
   // Create the Hypothetical Mid Rap Calorimeter
   // eta = -1 -- 1
   Smear::Acceptance::Zone hMid(ThetaFromEta ( 1 ),ThetaFromEta ( -1 ));
   //Smear::Device hcalMid(Smear::kE, "sqrt(0.07*0.07*E*E + 0.85*0.85*E)"); // ~CMS
   Smear::Device hcalMid(Smear::kE, "sqrt(0.02*0.02*E*E + 0.35*0.35*E)"); // ~Zeus
   hcalMid.Accept.AddZone(hMid);
   hcalMid.Accept.SetGenre(Smear::kHadronic);
   det.AddDevice(hcalMid);


   // // Create PID based on Hermes RICH.
   // Smear::ParticleID rich("PIDMatrix.dat");
   // det.AddDevice(rich);

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
// -------------------------------------------------------------------

