// Adapted from https://gitlab.com/eic/escalate/ejana/-/blob/master/src/plugins/reco/eic_smear/SS_DetectorJleic_v1_0_2.cc
// Work in progress - no guarantees for correctness!
// Characteristics:
// -- Calorimetry for CHARGED and PHOTONS: sigma (E) = a *sqrt(E)  (+)  b*E (+) c
// ----  in (Eta() > -3.5 && Eta() < -1.1): a= 0.02, b=0.001, c = 0.005
// ----  in (Eta() > -1.1 && Eta() <  3.5): a= 0.1,  b=0.01,  c = 0.005
// -- Calorimetry for NEUTRALS:
//     100% /  sqrt(p.E()) with NO acceptance cuts!
//     technically, special case with the same 100% /sqrt(p.E()) for neutrons in Theta() < 0.01
// -- Tracking for all CHARGED sigma_PT
// ----  in (Eta() > -3.5 && Eta() < 3.5): sigma = 0.01 Pt()^2 (+) 0.005 Pt()
//        with angular resolution  = 0.001 in phi and theta
// ----  in Theta() < 0.01 : sigma = 0.02
// ----  in 0.01 < Theta()< 0.05 : sigma = 0.05




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

using std::cout;
using std::endl;

// declare static --> local to this file, won't clash with others
static double ThetaFromEta( const double eta );
static double EtaFromTheta( const double theta );

Smear::Detector BuildJLEIC() {

  gSystem->Load("libeicsmear");

  // Create the detector object to hold all devices
  Smear::Detector det;
  det.SetEventKinematicsCalculator("NM JB DA"); // The detector will calculate event kinematics from smeared values

  // Phi and theta smearing for all particles
  // ----------------------------------------
  // Eic-smear rule: if it is smeared, it needs phi and theta information
  // Here: Finite resolution only given for charged particles in +/-3.5
  // But neutral hadrons are smeared everywhere, so endow them with perfect resolution

  // Neutral hadrons - assume perfect resolution
  Smear::Acceptance::Zone AngleZoneNeutralHadrons (ThetaFromEta(15.),ThetaFromEta(-15.));
  Smear::Device SmearThetaNeutralHadrons(Smear::kTheta, "0");
  SmearThetaNeutralHadrons.Accept.AddZone(AngleZoneNeutralHadrons);
  SmearThetaNeutralHadrons.Accept.SetGenre(Smear::kHadronic);
  SmearThetaNeutralHadrons.Accept.SetCharge(Smear::kNeutral);
  det.AddDevice(SmearThetaNeutralHadrons);

  Smear::Device SmearPhiNeutralHadrons(Smear::kPhi, "0");
  SmearPhiNeutralHadrons.Accept.AddZone(AngleZoneNeutralHadrons);
  SmearPhiNeutralHadrons.Accept.SetGenre(Smear::kHadronic);
  SmearPhiNeutralHadrons.Accept.SetCharge(Smear::kNeutral);
  det.AddDevice(SmearPhiNeutralHadrons);

  // Charged particles, barrel+endcap
  Smear::Acceptance::Zone AngleZoneChargedBarrel(ThetaFromEta ( 3.5 ),ThetaFromEta ( -3.5 ));
  Smear::Device SmearThetaChargedBarrel(Smear::kTheta, "0.001");
  SmearThetaChargedBarrel.Accept.AddZone(AngleZoneChargedBarrel);
  SmearThetaChargedBarrel.Accept.SetCharge(Smear::kCharged);
  det.AddDevice(SmearThetaChargedBarrel);

  Smear::Device SmearPhiChargedBarrel(Smear::kPhi, "0.001");
  SmearPhiChargedBarrel.Accept.AddZone(AngleZoneChargedBarrel);
  SmearPhiChargedBarrel.Accept.SetCharge(Smear::kCharged);
  det.AddDevice(SmearPhiChargedBarrel);

  // barrel + endcap calorimetry also accepts photons
  Smear::Device SmearThetaCalPhotons(Smear::kTheta, "0.001");
  SmearThetaCalPhotons.Accept.AddZone(AngleZoneChargedBarrel);
  SmearThetaCalPhotons.Accept.AddParticle(22);
  det.AddDevice(SmearThetaCalPhotons);

  Smear::Device SmearPhiCalPhotons(Smear::kPhi, "0.001");
  SmearPhiCalPhotons.Accept.AddZone(AngleZoneChargedBarrel);
  SmearPhiCalPhotons.Accept.AddParticle(22);
  det.AddDevice(SmearPhiCalPhotons);


  // Charged particles, far forward - assume perfect resolution
  // keep theta finite to avoid division by zero
  Smear::Acceptance::Zone AngleZoneChargedFarForward( 1e-12, 0.05);
  Smear::Device SmearThetaChargedFarForward(Smear::kTheta, "0");
  SmearThetaChargedFarForward.Accept.AddZone(AngleZoneChargedFarForward);
  SmearThetaChargedFarForward.Accept.SetCharge(Smear::kCharged);
  det.AddDevice(SmearThetaChargedFarForward);

  Smear::Device SmearPhiChargedFarForward(Smear::kPhi, "0");
  SmearPhiChargedFarForward.Accept.AddZone(AngleZoneChargedFarForward);
  SmearPhiChargedFarForward.Accept.SetCharge(Smear::kCharged);
  det.AddDevice(SmearPhiChargedFarForward);

  // Calorimetry, Neutral hadrons
  // ----------------------------
  // Not separating out the ZDC since it logically does nothing on top of this.
  Smear::Device CalNeutralHadrons(Smear::kE, "1*sqrt(E)" );
  CalNeutralHadrons.Accept.AddZone(AngleZoneNeutralHadrons);
  CalNeutralHadrons.Accept.SetGenre(Smear::kHadronic);
  CalNeutralHadrons.Accept.SetCharge(Smear::kNeutral);
  det.AddDevice( CalNeutralHadrons );

  // Calorimetry, Charged and photons
  // --------------------------------
  // There's no elegant way to accept charged hadrons, leptons, and photons
  // so we make two devices
  // eta = -3.5 -- -1.1
  // sigma^2 = (0.02*sqrt(E))^2 + (0.001*E)^2 + 0.005^2;
  Smear::Acceptance::Zone CalBackZone(ThetaFromEta ( -1.1 ), ThetaFromEta ( -3.5 ));
  Smear::Device CalBack(Smear::kE, "sqrt( pow(0.02*sqrt(E),2) + pow(0.001*E,2) + pow(0.005,2) )" );
  CalBack.Accept.AddZone(CalBackZone);
  CalBack.Accept.SetCharge(Smear::kCharged);
  CalBack.Accept.SetGenre(Smear::kHadronic);
  det.AddDevice(CalBack);

  // eta = -1.1 -- 3.5
  // sigma^2 = (0.1*sqrt(E))^2 + (0.01*E)^2 + 0.005^2;
  Smear::Acceptance::Zone CalForwardZone(ThetaFromEta ( 3.5 ), ThetaFromEta ( -1.1 ));
  Smear::Device CalForward(Smear::kE, "sqrt( pow(0.1*sqrt(E),2) + pow(0.01*E,2) + pow(0.005,2) )" );
  CalForward.Accept.AddZone(CalForwardZone);
  CalForward.Accept.SetCharge(Smear::kCharged);
  CalForward.Accept.SetGenre(Smear::kHadronic);
  det.AddDevice(CalForward);

  // eta = -3.5 -- -1.1
  // sigma^2 = (0.02*sqrt(E))^2 + (0.001*E)^2 + 0.005^2;
  Smear::Device PhotonCalBack(Smear::kE, "sqrt( pow(0.02*sqrt(E),2) + pow(0.001*E,2) + pow(0.005,2) )" );
  PhotonCalBack.Accept.AddZone(CalBackZone);
  PhotonCalBack.Accept.SetGenre(Smear::kElectromagnetic);
  det.AddDevice(PhotonCalBack);

  // eta = -1.1 -- 3.5
  // sigma^2 = (0.1*sqrt(E))^2 + (0.01*E)^2 + 0.005^2;
  Smear::Device PhotonCalForward(Smear::kE, "sqrt( pow(0.1*sqrt(E),2) + pow(0.01*E,2) + pow(0.005,2) )" );
  PhotonCalForward.Accept.AddZone(CalForwardZone);
  PhotonCalForward.Accept.SetGenre(Smear::kElectromagnetic);
  det.AddDevice(PhotonCalForward);

  // Tracking
  // --------
  // Note: Original implementation has vertex smearing,
  // which (at least for now) is not supported here

  // eta = -3.5 --  3.5
  Smear::Acceptance::Zone TrackBarrelZone(ThetaFromEta ( 3.5 ),ThetaFromEta ( -3.5 ));
  Smear::Device TrackBarrel(Smear::kPt, "sqrt ( pow( 0.01 * pow(pT,2), 2) + pow(0.005 * pT,2))");
  TrackBarrel.Accept.AddZone(TrackBarrelZone);
  TrackBarrel.Accept.SetCharge(Smear::kCharged);
  det.AddDevice(TrackBarrel);


  // Roman Pots, Theta < 0.01
  Smear::Acceptance::Zone TrackRPZone(1e-7,0.01);
  Smear::Device TrackRP(Smear::kPt, "0.02");
  TrackRP.Accept.AddZone(TrackRPZone);
  TrackRP.Accept.SetCharge(Smear::kCharged);
  det.AddDevice(TrackRP);

  // D1, 0.01 < Theta < 0.05
  Smear::Acceptance::Zone TrackD1Zone(0.01,0.05);
  Smear::Device TrackD1(Smear::kPt, "0.01");
  TrackD1.Accept.AddZone(TrackD1Zone);
  TrackD1.Accept.SetCharge(Smear::kCharged);
  det.AddDevice(TrackD1);

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
