/**
 \file eSTARDetector.cpp
 Example smearing script for the eSTAR detector

 \author    Thomas Burton
 \date      2014-03-21
 \copyright 2014 Brookhaven National Lab
 */

/**
 Helper function to convert eta to theta (radians)
 
 Detector acceptances require theta, not eta
 */
double etaToTheta(double eta) {
  return 2. * std::atan(std::exp(-eta));
}

/**
 Helper function producing a zone in eta
 */
Smear::Acceptance::Zone zoneEta(double etamin, double etamax) {
  // First two arguments to Zone constructor are (theta-min, theta-max)
  // Note that eta-min --> theta-max and eta-max --> theta-min
  return Smear::Acceptance::Zone(etaToTheta(etamax), etaToTheta(etamin));
}

/**
 Smearing parameterisations for the eSTAR detector.

 Based on parameterisations given in Zhangbu Xu's talk here (slide 5):
 https://wiki.bnl.gov/conferences/index.php/January_2014
 This includes momentum and energy resolutions, but no particle identification.
 
 Note: you must gSystem->Load("libeicsmear") BEFORE loading this script,
 as ROOT needs to understand what a Smear::Detector is.
 */
Smear::Detector BuildDetector() {
  // Recall that the parameterisations passed to Smear::Device objects
  // give sigma(X) *not* sigma(X) / X
  // Electromagnetic calorimeter in the electron (east) direction, -4 < eta < -2
  // From Zhangbu's talk: sigma(E) / E = 2% / sqrt(E) \oplus 0.75%
  Smear::Device eastEcal("E", "0.02 * sqrt(E) + 0.0075 * E",
                         Smear::kElectromagnetic);
  eastEcal.Accept.AddZone(zoneEta(-4., -2.));
  // Tracking in the electron (east) direction, -2 < eta < -1
  // From Zhangbu's talk: sigma(p)/p = 1/(pT/pZ-1/6) x (0.45%pT \oplus 0.3%)
  //                                   \oplus (pZ/pT)x(0.2%/p/beta)
  // beta = p/E so 0.2%/p/beta = 0.2%E/p^2
  Smear::Device eastTracking("P",
    "P/(pT/pZ-1/6)*(0.0045*pT+0.003)+pZ/pT*0.002*E/P");
  eastTracking.Accept.AddZone(zoneEta(-2., -1.));
  eastTracking.Accept.SetCharge(Smear::kCharged);
  // Electromagnetic calorimeter in the barrel region, -1 < eta < 1
  // From Zhangbu's talk: sigma(E)/E = 14%/sqrt(E) \oplus 2%
  Smear::Device barrelEcal("E", "0.14*sqrt(E)+0.02*E", Smear::kElectromagnetic);
  barrelEcal.Accept.AddZone(zoneEta(-1., 1.));
  // Tracking in the barrel region, -1 < eta < 1
  // From Zhangbu's talk: sigma(p)/p = 0.45%pT \oplus 0.3% \oplus 0.2%/p/beta
  Smear::Device barrelTracking("P", "0.0045*pT*P+0.003*P+0.002*E/P");
  barrelTracking.Accept.AddZone(zoneEta(-1., 1.));
  barrelTracking.Accept.SetCharge(Smear::kCharged);
  // Tracking in hadron (west) direction, 1 < eta < 1.7
  // From Zhangbu's talk: as for east tracking, except 1/4 instead of 1/6
  Smear::Device westTracking("P",
    "P/(pT/pZ-1/4)*(0.0045*pT+0.003)+pZ/pT*0.002*E/P");
  westTracking.Accept.SetCharge(Smear::kCharged);
  westTracking.Accept.AddZone(zoneEta(1., 1.7));
  // Electromagnetic calorimeter in the hadron (west) direction, 1 < eta < 2
  // From Zhangbu's talk: sigma(E)/E = 16%/sqrt(E) \oplus 2%
  Smear::Device westEcal("E", "0.16*sqrt(E)+0.02*E", Smear::kElectromagnetic);
  westEcal.Accept.AddZone(zoneEta(1., 2.));
  // Electromagnetic calorimeter in the far-hadron (west) direction
  // 2.5 < eta < 5
  // From Zhangbu's talk: sigma(E)/E = 12%/sqrt(E) \oplus 1.4%
  Smear::Device westEcal2("E", "0.12*sqrt(E)+0.014*E", Smear::kElectromagnetic);
  westEcal2.Accept.AddZone(zoneEta(2.5, 5.));
  // Hadronic calorimeter in the far-hadron (west) direction, 2.5 < eta < 5
  // From Zhangbu's talk: sigma(E)/E = 38%/sqrt(E) \oplus 3%
  Smear::Device westHcal("E", "0.38*sqrt(E)+0.03*E", Smear::kHadronic);
  westHcal.Accept.AddZone(zoneEta(2.5, 5));
  // Assume perfect theta and phi performance i.e. momentum resolution
  // dominates over angular resolution
  Smear::Device theta("theta", "0");
  Smear::Device phi("phi", "0");
  // PID performance is unparameterised as of now
  Smear::PerfectID pid;
  // Combine the devices into a detector.
	Smear::Detector estar;
  estar.AddDevice(eastEcal);
  estar.AddDevice(barrelEcal);
  estar.AddDevice(westEcal);
  estar.AddDevice(westEcal2);
  estar.AddDevice(eastTracking);
  estar.AddDevice(barrelTracking);
  estar.AddDevice(westTracking);
  estar.AddDevice(westHcal);
  estar.AddDevice(theta);
  estar.AddDevice(phi);
  estar.AddDevice(pid);
  estar.SetEventKinematicsCalculator("NM JB DA");
  return estar;
}