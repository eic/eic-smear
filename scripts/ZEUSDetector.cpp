/**
 \file ZEUSDetector.cpp
 Example smearing script for the ZEUS detector

 \author    Thomas Burton
 \date      2014-01-10
 \copyright 2014 Brookhaven National Lab
 */

/**
 Convert pseudorapidity (eta) to polar angle (theta) in radians.
 Make use of TLorentzVector to do eta-to-theta conversion.
 */
double etaToTheta(const double eta) {
  TLorentzVector v;
  v.SetPtEtaPhiM(1., eta, 0., 0.);
  return v.Theta();
}

/**
 Convert and angle in degrees to one in radians.
 */
double degreesToRadians(double degrees) {
  return degrees / 180. * TMath::Pi();
}

/**
 Smearing parameterisations for the ZEUS detector.

 See JHEP05 (2009) 108.

 Note: you must gSystem->Load("libeicsmear") BEFORE loading this script.
 */
Smear::Detector BuildDetector() {
  // The central calorimeter, providing both electromagnetic and hadronic
  // calorimetry, but with different resolutions.
  // Note that this excludes the forward calorimeter.
	Smear::Device emCal(Smear::kE, "0.18*sqrt(E)", Smear::kElectromagnetic);
	Smear::Device hCal(Smear::kE, "0.35*sqrt(E)", Smear::kHadronic);
  // Calorimeter acceptance is +/- 4 in pseudorapidity.
  // Note that 4 is before -4 as eta-max corresponds to theta-min.
  Smear::Acceptance::Zone cal(etaToTheta(4.), etaToTheta(-4.));
  emCal.Accept.AddZone(cal);
  hCal.Accept.AddZone(cal);
  // ZEUS tracking acceptance covers 15 to 164 degrees in polar angle.
  // This is approximately +/- 2 in pseudorapidity.
  // Define with two volumes, for momentum and theta.
	Smear::Device theta(Smear::kTheta, "0.0005*P + 0.003");
	Smear::Device momentum(Smear::kP, "0.0085*P + 0.0025*P*P");
  Smear::Acceptance::Zone tracking(degreesToRadians(15.),
                                   degreesToRadians(164.));
  theta.Accept.AddZone(tracking);
  momentum.Accept.AddZone(tracking);
  // The acceptance of the Roman pots is a complicated function of
  // angle and momentum, and not easy to define using the available
  // Acceptance implementation. For completeness the parameterisations
  // are included below in case they are of use, and users wish to
  // define their own acceptance some way.
	// Smear::Device romanTran(Smear::kPt, "0.005");
	// Smear::Device romanLong(Smear::kPz, "0.005 * pZ");
  // There is no parameterisation for phi, so add a dummy device:
	Smear::Device phi(Smear::kPhi);
  // Add all the devices to the detector.
	Smear::Detector zeus;
  zeus.AddDevice(emCal);
  zeus.AddDevice(hCal);
  zeus.AddDevice(theta);
  zeus.AddDevice(momentum);
  zeus.AddDevice(phi);
  // zeus.AddDevice(romanTran);
  // zeus.AddDevice(romanLong);
  zeus.SetEventKinematicsCalculator("NM JB DA");
	return zeus;
}
