/**
 \file STARDetector.cpp
 Example smearing script for the STAR detector

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
 Returns an Acceptance::Zone spanning a range in eta.
 */
Smear::Acceptance::Zone makeZone(double etaMin, double etaMax) {
  // Note that we need to flip the order of arguments because the numerically
  // larger eta corresponds to the minimum angle and the smaller eta to the
  // maximum angle.
  // e.g. for eta -1 to 1, -1 --> theta=2.4 and 1 --> theta=0.7
  return Smear::Acceptance::Zone(etaToTheta(etaMax), etaToTheta(etaMin));
}

/**
 Smearing parameterisations for the STAR detector.
 
 These parameterisations are non-exhaustive: they do not cover elements such
 as particle identification, and they are only for the central elements of
 the detector - essentially, just TPC and B/EEMC.

 Note: you must gSystem->Load("libeicsmear") BEFORE loading this script,
 as ROOT needs to understand what a Smear::Detector is.
 */
Smear::Detector BuildDetector() {
  // Electromagnetic calorimeter.
  // Acceptance is valid for the BEMC (-1 < eta < 1) and EEMC (1 < eta < 2)
  // so define a single acceptance zone -1 < eta < 2.
  // Acceptance is defined in terms of angle (theta) not pseudorapidity (eta)
  // so we need to convert eta -> theta.
  Smear::Device emCal(Smear::kE,
                      "0.015*E+0.14*sqrt(E)",
                      Smear::kElectromagnetic);
  emCal.Accept.AddZone(makeZone(-1., 2.));
  // Tracking, comprising momentum and theta.
  // These just span the TPC, -1 < eta < 1.
  Smear::Device momentum(Smear::kP, "0.005*P+0.004*P*P");
  momentum.Accept.SetCharge(Smear::kCharged);
  momentum.Accept.AddZone(makeZone(-1., 1.));
  // See comments at the end of the file for Zhangbu's email about
  // theta resolution.
  Smear::Device theta(Smear::kTheta,
                      "sqrt(pow(3e-4, 2) + pow(9e-4 / P, 2) / sin(theta))");
  theta.Accept.SetCharge(Smear::kCharged);
  theta.Accept.AddZone(makeZone(-1., 1.));
  // We don't have a parameterisation for phi, so just use a
  // device that gives perfect resolution.
	Smear::Device phi(Smear::kPhi, "0");
	phi.Accept.SetCharge(Smear::kCharged);
  // Combine the devices into a detector.
	Smear::Detector star;
	star.SetLegacyMode ( true ); // turn off checks and enhanced momentum consistency
  star.AddDevice(emCal);
  star.AddDevice(momentum);
  star.AddDevice(theta);
  star.AddDevice(phi);
  star.SetEventKinematicsCalculator("NM JB DA");
  return star;
}

/*
Here is the email from Zhangbu giving the parameterisations
used in the above.


Date: Thu, 2 Jun 2011 17:35:45 -0400
From: "Xu, Zhangbu" <xzb@bnl.gov>
To: "Aschenauer, Elke" <elke@bnl.gov>
Subject: theta angular resolution

Hi, Elke:

Talk to a couple of people doing the TPC calibration,
from drift velocity and TPC cluster uncertainties, the
d (theta) ~=3xe-4

I am working on the MC simulation to get the resolution
due to multiple scattering, this will be momentum dependent.
A hand calculation of the theta angular spread due to multiple
scattering for STAR TPC material (mainly beam pipe and inner field
cage ~=0.5%X0) is
theta0 = 13.6MeV/beta/p*sqrt(x/X0)~=0.9MeV/p,
So for a 2GeV/c electron, the theta angular resolution will be
5e-4, close to the detector resolution.

So it is probably reasonable to assume the theta angular resolution
to be: d(theta) = sqrt((3xe-4)**2+(0.9MeV/p)**2/sin(theta))


Thanks!

Zhangbu
*/
