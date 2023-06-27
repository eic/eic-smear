/**
 \file
 Implementations of kinematic calculation classes.
 
 \author    Thomas Burton
 \date      2011-07-07
 \copyright 2011 Brookhaven National Lab
 */

#include "eicsmear/erhic/Kinematics.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <list>
#include <memory>
#include <numeric>  // For std::accumulate
#include <stdexcept>
#include <utility>
#include <vector>

#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <TVector3.h>

#include "eicsmear/erhic/EventDis.h"
#include "eicsmear/erhic/ParticleIdentifier.h"

bool erhic::DisKinematics::BoundaryWarning=true;
namespace {

const double chargedPionMass =
TDatabasePDG::Instance()->GetParticle(211)->Mass();

typedef erhic::VirtualParticle Particle;

// ==========================================================================
// Helper function returning W^2 from x, Q^2 and mass.
// ==========================================================================
double computeW2FromXQ2M(double x, double Q2, double m) {
  if (x > 0.) {
    return std::pow(m, 2.) + (1. - x) * Q2 / x;
  }  // if
  return 0.;
}

// ==========================================================================
// Returns the value x bounded by [minimum, maximum].
// ==========================================================================
double bounded(double x, double minimum, double maximum) {
  // KK 12/10/19: x>1 can be physically possible (in eA).
  // In general, silently cutting off values is dangerous practice.
  // Instead, we will optionally issue a warning but accept
  // Note: The warning is on by default, because this (local) function does
  // often get called with values for minimum and maximum,
  // and we shouldn't just silently ignore that.
  // It may be better to remove "bounded" entirely though.
  
  if ( erhic::DisKinematics::BoundaryWarning && ( x< minimum || x > maximum ) ){
    std::cerr << "Warning in Kinematics, bounded(): x (or y) = " << x
	      << " is outside [" << minimum << "," << maximum << "]" << std::endl;
    std::cerr << "To disable this warning, set erhic::DisKinematics::BoundaryWarning=false;" << std::endl;
  }
  return x; 
  // return std::max(minimum, std::min(x, maximum));
}

/*
 Calculates kinematic information that is knowable by a detector.
 
 This takes account whether p, E and particle id are known or not, and
 computes assumed values for variables that aren't explicitly known.
 As a simple example: if p is not known, but E is, assume p == E.
 Returns a ParticleMC with these "best-guess" values, specifically p, E,
 ID, mass, pz, px, py, pt.
 */
using erhic::ParticleMC;
using erhic::VirtualParticle;
typedef std::vector<const VirtualParticle*>::iterator VirtPartIter;
typedef std::vector<const VirtualParticle*>::const_iterator cVirtPartIter;

class MeasuredParticle {
 public:
  static ParticleMC* Create(const erhic::VirtualParticle* particle) {
    if (!particle) {
      throw std::invalid_argument("MeasuredParticle given NULL pointer");
    }  // if
    ParticleMC* measured = new ParticleMC;
    // Copy ID from the input particle or guess it if not known.
    measured->SetId(CalculateId(particle));
    // Set mass from known/guessed ID.
    TParticlePDG* pdg = measured->Id().Info();
    if (pdg) {
      measured->SetM(pdg->Mass());
    }  // if
    std::pair<double, double> ep =
      CalculateEnergyMomentum(particle, pdg->Mass());
    TLorentzVector vec(0., 0., ep.second, ep.first);
    vec.SetTheta(particle->GetTheta());
    vec.SetPhi(particle->GetPhi());
    measured->Set4Vector(vec);
    return measured;
  }
  /*
   Determine the particle ID.
   If the input particle ID is known, returns the same value.
   If the input ID is not known, return an assumed ID based on available
   information.
   */
  static int CalculateId(const erhic::VirtualParticle* particle) {
    int id(0);
    TParticlePDG* pdg = particle->Id().Info();
    // Skip pid of 0 (the default) as this is a dummy "ROOTino".
    if (pdg && particle->Id() != 0) {
      id = particle->Id();
    } else if (particle->GetP() > 0.) {
      // The particle ID is unknown.
      // If there is momentum information, it must be charged so
      // assume it is a pi+. We don't currently track whether charge info
      // is known so we can't distinguish pi+ and pi- here.
      // If there is no momentum information, assume it is a photon.
      // *This is not really correct, but currently we don't track
      // whether a particle is known to be EM or hadronic, so we can't
      // e.g. assume neutron for hadronic particles.
      id = 211;
    } else {
      id = 22;
    }  // if
    return id;
  }
  /*
   Calculate energy, using momentum and mass information if available.
   If mass is greater than zero, use this as the assumed mass when
   calculating via momentum, otherwise calculate
   mass from scratch via either the known or assumed ID returned by
   CalculateID(). Use this argument if you already know this mass and
   don't want to repeat the calculation.
   */
  static std::pair<double, double> CalculateEnergyMomentum(
      const erhic::VirtualParticle* particle, double mass = -1.) {
    if (mass < 0.) {
      int id = CalculateId(particle);
      TParticlePDG* pdg = TDatabasePDG::Instance()->GetParticle(id);
      if (pdg) {
        mass = pdg->Mass();
      } else {
        mass = 0.;
      }  // if
    }  // if
    // If momentum greater than zero, we assume the particle represents
    // data where tracking information was available, and we use the
    // momentum to compute energy.
    std::pair<double, double> ep(0., 0.);
    if (particle->GetP() > 0.) {
      ep.first = sqrt(pow(particle->GetP(), 2.) + pow(mass, 2.));
      ep.second = particle->GetP();
    } else if (particle->GetE() > 0.) {
      ep.first = particle->GetE();
      // need to catch cases where E was smeared < E
      // Assign P=0 in this case
      auto E = particle->GetE();
      if ( E >= mass ){
	ep.second = sqrt(pow(E, 2.) - pow(mass, 2.));
      } else {
	ep.second = 0;
      }
    }  // if
    return ep;
  }
};

}  // anonymous namespace

namespace erhic {

DisKinematics::DisKinematics()
: mX(0.)
, mQ2(0.)
, mW2(0.)
, mNu(0.)
, mY(0) {
}

DisKinematics::DisKinematics(double x, double y, double nu,
                             double Q2, double W2)
: mX(x)
, mQ2(Q2)
, mW2(W2)
, mNu(nu)
, mY(y) {
}

// ==========================================================================
// ==========================================================================
LeptonKinematicsComputer::LeptonKinematicsComputer(const EventDis& event) {
  //   ParticleIdentifier::IdentifyBeams(event, mBeams);
  mBeams.push_back(event.BeamLepton());
  mBeams.push_back(event.BeamHadron());
  mBeams.push_back(event.ExchangeBoson());
  mBeams.push_back(event.ScatteredLepton());
}

// ==========================================================================
// ==========================================================================
DisKinematics* LeptonKinematicsComputer::Calculate() {
  // Create kinematics with default values.
  DisKinematics* kin = new DisKinematics;
  try {
    // Use E, p and ID of scattered lepton to create "best-guess" kinematics.
    // MeasuredParticle::Create will throw an exception in case of a NULL
    // pointer argument.
    std::unique_ptr<const VirtualParticle> scattered( MeasuredParticle::Create(mBeams.at(3)));
    // If there is no measurement of theta of the scattered lepton we
    // cannot calculate kinematics. Note that via MeasuredParticle the momentum
    // may actually be derived from measured energy instead of momentum.
    if (scattered->GetTheta() > 0. && scattered->GetP() > 0.) {
      const TLorentzVector& l = mBeams.at(0)->Get4Vector();
      const TLorentzVector& h = mBeams.at(1)->Get4Vector();
      // Calculate kinematic quantities, making sure to bound
      // the results by physical limits.
      // First, Q-squared.
      // double Q2 = 2. * l.E() * scattered->GetE() * (1. + cos(scattered->GetTheta()));
      double Q2 = 2. * ( l.E() * scattered->GetE() +  scattered->GetP()*l.P()*cos(scattered->GetTheta()) - l.M2()*l.M2() );
      kin->mQ2 = std::max(0., Q2);
      // Find scattered lepton energy boosted to the rest
      // frame of the hadron beam. Calculate nu from this.
      double gamma = mBeams.at(1)->GetE() / mBeams.at(1)->GetM();
      double beta = mBeams.at(1)->GetP() / mBeams.at(1)->GetE();
      double ELeptonInNucl = gamma * (l.P() - beta * l.Pz());
      double ELeptonOutNucl = gamma *
      (scattered->GetP() - beta * scattered->GetPz());
      kin->mNu = std::max(0., ELeptonInNucl - ELeptonOutNucl);
      // Calculate y using the exchange boson.
      // Determine the exchange boson 4-vector from the scattered lepton, as
      // this will then be valid for smeared event input also (where the
      // exchange boson is not recorded).
      const TLorentzVector q = l - scattered->Get4Vector();
      const double ldoth = l.Dot(h);
      double y = q.Dot(h) / ldoth ;
      if ( y<0) y=0; // kk: catching unphysical negative values.
      kin->mY = bounded(y, 0., 1.);
      double x = 0;
      // Calculate x from Q^2 = sxy
      // double cme = (l + h).M2();
      // if ( y>0 && cme>0 ) x = kin->mQ2 / kin->mY / cme;
      // Calculate from x = Q^2 / ( 2 * h*q ), y = h*q / (l * h)
      if ( y>0 && std::abs(ldoth) > 0) x = kin->mQ2 / 2 / y / ldoth;	    
      if ( x<0 ) x=0; // kk: catching unphysically negative values.
      kin->mX = bounded(x, 0., 1.);
      kin->mW2 = computeW2FromXQ2M(kin->mX, kin->mQ2, h.M());
      // if ( y > 1 ) std::cout << " NM y = " << y << std::endl;
      // if ( x > 1 ) std::cout << " NM x = " << x << "    y = " << y << "    Q2 = " << kin->mQ2 << std::endl;
    }  // if
  }  // try
  catch(std::invalid_argument& except) {
    // In case of no scattered lepton return default values.
    std::cerr << "No lepton for kinematic calculations" << std::endl;
  }  // catch
  return kin;
}

// ==========================================================================
// ==========================================================================
JacquetBlondelComputer::~JacquetBlondelComputer() {
  // Delete all "measureable" particles.
  for (VirtPartIter i = mParticles.begin(); i != mParticles.end(); ++i) {
    if (*i) {
      delete *i;
      *i = NULL;
    }  // if
  }  // for
  mParticles.clear();
}

// ==========================================================================
// ==========================================================================
JacquetBlondelComputer::JacquetBlondelComputer(const EventDis& event)
: mEvent(event) {
  // Get the full list of final-state particles in the event.
  // including decay bosons and leptons 
  std::vector<const erhic::VirtualParticle*> final;
  mEvent.HadronicFinalState(final);
  // Populate the stored particle list with "measurable" versions of
  // each final-state particle.
  // std::transform(final.begin(), final.end(), std::back_inserter(mParticles),
  //                std::ptr_fun(&MeasuredParticle::Create));
  // transform applies the function MeasuredParticle::Create to each element of
  // final and stores the result in mParticles.
  // Easier, and compatible with C++98 and C++17:
  for (VirtPartIter it = final.begin() ; it != final.end(); ++it){
    mParticles.push_back ( MeasuredParticle::Create ( *it ) );
  }
}

// ==========================================================================
// ==========================================================================
DisKinematics* JacquetBlondelComputer::Calculate() {
  // Get all the final-state particles except the scattered lepton.
  DisKinematics* kin = new DisKinematics;
  kin->mY = ComputeY();
  kin->mQ2 = ComputeQSquared();
  kin->mX = ComputeX();
  kin->mW2 = computeW2FromXQ2M(kin->mX, kin->mQ2,
                               mEvent.BeamHadron()->GetM());
  return kin;
}

// ==========================================================================
// ==========================================================================
Double_t JacquetBlondelComputer::ComputeY() const {
  double y(0.);
  // Calculate y, as long as we have beam information.
  const VirtualParticle* hadron = mEvent.BeamHadron();
  const VirtualParticle* lepton = mEvent.BeamLepton();
  if (hadron && lepton) {
    // Sum the energies of the final-state hadrons
    std::vector<double> E;
    for (cVirtPartIter it = mParticles.begin() ; it != mParticles.end(); ++it){
      E.push_back ( (*it)->GetE() );
    }
    const double sumEh = std::accumulate(E.begin(), E.end(), 0.);

    // Sum the pz of the final-state hadrons
    std::vector<double> pz;
    for (cVirtPartIter it = mParticles.begin() ; it != mParticles.end(); ++it){
      pz.push_back ( (*it)->GetPz() );
    }
    const double sumPzh = std::accumulate(pz.begin(), pz.end(), 0.);
    // Compute y.
    // This expression seems more accurate at small y than the usual
    // (sumE - sumPz) / 2E_lepton.
    // Tested by Xiaoxuan Chu and Elke Aschenauer, this formula
    // is more appropriate in the presence of radiative effects.
    y = (hadron->GetE() * sumEh -
         hadron->GetPz() * sumPzh -
         pow(hadron->GetM(), 2.)) /
    (lepton->GetE() * hadron->GetE() - lepton->GetPz() * hadron->GetPz());

    // Can result in negative y. Catch that here since bounded() below is
    // for legacy reasons only and produces warnings instead of cutoffs
    if (y<0) y=0;
  }  // if

  // Return y, bounding it by the range [0, 1]
  // return bounded(y, 0., 1.);
  // changed: y_jb can be slightly >1
  // we're returning y as is and circumvent the warnning
  //std::cout << " JB y = " << y << std::endl;
  return y;
}

// ==========================================================================
// We use the following exact expression:
// 2(E_p * sum(E_h) - pz_p * sum(pz_h)) - sum(E_h)^2 + sum(p_h)^2 - M_p^2
// ==========================================================================
Double_t JacquetBlondelComputer::ComputeQSquared() const {
  double Q2(0.);
  // Calculate Q^2, as long as we have beam information.
  const VirtualParticle* hadron = mEvent.BeamHadron();
  if (hadron) {
    // Get the px of each particle:
    std::vector<double> px;
    for (cVirtPartIter it = mParticles.begin() ; it != mParticles.end(); ++it){
      px.push_back ( (*it)->GetPx() );
    }

    // Get the py of each particle:
    std::vector<double> py;
    for (cVirtPartIter it = mParticles.begin() ; it != mParticles.end(); ++it){
      py.push_back ( (*it)->GetPy() );
    }

    double sumPx = std::accumulate(px.begin(), px.end(), 0.);
    double sumPy = std::accumulate(py.begin(), py.end(), 0.);
    double y = ComputeY();
    if (y < 1.) {
      Q2 = (pow(sumPx, 2.) + pow(sumPy, 2.)) / (1. - y);
    }  // if
  }  // if
  return std::max(0., Q2);
}

// ==========================================================================
// x_JB = Q2_JB / (y_JB * s)
// ==========================================================================
Double_t JacquetBlondelComputer::ComputeX() const {
  double x(0.);
  // Calculate x, as long as we have beam information.
  const VirtualParticle* hadron = mEvent.BeamHadron();
  const VirtualParticle* lepton = mEvent.BeamLepton();
  if (hadron && lepton) {
    double y = ComputeY();
    double s = (hadron->Get4Vector() + lepton->Get4Vector()).M2();
    // Use the approximate relation Q^2 = sxy to calculate x.
    if (y > 0.0) {
      x = ComputeQSquared() / y / s;
    }  // if
  }  // if

  // return bounded(x, 0., 1.);
  // changed: x_jb can be slightly >1 when y_jb is very small
  // we're returning x as is and circumvent the warning (but keep warning for x<0)
  // std::cout << " JB x = " << x << std::endl;
  return bounded(x, 0., 100.);
}

// ==========================================================================
// ==========================================================================
DoubleAngleComputer::~DoubleAngleComputer() {
  // Delete all "measureable" particles.
  typedef std::vector<const VirtualParticle*>::iterator Iter;
  for (Iter i = mParticles.begin(); i != mParticles.end(); ++i) {
    if (*i) {
      delete *i;
      *i = NULL;
    }  // if
  }  // for
  mParticles.clear();
}

// ==========================================================================
// ==========================================================================
DoubleAngleComputer::DoubleAngleComputer(const EventDis& event)
: mEvent(event) {
  // Get the full list of final-state particles in the event.
  // including decay bosons and leptons 
  std::vector<const erhic::VirtualParticle*> final;
  mEvent.HadronicFinalState(final);
  // Populate the stored particle list with "measurable" versions of
  // each final-state particle.
  for (VirtPartIter it = final.begin() ; it != final.end(); ++it){
    mParticles.push_back ( MeasuredParticle::Create ( *it ) );
  }
}

// ==========================================================================
// Formulae used are from F.D. Aaron et al., JHEP01 (2010) 109.
// ==========================================================================
DisKinematics* DoubleAngleComputer::Calculate() {
  mHasChanged = true;
  DisKinematics* kin = new DisKinematics;
  kin->mQ2 = ComputeQSquared();
  kin->mX = ComputeX();
  kin->mY = ComputeY();
  // Calculate W^2 from M^2 + (1 - x) * Q^2 / x
  kin->mW2 = computeW2FromXQ2M(kin->mX, kin->mQ2, mEvent.BeamHadron()->GetM());
  return kin;
}

// ==========================================================================
// Scattering angle of struck quark
// cos(angle) = A / B
// where A = (sum of px_h)^2 + (sum of py_h)^2 - (sum of [E_h - pz_h])^2
// and   B = (sum of px_h)^2 + (sum of py_h)^2 + (sum of [E_h - pz_h])^2
// This is a computationally expensive operation that is called a lot,
// so cashe the result until the particle list changes.
// ==========================================================================
Double_t DoubleAngleComputer::ComputeQuarkAngle() const {
  // Return the cached value if no changes have occurred since
  // the last call to computeQuarkAngle().
  if (!mHasChanged) {
    return mAngle;
  }  // if
  std::vector<TLorentzVector> hadrons;
  // std::transform(mParticles.begin(),
  //                mParticles.end(),
  //                std::back_inserter(hadrons),
  //                std::mem_fun(&VirtualParticle::Get4Vector));
  for (cVirtPartIter it = mParticles.begin() ; it != mParticles.end(); ++it){
    hadrons.push_back ( (*it)->Get4Vector() );
  }
  TLorentzVector h = std::accumulate(hadrons.begin(),
                                     hadrons.end(),
                                     TLorentzVector(0., 0., 0., 0.));
  mAngle = 2. * TMath::ATan((h.E() - h.Pz()) / h.Pt());
  mHasChanged = false;
  return mAngle;
}

// ==========================================================================
// ==========================================================================
Double_t DoubleAngleComputer::ComputeY() const {
  double y(0.);
  const VirtualParticle* scattered = mEvent.ScatteredLepton();
  if (scattered) {
    double theta = scattered->GetTheta();
    double gamma = ComputeQuarkAngle();
    double denominator = tan(theta / 2.) + tan(gamma / 2.);
    if (denominator > 0.) {
      y = tan(gamma / 2.) / denominator;
    }  // if
  }  // if
     // Return y bounded by the range [0, 1].
  return bounded(y, 0., 1.);
}

// ==========================================================================
// ==========================================================================
Double_t DoubleAngleComputer::ComputeQSquared() const {
  double Q2(0.);
  const VirtualParticle* lepton = mEvent.BeamLepton();
  const VirtualParticle* scattered = mEvent.ScatteredLepton();
  if (lepton && scattered) {
    double theta = scattered->GetTheta();
    double gamma = ComputeQuarkAngle();
    double denominator = tan(theta / 2.) + tan(gamma / 2.);
    // kk: could catch theta=0, but it's a canary for a faulty detector   if (denominator > 0. && tan(theta / 2.) !=0 ) {
    if (denominator > 0. ) {
      Q2 = 4. * pow(lepton->GetE(), 2.) / tan(theta / 2.) / denominator;
    }  // if
  }  // if
  // Return Q^2, requiring it to be positive.
  return std::max(0., Q2);
}

// ==========================================================================
// ==========================================================================
Double_t DoubleAngleComputer::ComputeX() const {
  double x(0.);
  const VirtualParticle* lepton = mEvent.BeamLepton();
  const VirtualParticle* hadron = mEvent.BeamHadron();
  if (lepton && hadron) {
    double s = (lepton->Get4Vector() + hadron->Get4Vector()).M2();
    double y = ComputeY();
    if (s > 0. && y > 0.) {
      x = ComputeQSquared() / y / s;
    }  // if
  }  // if
  // Return x bounded by the range [0, 1].
  return bounded(x, 0., 1.);
}

}  // namespace erhic
