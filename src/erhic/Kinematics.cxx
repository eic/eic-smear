//
// Kinematics.cxx
//
// Created by TB on 7/7/11.
// Copyright 2011 BNL. All rights reserved.
//

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <numeric>

#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <TVector3.h>

#include "eicsmear/erhic/BeamParticles.h"
#include "eicsmear/erhic/EventDis.h"
#include "eicsmear/erhic/Kinematics.h"
#include "eicsmear/erhic/ParticleIdentifier.h"

namespace {
   typedef erhic::VirtualParticle Particle;
   //	==========================================================================
   // Helper function returning W^2 from x, Q^2 and mass.
   //	==========================================================================
   double computeW2FromXQ2M(double x, double Q2, double m) {
      return std::pow(m, 2.) + (1. - x) * Q2 / x;
   }
   //	==========================================================================
   // Returns the energy of a particle based on either its stored energy or
   // its total momentum and its mass.
   // If the particle's momentum is greater than zero, tracking information is
   // assumed to have been present and the particle's momentum is used to
   // calculate energy, via E^2 = p^2 + m^2.
   // For mass, if the particle ID is available, use the mass of that particle.
   // If not, assume the charged pion mass.
   // If momentum is not greater than zero, returns the stored energy of the
   // particle (only calorimeter information is assumed to be available).
   //	==========================================================================
   struct EnergyFromMomentumAndId :
   public std::unary_function<const Particle*, double> {
      double operator()(const Particle* p) const {
         if(not p) {
            return 0.;
         } // if
         // Use the stored energy as the default.
         double energy = p->GetE();
         // If momentum greater than zero, we assume the particle represents
         // data where tracking information was available, and we use the
         // momentum to compute energy.
         if(p->GetP() > 0.) {
            double mass(0.);
            TParticlePDG* pdg = p->Id().Info();
            if(pdg) {
               mass = pdg->Mass();
            } // if
            // Compute energy from p and m
            energy = sqrt(pow(p->GetP(), 2.) + pow(mass, 2.));
         } // if
         return energy;
      }
   };
   //	==========================================================================
   // Helper functor for calculating E - p_z of a particle.
   //	==========================================================================
   struct EMinusPz : public std::unary_function<const Particle*, double> {
      double operator()(const Particle* p) const {
         TLorentzVector fourMomentum(0., 0., 0., 0.);
         if(p) {
            fourMomentum = p->Get4Vector();
         } // if
         return fourMomentum.E() - fourMomentum.Pz();
      }
   };
   //	==========================================================================
   // Returns the energy-momentum 4-vector of a particle.
   // The returned energy is that computed using EnergyFromMomentumAndId, not
   // the stored energy of the particle.
   //	==========================================================================
   struct EnergyMomentum4Vector :
   public std::unary_function<const Particle*, TLorentzVector> {
      TLorentzVector operator()(const Particle* p) const {
         TLorentzVector ep;
         if(p) {
            ep = p->Get4Vector();
            // Attempt to calculate a (hopefully more precise) energy
            // using momentum and mass.
            ep.SetE(EnergyFromMomentumAndId()(p));
         } // if
         return ep;
      }
   };
} // anonymous namespace

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
   //	==========================================================================
   //	==========================================================================
   //	==========================================================================
   //	==========================================================================
   LeptonKinematicsComputer::LeptonKinematicsComputer(const EventDis& event) {
      ParticleIdentifier::IdentifyBeams(event, mBeams);
   }
   LeptonKinematicsComputer::LeptonKinematicsComputer(const BeamParticles& b)
   : mBeams(b) {
   }
   DisKinematics* LeptonKinematicsComputer::Calculate() {
      DisKinematics* kin = new DisKinematics(-1., -1., -1., -1., -1.);
      const TLorentzVector& s = mBeams.ScatteredLepton();
      // If there is no measurement of theta of the scattered lepton we
      // cannot calculate kinematics. If we have theta, we can calculate
      // using either energy or momentum (effectively energy, but generally
      // has better resolution), so use momentum if possible, energy if not.
      if(s.Theta() > 0. and (s.E() > 0 or s.P() > 0.)) {
         const TLorentzVector& l = mBeams.BeamLepton();
         const TLorentzVector& h = mBeams.BeamHadron();
         double cme = (l + h).M2();
         // Find scattered lepton energy boosted to the rest
         // frame of the hadron beam.
         double ELeptonInNucl = h.Gamma() * (l.P() - h.Beta() * l.Pz());
         double ELeptonOutNucl = h.Gamma() * (s.P() - h.Beta() * s.Pz());
         double theta = s.Theta();
         // Use momentum if available, energy if not.
         double energy = (s.P() > 0. ? s.P() : s.E());
         kin->mQ2 = 2. * l.P() * energy * (1. + cos(theta));
         kin->mNu = ELeptonInNucl - ELeptonOutNucl;
         kin->mY = kin->mNu * 2. * h.M() / cme;
         kin->mX = kin->mQ2 / kin->mY / cme;
         kin->mW2 = computeW2FromXQ2M(kin->mX, kin->mQ2, h.M());
      } // if
      return kin;
   }
   //	==========================================================================
   //	==========================================================================
   //	==========================================================================
   //	==========================================================================
   JacquetBlondelComputer::JacquetBlondelComputer(const EventDis& event)
   : mEvent(event) {
   }
   DisKinematics* JacquetBlondelComputer::Calculate() {
      // Get all the final-state particles except the scattered lepton.
      mParticles.clear();
      mEvent.HadronicFinalState(mParticles);
      DisKinematics* kin = new DisKinematics;
      kin->mY = ComputeY();
      kin->mQ2 = ComputeQSquared();
      kin->mX = ComputeX();
      kin->mW2 = computeW2FromXQ2M(kin->mX, kin->mQ2,
                                   mEvent.BeamHadron()->GetM());
      return kin;
   }
   //	==========================================================================
   //	==========================================================================
   Double_t JacquetBlondelComputer::ComputeY() const {
      double y(0.);
      // Calculate y, as long as we have beam information.
      const VirtualParticle* hadron = mEvent.BeamHadron();
      const VirtualParticle* lepton = mEvent.BeamLepton();
      if(hadron and lepton) {
         // Sum the energies of the final-state hadrons
         std::list<double> E;
         std::transform(mParticles.begin(),
                        mParticles.end(),
                        std::back_inserter(E),
                        EnergyFromMomentumAndId());
         const double sumEh = std::accumulate(E.begin(), E.end(), 0.);
         // Sum the pz of the final-state hadrons
         std::list<double> pz;
         std::transform(mParticles.begin(),
                        mParticles.end(),
                        std::back_inserter(pz),
                        std::mem_fun(&VirtualParticle::GetPz));
         const double sumPzh = std::accumulate(pz.begin(), pz.end(), 0.);
         // Compute y.
         // This expression seems more accurate at small y than the usual
         // (sumE - sumPz) / 2E_lepton.
         y = (hadron->GetE() * sumEh -
              hadron->GetPz() * sumPzh -
              pow(hadron->GetM(), 2.)) /
         (lepton->GetE() * hadron->GetE() - lepton->GetPz() * hadron->GetPz());
      } // if
      // Return y, bounding it by the range [0, 1]
      return std::max(0., std::min(y, 1.));
   }
   //	==========================================================================
   // We use the following exact expression:
   // 2(E_p * sum(E_h) - pz_p * sum(pz_h)) - sum(E_h)^2 + sum(p_h)^2 - M_p^2
   //	==========================================================================
   Double_t JacquetBlondelComputer::ComputeQSquared() const {
      double Q2(0.);
      // Calculate Q^2, as long as we have beam information.
      const VirtualParticle* hadron = mEvent.BeamHadron();
      if(hadron) {
         // Get the px of each particle:
         std::list<double> px;
         std::transform(mParticles.begin(),
                        mParticles.end(),
                        std::back_inserter(px),
                        std::mem_fun(&VirtualParticle::GetPx));
         // Get the py of each particle:
         std::list<double> py;
         std::transform(mParticles.begin(),
                        mParticles.end(),
                        std::back_inserter(py),
                        std::mem_fun(&VirtualParticle::GetPy));
         double sumPx = std::accumulate(px.begin(), px.end(), 0.);
         double sumPy = std::accumulate(py.begin(), py.end(), 0.);
         Q2 = (pow(sumPx, 2.) + pow(sumPy, 2.)) / (1. - ComputeY());
      } // if
      return Q2;
   }
   //	==========================================================================
   // x_JB = Q2_JB / (y_JB * s)
   //	==========================================================================
   Double_t JacquetBlondelComputer::ComputeX() const {
      double x(0.);
      // Calculate x, as long as we have beam information.
      const VirtualParticle* hadron = mEvent.BeamHadron();
      const VirtualParticle* lepton = mEvent.BeamLepton();
      if(hadron and lepton) {
         double y = ComputeY();
         double s = (hadron->Get4Vector() + lepton->Get4Vector()).M2();
         // Use the approximate relation Q^2 = sxy to calculate x.
         if(y > 0.0) {
            x = ComputeQSquared() / y / s;
         } // if
      } // if
      return x;
   }
   //	==========================================================================
   //	==========================================================================
   //	==========================================================================
   //	==========================================================================
   DoubleAngleComputer::DoubleAngleComputer(const EventDis& event)
   : mEvent(event) {
   }
   DisKinematics* DoubleAngleComputer::Calculate() {
      mParticles.clear();
      mEvent.HadronicFinalState(mParticles);
      mHasChanged = true;
      DisKinematics* kin = new DisKinematics;
      kin->mY = computeY();
      kin->mQ2 = computeQSquared();
      kin->mX = computeX();
      // Calculate W^2 from M^2 + (1 - x) * Q^2 / x
      kin->mW2 = computeW2FromXQ2M(kin->mX, kin->mQ2,
                                   mEvent.BeamHadron()->GetM());
      return kin;
   }
   //	==========================================================================
   // Scattering angle of struck quark
   // cos(angle) = A / B
   // where A = (sum of px_h)^2 + (sum of py_h)^2 - (sum of [E_h - pz_h])^2
   // and   B = (sum of px_h)^2 + (sum of py_h)^2 + (sum of [E_h - pz_h])^2
   // This is a computationally expensive operation that is called a lot,
   // so cashe the result until the particle list changes.
   //	==========================================================================
   Double_t DoubleAngleComputer::computeQuarkAngle() const {
      using namespace std;
      // Return the cached value if no changes have occurred since
      // the last call to computeQuarkAngle().
      if(not mHasChanged) {
         return mAngle;
      } // if
      list<double> px, py, eMinusPZ;
      // Get the px of each particle:
      transform(mParticles.begin(), mParticles.end(), back_inserter(px),
                mem_fun(&VirtualParticle::GetPx));
      // Get the py of each particle:
      transform(mParticles.begin(), mParticles.end(), back_inserter(py),
                mem_fun(&VirtualParticle::GetPy));
      // Get the (E-pz) of each particle:
      transform(mParticles.begin(), mParticles.end(), back_inserter(eMinusPZ),
                EMinusPz());
      // Square the sum of each momentum component
      struct SquareOfSum {
         double operator()(const list<double>& l) const {
            return pow(accumulate(l.begin(), l.end(), 0.0), 2.0);
         }
      } squareOfSum;
      double sumOfPx2 = squareOfSum(px);
      double sumOfPy2 = squareOfSum(py);
      double sumOfEMinusPz = squareOfSum(eMinusPZ);
      // Calculate the angle:
      double denominator = (sumOfPx2 + sumOfPy2 + sumOfEMinusPz);
      // Set cos(angle) to 1 by default. If denominator is zero or less,
      // angle will then evaluate to acos(1) = 0.
      double cosAngle(1.0);
      if(denominator > 0.0) {
         cosAngle = (sumOfPx2 + sumOfPy2 - sumOfEMinusPz)
         / denominator;
      } // if
      mAngle = TMath::ACos(cosAngle);
      mHasChanged = false;
      return mAngle;
   }
   //	==========================================================================
   // y_DA = Q2_DA / (x_DA * s)
   //	==========================================================================
   Double_t DoubleAngleComputer::computeY() const {
      double x = computeX();
      double s = 4. * mEvent.BeamLepton()->GetE() * mEvent.BeamHadron()->GetE();
      if(x > 0.0) {
         return computeQSquared() / x / s;
      } // if
      return 0.0;
   }
   //	==========================================================================
   // a_q = quark scattering angle, a_e = electron scattering angle
   // Q2_DA = 4 * E_e^2 * sin(a_q) (1 + cos(a_e)) /
   //         [ sin(a_q) + sin(a_e) - sin(a_q + a_e) ]
   //	==========================================================================
   Double_t DoubleAngleComputer::computeQSquared() const {
      double theta = mEvent.ScatteredLepton()->GetTheta();
      double numerator = 4.0
      * pow(mEvent.BeamLepton()->GetE(), 2.0)
      * sin(computeQuarkAngle())
      * (1.0 + cos(theta));
      double denominator = sin(computeQuarkAngle())
      + sin(theta)
      - sin(computeQuarkAngle() + theta);
      if(denominator > 0.0) {
         return numerator / denominator;
      } // if
      return 0.0;
   }
   //	==========================================================================
   // x_DA = (E_e/E_p)[sin(a_q)+sin(a_e)+sin(a_q+a_e)]/
   //        [sin(a_q)+sin(a_e)-sin(a_q+a_e)]
   //	==========================================================================
   Double_t DoubleAngleComputer::computeX() const {
      double theta = mEvent.ScatteredLepton()->GetTheta();
      double numerator = sin(computeQuarkAngle())
      + sin(theta)
      + sin(computeQuarkAngle() + theta);
      double denominator = sin(computeQuarkAngle())
      + sin(theta)
      - sin(computeQuarkAngle() + theta);
      if(denominator > 0.0) {
         return mEvent.BeamLepton()->GetE() / mEvent.BeamHadron()->GetE() *
         numerator / denominator;
      } // if
      return 0.0;
   }
} // namespace erhic
