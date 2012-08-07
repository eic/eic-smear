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
#include <numeric> // For std::accumulate

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
      if(x > 0.) {
         return std::pow(m, 2.) + (1. - x) * Q2 / x;
      } // if
      return 0.;
   }
   //	==========================================================================
   // Returns the value x bounded by [minimum, maximum].
   //	==========================================================================
   double bounded(double x, double minimum, double maximum) {
      return std::max(minimum, std::min(x, maximum));
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
   //	==========================================================================
   // Helper functor for calculating E - p_z of a particle.
   //	==========================================================================
   struct EMinusPz : public std::unary_function<const Particle*, double> {
      double operator()(const Particle* p) const {
         TLorentzVector fourMomentum(0., 0., 0., 0.);
         if(p) {
            fourMomentum = EnergyMomentum4Vector()(p);
         } // if
         return fourMomentum.E() - fourMomentum.Pz();
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
      DisKinematics* kin = new DisKinematics(0., 0., 0., 0., 0.);
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
         // Calculate kinematic quantities, making sure to bound
         // the results by physical limits.
         double Q2 = 2. * l.P() * energy * (1. + cos(theta));
         kin->mQ2 = std::max(0., Q2);
         kin->mNu = std::max(0., ELeptonInNucl - ELeptonOutNucl);
         double y = kin->mNu * 2. * h.M() / cme;
         kin->mY = bounded(y, 0., 1.);
         double x = kin->mQ2 / kin->mY / cme;
         kin->mX = bounded(x, 0., 1.);
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
      return bounded(y, 0., 1.);
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
         double y = ComputeY();
         if(y < 1.) {
             Q2 = (pow(sumPx, 2.) + pow(sumPy, 2.)) / (1. - y);
         } // if
      } // if
      return std::max(0., Q2);
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
      return bounded(x, 0., 1.);
   }
   //	==========================================================================
   //	==========================================================================
   //	==========================================================================
   //	==========================================================================
   DoubleAngleComputer::DoubleAngleComputer(const EventDis& event)
   : mEvent(event) {
   }
   //	==========================================================================
   // Formulae used are from F.D. Aaron et al., JHEP01 (2010) 109.
   //	==========================================================================
   DisKinematics* DoubleAngleComputer::Calculate() {
      mParticles.clear();
      mEvent.HadronicFinalState(mParticles);
      mHasChanged = true;
      DisKinematics* kin = new DisKinematics;
      kin->mQ2 = ComputeQSquared();
      kin->mX = ComputeX();
      kin->mY = ComputeY();
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
   Double_t DoubleAngleComputer::ComputeQuarkAngle() const {
      // Return the cached value if no changes have occurred since
      // the last call to computeQuarkAngle().
      if(not mHasChanged) {
         return mAngle;
      } // if
      std::list<TLorentzVector> hadrons;
      std::transform(mParticles.begin(),
                     mParticles.end(),
                     std::back_inserter(hadrons),
                     EnergyMomentum4Vector());
      TLorentzVector h = std::accumulate(hadrons.begin(),
                                         hadrons.end(),
                                         TLorentzVector(0., 0., 0., 0.));
      mAngle = 2. * TMath::ATan((h.E() - h.Pz()) / h.Pt());
      mHasChanged = false;
      return mAngle;
   }
   //	==========================================================================
   //	==========================================================================
   Double_t DoubleAngleComputer::ComputeY() const {
      double y(0.);
      const VirtualParticle* scattered = mEvent.ScatteredLepton();
      if(scattered) {
         double theta = scattered->GetTheta();
         double gamma = ComputeQuarkAngle();
         double denominator = tan(theta / 2.) + tan(gamma / 2.);
         if(denominator > 0.) {
            y = tan(gamma / 2.) / denominator;
         } // if
      } // if
      // Return y bounded by the range [0, 1].
      return bounded(y, 0., 1.);
   }
   //	==========================================================================
   //	==========================================================================
   Double_t DoubleAngleComputer::ComputeQSquared() const {
      double Q2(0.);
      const VirtualParticle* lepton = mEvent.BeamLepton();
      const VirtualParticle* scattered = mEvent.ScatteredLepton();
      if(lepton and scattered) {
         double theta = scattered->GetTheta();
         double gamma = ComputeQuarkAngle();
         double denominator = tan(theta / 2.) + tan(gamma / 2.);
         if(denominator > 0.) {
            Q2 = 4. * pow(lepton->GetE(), 2.) / tan(theta / 2.) / denominator;
         } // if
      } // if
      // Return Q^2, requiring it to be positive.
      return std::max(0., Q2);
   }
   //	==========================================================================
   //	==========================================================================
   Double_t DoubleAngleComputer::ComputeX() const {
      double x(0.);
      const VirtualParticle* lepton = mEvent.BeamLepton();
      const VirtualParticle* hadron = mEvent.BeamHadron();
      if(lepton and hadron) {
         double s = (lepton->Get4Vector() + hadron->Get4Vector()).M2();
         double y = ComputeY();
         if(s > 0. and y > 0.) {
            x = ComputeQSquared() / y / s;
         } // if
      } // if
      // Return x bounded by the range [0, 1].
      return bounded(x, 0., 1.);
   }
} // namespace erhic
