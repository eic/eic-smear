//
// Kinematics.cxx
//
// Created by TB on 7/7/11.
// Copyright 2011 BNL. All rights reserved.
//

#include <cmath>
#include <functional>
#include <iostream>
#include <numeric>

#include <TDatabasePDG.h>
#include <TVector3.h>

#include "eicsmear/erhic/BeamParticles.h"
#include "eicsmear/erhic/Kinematics.h"
#include "eicsmear/erhic/ParticleIdentifier.h"

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
   LeptonKinematicsComputer::LeptonKinematicsComputer(const EventDis& event) {
      ParticleIdentifier::IdentifyBeams(event, mBeams);
   }
   LeptonKinematicsComputer::LeptonKinematicsComputer(const BeamParticles& b)
   : mBeams(b) {
   }
   DisKinematics* LeptonKinematicsComputer::Calculate() {
      DisKinematics* kin(NULL);
      const TLorentzVector& l = mBeams.BeamLepton();
      const TLorentzVector& h = mBeams.BeamHadron();
      const TLorentzVector& s = mBeams.ScatteredLepton();
      kin = new DisKinematics;
      kin->mQ2 = 2. * l.E() * s.E() * (1. + s.CosTheta());
      double ELeptonInNucl = h.Gamma() * (l.E() - h.Beta() * l.Pz());
      double ELeptonOutNucl = h.Gamma() * (s.E() - h.Beta() * s.Pz());
      kin->mNu = ELeptonInNucl - ELeptonOutNucl;
      kin->mX = kin->mQ2 / (2. * h.M() * kin->mNu);
      kin->mY = kin->mNu / ELeptonInNucl;
      kin->mW2 = h.M2() + (1. - kin->mX) * kin->mQ2 / kin->mX;
      return kin;
   }
   JacquetBlondelComputer::JacquetBlondelComputer(const EventDis& event,
                                                  const BeamParticles* beams)
   : mEvent(event) {
      if(beams) {
         mBeams = *beams;
      } // if
      else {
         ParticleIdentifier::IdentifyBeams(event, mBeams);
      } // else
   }
   DisKinematics* JacquetBlondelComputer::Calculate() {
      DisKinematics* kin(NULL);
      ::JacquetBlondel jacquetBlondel;
      jacquetBlondel.setBeamLepton(mBeams.BeamLepton());
      jacquetBlondel.setBeamHadron(mBeams.BeamHadron());
      std::vector<const VirtualParticle*> final;
      mEvent.HadronicFinalState(final);
      for(unsigned i(0); i < final.size(); ++i ) {
         jacquetBlondel.addParticle(final.at(i)->Get4Vector() );
      } // for
      kin = new DisKinematics;
      kin->mY = jacquetBlondel.computeY();
      kin->mQ2 = jacquetBlondel.computeQSquared();
      kin->mX = jacquetBlondel.computeX();
      kin->mW2 = mBeams.BeamHadron().M2() + (1. - kin->mX) *
                 kin->mQ2 / kin->mX;
      return kin;
   }
   DoubleAngleComputer::DoubleAngleComputer(const EventDis& event,
                                            const BeamParticles* beams)
   : mEvent(event) {
      if(beams) {
         mBeams = *beams;
      } // if
      else {
      ParticleIdentifier::IdentifyBeams(event, mBeams);
      } // else
   }
   DisKinematics* DoubleAngleComputer::Calculate() {
      DisKinematics* kin(NULL);
      ::DoubleAngle doubleAngle;
      doubleAngle.setBeamLepton(mBeams.BeamLepton());
      doubleAngle.setBeamHadron(mBeams.BeamHadron());
      doubleAngle.setLeptonAngle(mBeams.ScatteredLepton().Theta());
      std::vector<const VirtualParticle*> final;
      mEvent.HadronicFinalState(final);
      for(unsigned i(0); i < final.size(); ++i ) {
         doubleAngle.addParticle(final.at(i)->Get4Vector() );
      } // for
      kin = new DisKinematics;
      kin->mY = doubleAngle.computeY();
      kin->mQ2 = doubleAngle.computeQSquared();
      kin->mX = doubleAngle.computeX();
      kin->mW2 = mBeams.BeamHadron().M2() + (1. - kin->mX) *
                 kin->mQ2 / kin->mX;
      return kin;
   }
} // namespace erhic

// -------------------------------------------------------------------------------------------------
// Implementation of the KinematicsFromHadrons class follows
// -------------------------------------------------------------------------------------------------

//	================================================================================================
//	================================================================================================
KinematicsFromHadrons::KinematicsFromHadrons()
{ }

//	================================================================================================
//	================================================================================================
void KinematicsFromHadrons::addParticle( const TLorentzVector& p ) {
   mParticles.push_back( p );
}

//	================================================================================================
//	================================================================================================
void KinematicsFromHadrons::clearParticles() {
   mParticles.clear();
}

//	================================================================================================
//	================================================================================================
Double_t KinematicsFromHadrons::getLeptonEnergy() const {
   return mBeamLepton.E();
}

//	================================================================================================
//	================================================================================================
Double_t KinematicsFromHadrons::getMandelstamS() const {
   return 4. * mBeamLepton.E() * mBeamHadron.E();
}

// -------------------------------------------------------------------------------------------------
// Implementation of the JacquetBlondel class follows
// -------------------------------------------------------------------------------------------------

//	================================================================================================
// Helper functor for calculating E_h - pZ from a TLorentzVector.
//	================================================================================================
struct EMinusPz : public std::unary_function<TLorentzVector, double> {
   double operator()( const TLorentzVector& fourMomentum ) const {
      return fourMomentum.E() - fourMomentum.Pz();
   }
};

JacquetBlondel::JacquetBlondel()
{}


//	================================================================================================
// y_JB = sum_h{ E_h - pz_h } / 2E_e
//	================================================================================================
Double_t JacquetBlondel::computeY() const {
   return computeYExact();
}

Double_t JacquetBlondel::computeYExact() const {
   using namespace std;
   // Sum the energies of the hadrons
   list<double> E;
   transform(mParticles.begin(),
             mParticles.end(),
             back_inserter(E),
             mem_fun_ref(&TLorentzVector::E));
             //             mem_fun_ref(&TLorentzVector::P));
   const double sumEh = accumulate(E.begin(), E.end(), 0.);
   
   // Sum the pz of the hadrons
   list<double> pz;
   transform(mParticles.begin(),
             mParticles.end(),
             back_inserter(pz),
             mem_fun_ref(&TLorentzVector::Pz));
   const double sumPzh = accumulate(pz.begin(), pz.end(), 0.);
   return
   (mBeamHadron.E() * sumEh -
    mBeamHadron.Pz() * sumPzh -
    mBeamHadron.M2())
   /
   (mBeamLepton.E() * mBeamHadron.E() -
    mBeamLepton.Pz() * mBeamHadron.Pz());
}

//	================================================================================================
// Q2_JB = [ ( sum_h{ px_h } )^2 + ( sum_h{ py_h } )^2 ] / [ 1 - y_JB ]
// Q2_JB = [ ( sum_h{ pT_h } )^2 ] / [ 1 - y_JB ]
//	================================================================================================
Double_t JacquetBlondel::computeQSquared() const {
   using namespace std;
   return computeQSquaredExact();
}

// We use the following exact expression:
// 2(E_p * sum(E_h) - pz_p * sum(pz_h)) - sum(E_h)^2 + sum(p_h)^2 - M_p^2
Double_t JacquetBlondel::computeQSquaredExact() const {
   using namespace std;
   TLorentzVector total = accumulate(mParticles.begin(),
                                     mParticles.end(),
                                     TLorentzVector());
   return 2. * (mBeamHadron.E() * total.E() - mBeamHadron.Pz() * total.Pz())
   - total.M2() - mBeamHadron.M2();
}

//	================================================================================================
// x_JB = Q2_JB / ( y_JB * s )
//	================================================================================================
Double_t JacquetBlondel::computeX() const {
   double y = computeY();
   if( y > 0.0 ) return computeQSquared() / y / getMandelstamS();
   return 0.0;
}

// -------------------------------------------------------------------------------------------------
// Implementation of the DoubleAngle class follows
// -------------------------------------------------------------------------------------------------

//	================================================================================================
// 
//	================================================================================================
DoubleAngle::DoubleAngle()
: mHasChanged(false)
, mAngle(NAN)
{ }

//	================================================================================================
// 
//	================================================================================================
void DoubleAngle::addParticle( const TLorentzVector& p ) {
   KinematicsFromHadrons::addParticle( p );
   mHasChanged = true;
}

//	================================================================================================
// Scattering angle of struck quark
// cos(angle) = A / B
// where A = (sum of px_h)^2 + (sum of py_h)^2 - (sum of [E_h - pz_h])^2
// and   B = (sum of px_h)^2 + (sum of py_h)^2 + (sum of [E_h - pz_h])^2
// This is a computationally expensive operation that is called a lot, so cashe the result
// until the particle list changes.
//	================================================================================================
Double_t DoubleAngle::computeQuarkAngle() const {
   
   using namespace std;
   
   if( not mHasChanged ) return mAngle;
   
   list<double> px, py, eMinusPZ;
   
   // Get the px of each particle:
   transform( mParticles.begin(), mParticles.end(), back_inserter( px ),
              mem_fun_ref( &TLorentzVector::Px ) );
   // Get the py of each particle:
   transform( mParticles.begin(), mParticles.end(), back_inserter( py ),
              mem_fun_ref( &TLorentzVector::Py ) );
   // Get the (E-pz) of each particle:
   transform( mParticles.begin(), mParticles.end(), back_inserter( eMinusPZ ),
              EMinusPz() );
   
   // Square the sum of each momentum component
   struct SquareOfSum {
      double operator()( const list<double>& l ) const {
         return pow( accumulate( l.begin(), l.end(), 0.0 ), 2.0 );
      }
   } squareOfSum;
   
   double sumOfPx2 = squareOfSum( px );
   double sumOfPy2 = squareOfSum( py );
   double sumOfEMinusPz = squareOfSum( eMinusPZ );
   
   // Calculate the angle:
   double denominator = ( sumOfPx2 + sumOfPy2 + sumOfEMinusPz );
   // Set cos(angle) to 1 by default. If denominator is zero or less, angle will
   // then evaluate to acos(1) = 0.
   double cosAngle( 1.0 );
   if( denominator > 0.0 ) {
      cosAngle = ( sumOfPx2 + sumOfPy2 - sumOfEMinusPz )
      / denominator;
   } // if
   
   mAngle = TMath::ACos( cosAngle );
   
   mHasChanged = false;
   return mAngle;
}

//	================================================================================================
// y_DA = Q2_DA / ( x_DA * s )
//	================================================================================================
Double_t DoubleAngle::computeY() const {
   double x = computeX();
   if( x > 0.0 ) return computeQSquared() / x / getMandelstamS();
   return 0.0;
}

//	================================================================================================
// a_q = quark scattering angle, a_e = electron scattering angle
// Q2_DA = 4 * E_e^2 * sin(a_q) ( 1 + cos(a_e) ) / [ sin(a_q) + sin(a_e) - sin( a_q + a_e ) ]
//	================================================================================================
Double_t DoubleAngle::computeQSquared() const {
   
   double numerator = 4.0
   * pow( getLeptonEnergy(), 2.0 )
   * sin( computeQuarkAngle() )
   * ( 1.0 + cos( getLeptonAngle() ) );
   
   double denominator = sin( computeQuarkAngle() )
   + sin( getLeptonAngle() )
   - sin( computeQuarkAngle() + getLeptonAngle() );
   
   if( denominator > 0.0 ) return numerator / denominator;
   return 0.0;
}

//	================================================================================================
// x_DA = (E_e/E_p)[sin(a_q)+sin(a_e)+sin(a_q+a_e)]/[sin(a_q)+sin(a_e)-sin(a_q+a_e)]
//	================================================================================================
Double_t DoubleAngle::computeX() const {
   double numerator = sin( computeQuarkAngle() )
   + sin( getLeptonAngle() )
   + sin( computeQuarkAngle() + getLeptonAngle() );
   
   double denominator = sin( computeQuarkAngle() )
   + sin( getLeptonAngle() )
   - sin( computeQuarkAngle() + getLeptonAngle() );
   
   if( denominator > 0.0 ) {
      return getLeptonEnergy() / getProtonEnergy() * numerator / denominator;
   } // if
   return 0.0;
}
