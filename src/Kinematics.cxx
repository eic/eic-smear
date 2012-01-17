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

#include "Kinematics.h"

// -------------------------------------------------------------------------------------------------
// Implementation of the KinematicsFromHadrons class follows
// -------------------------------------------------------------------------------------------------

//	================================================================================================
//	================================================================================================
KinematicsFromHadrons::KinematicsFromHadrons()
//: mLeptonEnergy( 0.0 )
//, mMandelstamS( 0.0 )
{ }

//	================================================================================================
//	================================================================================================
//void KinematicsFromHadrons::setLeptonEnergy( const Double_t&){// e ) {
//   mLeptonEnergy = e;
//   std::cerr << "Please use setBeamLepton()" << std::endl;
//}

//	================================================================================================
//	================================================================================================
//void KinematicsFromHadrons::setMandelstamS( const Double_t&) {//s ) {
//   std::cerr << "Please use setBeamLepton() and setBeamHadron()" << std::endl;
   //   mMandelstamS = s;
   //}

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
//   return mLeptonEnergy;
   return mBeamLepton.E();
}

//	================================================================================================
//	================================================================================================
Double_t KinematicsFromHadrons::getMandelstamS() const {
   return 4. * mBeamLepton.E() * mBeamHadron.E();
//   return mMandelstamS;
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
// Using P instead of E here throws off Double angle method
//      return fourMomentum.GetP() - fourMomentum.GetPz();
   }
};

JacquetBlondel::JacquetBlondel()
//: mProtonBeam(NULL)
//, mElectronBeam(NULL) { }
//: mBeamLepton
{}


//	================================================================================================
// y_JB = sum_h{ E_h - pz_h } / 2E_e
//	================================================================================================
Double_t JacquetBlondel::computeY() const {
   
//   if(mProtonBeam and mElectronBeam) return computeYExact();
//   if(mExact and mProtonBeam and mElectronBeam) return computeYExact();
   /*
   using namespace std;
   
   list<double> eMinusPZ;
   
   transform(mParticles.begin(),
             mParticles.end(),
             back_inserter( eMinusPZ ),
             EMinusGetPz() );
   
   return accumulate(eMinusPZ.begin(), eMinusPZ.end(), 0.0)
          / 2.0 / getLeptonEnergy();
    */
   return computeYExact();
}

Double_t JacquetBlondel::computeYExact() const {
//   const TLorentzVector& p = mBeamHadron;//mProtonBeam->Get4Vector();
//   const TLorentzVector& k = mBeamLepton;//mElectronBeam->Get4Vector();
   /*
   TLorentzVector pPrime = std::accumulate(mParticles.begin(),
                                           mParticles.end(),
                                           TLorentzVector(0., 0., 0., 0.));
   return (p.Dot(pPrime) - p.Dot(p)) / p.Dot(k);
    */
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
             /*
   return
   accumulate(eMinusPZ.begin(),
              eMinusPZ.end(),
              -::pow(mBeamHadron.GetM(), 2.)) /
   (mBeamLepton.E() * mBeamHadron.E() +
    mElectronBeam->GetP() * mProtonBeam->GetP());
              */
}

//	================================================================================================
// Q2_JB = [ ( sum_h{ px_h } )^2 + ( sum_h{ py_h } )^2 ] / [ 1 - y_JB ]
// Q2_JB = [ ( sum_h{ pT_h } )^2 ] / [ 1 - y_JB ]
//	================================================================================================
Double_t JacquetBlondel::computeQSquared() const {
   
//   if(mProtonBeam and mElectronBeam) return computeQSquaredExact();
   using namespace std;
   /*
   
   list<double> px, py;
    
   transform( mParticles.begin(),
             mParticles.end(),
             back_inserter( px ),
             mem_fun_ref( &TLorentzVector::Px ) );
   
   transform( mParticles.begin(),
             mParticles.end(),
             back_inserter( py ),
             mem_fun_ref( &TLorentzVector::Py ) );
   
   double sumOfPx = accumulate( px.begin(), px.end(), 0.0 );
   double sumOfPy = accumulate( py.begin(), py.end(), 0.0 );
   
   return ( pow( sumOfPx, 2.0 ) + pow( sumOfPy, 2.0 ) ) / ( 1.0 - computeY() );
   */
   
   // TLorentzVector::Pt is overloaded, so specify which form:
   /*
   Double_t (TLorentzVector::*getPt)() const = &TLorentzVector::Pt;
   
   list<double> pt;
   transform(mParticles.begin(),
             mParticles.end(),
             back_inserter(pt),
             mem_fun_ref(getPt));
   double sumOfPt = accumulate(pt.begin(), pt.end(), 0.);
   
   return pow(sumOfPt, 2.) / (1. - computeY());
   */
   return computeQSquaredExact();
}

// We use the following exact expression:
// 2(E_p * sum(E_h) - pz_p * sum(pz_h)) - sum(E_h)^2 + sum(p_h)^2 - M_p^2
Double_t JacquetBlondel::computeQSquaredExact() const {
   using namespace std;
   
   /*
   // Sum the energies of the hadrons
   list<double> E;
   transform(mParticles.begin(),
             mParticles.end(),
             back_inserter(E),
             mem_fun_ref(&TLorentzVector::E));
             //             mem_fun_ref(&TLorentzVector::P));
   const double sumEh = accumulate(E.begin(), E.end(), 0.);
//   std::cout << "Sum of E = " << sumEh << std::endl;
   // Sum the pz of the hadrons
   list<double> pz;
   transform(mParticles.begin(),
             mParticles.end(),
             back_inserter(pz),
             mem_fun_ref(&TLorentzVector::Pz));
   const double sumPzh = accumulate(pz.begin(), pz.end(), 0.);
    */
//   std::cout << "Sum of pz = " << sumPzh   << std::endl;
   // Vector-sum the total momenta of the hadrons
   /*
   list<TVector3> ph;
   transform(mParticles.begin(),
             mParticles.end(),
             back_inserter(ph),
             mem_fun_ref(&TLorentzVector::Vect));
   const TVector3 sumPh = accumulate(ph.begin(),
                                     ph.end(),
                                     TVector3(0., 0., 0.));
    */
   TLorentzVector total = accumulate(mParticles.begin(),
                                     mParticles.end(),
                                     TLorentzVector());
   
//   return
//   2. * (mBeamHadron.GetP() * sumEh - mBeamHadron.GetPz() * sumPzh) -
//   total.M2() -
//   mBeamHadron.M2();
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
: mHasChanged( false )
, mAngle( NAN )
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

#if 0
namespace erhic {

   
   //=========================================================================
   //
   //=========================================================================
   LeptonKinematics::LeptonKinematics() { }
   
   LeptonKinematics::~LeptonKinematics() { }
   
   Double_t LeptonKinematics::ComputeQSquared(const EventBase& event
                                              , Getter getter
                                              ) const {
      // Locate the scattered lepton
      double Q2(NAN);
      if(event.BeamLepton() and event.ScatteredLepton()) {
         Q2 = 2. * (event.BeamLepton()->*getter)() *
         (event.ScatteredLepton()->*getter)() *
         (1. + event.ScatteredLepton()->Get4Vector().CosTheta());
      } // if
      return Q2;
   }
   
   Double_t LeptonKinematics::ComputeX(const EventBase& event
                                       , Getter getter
                                       ) const {
      double x(NAN);
      const ::Particle* p = event.BeamHadron();
      if(p) {
         x = ComputeQSquared(event) / (2. * p->GetGetM() * ComputeNu(event));
      } // for
      return x;
   }
   
   Double_t LeptonKinematics::ComputeY(const EventBase& event
                                       , Getter getter
                                       ) const {
      double y(NAN);
      const ::Particle* p = event.BeamHadron();
      const ::Particle* l = event.BeamLepton();
      if(l) {
         double gamma = p->Get4Vector().Gamma();
         double beta = p->Get4Vector().Beta();
         // Calculate incident lepton energy in the hadron rest frame.
         double EInNucl = gamma * ((l->*getter)() - beta * l->Get4Vector().GetPz());
         y = ComputeNu(event) / EInNucl;
      } // if
      return y;
   }
   
   Double_t LeptonKinematics::ComputeWSquared(const EventBase& event
                                              , Getter getter
                                              ) const {
      double W2(NAN);
      const ::Particle* p = event.BeamHadron();
      if(p) {
         return p->Get4Vector().M2() +
         (1. - ComputeX(event)) * ComputeQSquared(event) / ComputeX(event);
      } // if
      return W2;
   }
   
   Double_t LeptonKinematics::ComputeNu(const EventBase& event
                                        , Getter getter
                                        ) const {
      double nu(NAN);
      const ::Particle* p = event.BeamHadron();
      const ::Particle* l = event.BeamLepton();
      const ::Particle* s = event.ScatteredLepton();
      if(p and l and s) {
         double gamma = p->Get4Vector().Gamma();
         double beta = p->Get4Vector().Beta();
         double ELeptonInNucl = gamma * ((l->*getter)() - beta * l->Get4Vector().GetPz());
         double ELeptonOutNucl = gamma * ((s->*getter)() - beta * s->Get4Vector().GetPz());
         nu = ELeptonInNucl - ELeptonOutNucl;
      } // if
#if 0
      const ::Particle* p = event.BeamHadron();
      const ::Particle* q = event.ExchangeBoson();
      if(p and q) {
         nu = p->Get4Vector().Dot(q->Get4Vector()) / p->Get4Vector().GetM();
      } // if
#endif
      
      return nu;
   }
   
   
   //=========================================================================
   //
   //=========================================================================
   JacquetBlondel::JacquetBlondel() { }
   
   JacquetBlondel::~JacquetBlondel() { }
   
   Double_t JacquetBlondel::ComputeX(const EventBase& event, Getter) const {
      double y = ComputeY(event);
      double s = 4. * event.BeamHadron()->E() * event.BeamLepton()->E();
      if( y > 0.0 ) return ComputeQSquared(event) / y / s;
      return 0.0;
   }
   
   Double_t JacquetBlondel::ComputeQSquared(const EventBase& event, Getter) const {
      TLorentzVector total = event.HadronicFinalStateMomentum();
      TLorentzVector mBeamHadron = event.BeamHadron()->Get4Vector();
      TLorentzVector mBeamLepton = event.BeamLepton()->Get4Vector();
      return 2. * (mBeamHadron.E() * total.E() - mBeamHadron.GetPz() * total.GetPz())
      - total.M2() - mBeamHadron.M2();
   }
   
   Double_t JacquetBlondel::ComputeY(const EventBase& event, Getter) const {
      using namespace std;
//      cout<<"ComputeY(&event)"<<std::endl;
      vector<const Particle*> mParticles;
      event.HadronicFinalState(mParticles);
      // Sum the energies of the hadrons
      list<double> E;
      transform(mParticles.begin(),
                mParticles.end(),
                back_inserter(E),
                mem_fun(&Particle::GetE));
      //             mem_fun_ref(&TLorentzVector::P));
      const double sumEh = accumulate(E.begin(), E.end(), 0.);
      
      // Sum the pz of the hadrons
      list<double> pz;
      transform(mParticles.begin(),
                mParticles.end(),
                back_inserter(pz),
                mem_fun(&Particle::GetPz));
      const double sumPzh = accumulate(pz.begin(), pz.end(), 0.);
      
      TLorentzVector mBeamHadron = event.BeamHadron()->Get4Vector();
      TLorentzVector mBeamLepton = event.BeamLepton()->Get4Vector();
      return
      (mBeamHadron.E() * sumEh -
       mBeamHadron.GetPz() * sumPzh -
       mBeamHadron.M2())
      /
      (mBeamLepton.E() * mBeamHadron.E() -
       mBeamLepton.GetPz() * mBeamHadron.GetPz());
   }
   
   Double_t JacquetBlondel::ComputeWSquared(const EventBase& event, Getter) const {
      TLorentzVector p = event.BeamHadron()->Get4Vector();
      return p.M2() + ComputeQSquared(event) * (1. / ComputeX(event) - 1.);
   }
   
   Double_t JacquetBlondel::ComputeNu(const EventBase& event, Getter) const { return NAN; }
   
   
   //=========================================================================
   //
   //=========================================================================
   DoubleAngle::DoubleAngle() { }
   
   DoubleAngle::~DoubleAngle() { }
   
   Double_t DoubleAngle::ComputeX(const EventBase& event, Getter) const { return NAN; }
   Double_t DoubleAngle::ComputeQSquared(const EventBase& event, Getter) const { return NAN; }
   Double_t DoubleAngle::ComputeY(const EventBase& event, Getter) const { return NAN; }
   Double_t DoubleAngle::ComputeWSquared(const EventBase& event, Getter) const { return NAN; }
   Double_t DoubleAngle::ComputeNu(const EventBase& event, Getter) const { return NAN; }
   
} // namespace

#endif
