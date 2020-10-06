/**
 \file
 Implementation of class Smear::Detector.

 \author    Michael Savastio
 \date      2011-08-19
 \copyright 2011 Brookhaven National Lab
 */

#include "eicsmear/smear/Detector.h"

#include <algorithm>
#include <functional>
#include <list>
#include <memory>
#include <vector>

#include "eicsmear/erhic/EventDis.h"
#include "eicsmear/smear/EventSmear.h"
#include "eicsmear/erhic/Kinematics.h"
#include "eicsmear/smear/ParticleMCS.h"
#include "eicsmear/smear/Smearer.h"
#include "eicsmear/erhic/VirtualParticle.h"

using std::cout;
using std::cerr;
using std::endl;

namespace Smear {

Detector::Detector()
: useNM(false)
, useJB(false)
, useDA(false) {
}

Detector::Detector(const Detector& other)
: TObject(other) {
  useNM = other.useNM;
  useJB = other.useJB;
  useDA = other.useDA;
  Devices = other.CopyDevices();
  LegacyMode = other.GetLegacyMode();
}

Detector& Detector::operator=(const Detector& that) {
  if (this != &that) {
    useNM = that.useNM;
    useJB = that.useJB;
    useDA = that.useDA;
    Devices = that.CopyDevices();
    LegacyMode = that.GetLegacyMode();
  }  // if
  return *this;
}

Detector::~Detector() {
  DeleteAllDevices();
}

void Detector::DeleteAllDevices() {
  for (unsigned i(0); i < GetNDevices(); i++) {
    delete Devices.at(i);
    Devices.at(i) = NULL;
  }  // for
  Devices.clear();
}

void Detector::AddDevice(Smearer& dev) {
  Devices.push_back(dev.Clone());
}

void Detector::SetEventKinematicsCalculator(TString s) {
  s.ToLower();
  useNM = s.Contains("nm") || s.Contains("null");
  useJB = s.Contains("jb") || s.Contains("jacquet");
  useDA = s.Contains("da") || s.Contains("double");
}

Smearer* Detector::GetDevice(int n) {
  Smearer* smearer(NULL);
  if (unsigned(n) < Devices.size()) {
    smearer = Devices.at(n);
  }  // if
  return smearer;
}

void Detector::FillEventKinematics(Event* eventS) {
  if (!(useNM || useJB || useDA)) {
    return;
  }  // if
     // Need a bit of jiggery-pokery here, as the incident beam info isn't
     // associated with the smeared event.
     // So, get the beam info from the MC event, but replace the scattered
     // electron with the smeared version.
     // Then we can use the standard JB/DA algorithms on the smeared event.
  const ParticleMCS* scattered = eventS->ScatteredLepton();
  typedef std::unique_ptr<erhic::DisKinematics> KinPtr;
  if (useNM && scattered) {
    KinPtr kin(erhic::LeptonKinematicsComputer(*eventS).Calculate());
    if (kin.get()) {
      eventS->SetLeptonKinematics(*kin);
    }  // if
  } else {
    eventS->SetLeptonKinematics( erhic::DisKinematics(-1., -1., -1., -1., -1.));
  }  // if
  if (useJB) {
    KinPtr kin(erhic::JacquetBlondelComputer(*eventS).Calculate());
    if (kin.get()) {
      eventS->SetJacquetBlondelKinematics(*kin);
    }  // if
  }  // if
  if (useDA && scattered) {
    KinPtr kin(erhic::DoubleAngleComputer(*eventS).Calculate());
    if (kin.get()) {
      eventS->SetDoubleAngleKinematics(*kin);
    }  // if
  }  // if
}

std::list<Smearer*> Detector::Accept(const erhic::VirtualParticle& p) const {
  std::list<Smearer*> devices;
  // Only accept final-state particles, so skip the check against each
  // devices for non-final-state particles.
  if (p.GetStatus() == 1) {
    std::vector<Smearer*>::const_iterator iter;
    for (iter = Devices.begin(); iter != Devices.end(); ++iter) {
      // Store each device that accepts the particle.
      if ((*iter)->Accept.Is(p)) {
        devices.push_back(*iter);
      }  // if
    }  // for
  }  // if
  return devices;
}

ParticleMCS* Detector::Smear(const erhic::VirtualParticle& prt) const {
  // Does the particle fall in the acceptance of any device?
  // If so, we smear it, if not, we skip it (store a NULL pointer).
  std::list<Smearer*> devices = Accept(prt);
  ParticleMCS* prtOut(NULL);
  if (!devices.empty()) {
    // It passes through at least one device, so smear it.
    // Devices in which it doesn't pass won't smear it.
    prtOut = new ParticleMCS();
    prtOut->SetSmeared();
    std::list<Smearer*>::iterator iter;
    for (iter = devices.begin(); iter != devices.end(); ++iter) {
      (*iter)->Smear(prt, *prtOut);
    }  // for
    
    if (LegacyMode){
      // Compute derived momentum components.
      prtOut->SetPx( prtOut->GetP() * sin(prtOut->GetTheta()) * cos(prtOut->GetPhi()));
      prtOut->SetPy( prtOut->GetP() * sin(prtOut->GetTheta()) * sin(prtOut->GetPhi()));
      prtOut->SetPt( sqrt(pow(prtOut->GetPx(), 2.) + pow(prtOut->GetPy(), 2.)));
      prtOut->SetPz( prtOut->GetP() * cos(prtOut->GetTheta()));
      
    } else { // not in LegacyMode
      // Sanity check and computation of derived momentum components.
      // ---------------------------------------------------------------------
      // - this cannot happen...
      if ( !prtOut->IsSmeared() ) throw std::runtime_error ("particle seems to be not smeared?!");
      
      // Is momentum smeared at all?
      int MomComponentsChanged = prtOut->IsPSmeared() + prtOut->IsPtSmeared() + prtOut->IsPxSmeared() + prtOut->IsPySmeared() + prtOut->IsPzSmeared();
      int AngComponentsChanged = prtOut->IsPhiSmeared() + prtOut->IsThetaSmeared();
      
      if ( MomComponentsChanged==0 ){
	// Momentum is untouched (pure calo measurement)
	// all we have to do is ensure phi and theta are explicitly smeared
	if ( !prtOut->IsPhiSmeared() ) {
	  cerr << "Phi always needs to be smeared (at least with sigma=0)" << endl;
	  cerr << "For legacy smear scripts, use det.SetLegacyMode ( true );" << endl;
	  throw std::runtime_error ("Failed consistency check in Detector::Smear()");
	}
	if ( !prtOut->IsThetaSmeared() ){
	  cerr << "Theta always needs to be smeared (at least with sigma=0)" << endl;
	  cerr << "For legacy smear scripts, use det.SetLegacyMode ( true );" << endl;
	  throw std::runtime_error ("Failed consistency check in Detector::Smear()");
	}
      } else if ( MomComponentsChanged + AngComponentsChanged != 3){
	// - Momentum is exactly defined by three independent quantities, including theta and phi
	cerr << "Expected 0 (excluding phi, theta) or exactly 3 (excluding phi, theta) smeared momentum quantities." << endl;
	cerr << "For legacy smear scripts, use det.SetLegacyMode ( true );" << endl;
	cerr << "MomComponentsChanged = " << MomComponentsChanged << endl;
	cerr << prt.GetEta() << endl;
	cerr << prtOut->IsPSmeared() << endl;
	cerr << prtOut->IsPtSmeared() << endl;
	cerr << prtOut->IsPxSmeared() << endl;
	cerr << prtOut->IsPySmeared() << endl;
	cerr << prtOut->IsPzSmeared() << endl;	
	cerr << "AngComponentsChanged = " << AngComponentsChanged << endl;	
	cerr << " pid : " << prt.Id() << endl;
	throw std::runtime_error ("Failed consistency check in Detector::Smear()");
      } else {
	// We now have exactly three out of P, px, py, pz, pt, phi, theta. Compute the rest.
	// That's 35 total cases, but luckily many of them are not independent, or nonsensical in a detector.
	// NOTE: We need p^2 >= pT^2, pz^2.
	// Smearing (P and Pz) or (P and pT) is obscure enough to where we just make the ad-hoc decision
	// to adjust P in such a case.
	
	// - px, py, pt, phi are intimately related. Let's condense.
	if ( prtOut->IsPxSmeared() ^ prtOut->IsPySmeared() ) {// "^" = XOR
	  cerr << "Smearing only one out of px, py is not supported. Please contact the authors." << endl;
	  cerr << "For legacy smear scripts, use det.SetLegacyMode ( true );" << endl;
	  throw std::runtime_error ("Failed consistency check in Detector::Smear()");
	} // illegal px, py --> removes 2 * ( 5 choose 2 ) = 20 combinations
	
	if ( prtOut->IsPxSmeared() && prtOut->IsPySmeared() ) {
	  if ( prtOut->IsPhiSmeared() ) {
	    cerr << "Smearing px, py, and phi is inconsistent" << endl;
	    cerr << "For legacy smear scripts, use det.SetLegacyMode ( true );" << endl;
	    throw std::runtime_error ("Failed consistency check in Detector::Smear()");
	  }// 1, running 21
	  if ( prtOut->IsPtSmeared() ) {
	    cerr << "Smearing px, py, and pt is inconsistent" << endl;
	    cerr << "For legacy smear scripts, use det.SetLegacyMode ( true );" << endl;
	    throw std::runtime_error ("Failed consistency check in Detector::Smear()");	    
	  }  // 1, running 22
	  
	  // Can set phi and pt now
	  prtOut->SetPhi( std::atan2( prtOut->GetPy(),prtOut->GetPx() ) );
	  prtOut->SetPt( std::sqrt( std::pow( prtOut->GetPx(),2) + std::pow( prtOut->GetPy(),2) ) );
	  
	  // legal options remaining: pz, P, theta
	  if ( prtOut->IsPzSmeared() ){ // pz is smeared
	    prtOut->SetTheta( std::atan2(prtOut->GetPt() ,prtOut->GetPz() ) );
	    prtOut->SetP( std::sqrt( std::pow( prtOut->GetPt(),2) + std::pow( prtOut->GetPz(),2) ) );
	  } // 1, running 23
	  
	  if ( prtOut->IsPSmeared() ){ // p is smeared
	    if ( prtOut->GetP() < prtOut->GetPt() ) prtOut->SetP(prtOut->GetPt(), false);
	    prtOut->SetTheta( std::atan2(prtOut->GetPt() ,prtOut->GetPz() ) );
	    prtOut->SetPz( std::sqrt( std::pow(prtOut->GetP(), 2.) - std::pow(prtOut->GetPt(), 2.)) );
	  } // 1, running 24
	  
	  if ( prtOut->IsThetaSmeared() ){ // theta is smeared
	    // Note: Here (as in other cases), we rely on the device to deliver physically sound
	    // numbers (see and adapt HandleBogusValues). But for now we explicitly fault on division by zero
	    assert ( fabs(std::tan(prtOut->GetTheta())) > 1e-8 );
	    prtOut->SetPz( prtOut->GetPt() / std::tan(prtOut->GetTheta()) );
	    prtOut->SetP( std::sqrt( std::pow( prtOut->GetPt(),2) + std::pow( prtOut->GetPz(),2) ) );
	  } // 1, running 25
	  
	} // know both px and py --> knocked off 5 cases, 25 running
	
	// We're left with px and py unsmeared, and choose 3 out of P, pt, pz, phi, theta = 10 combinations.
	// But we can reject the case of unsmeared phi since without px, py we can't reconstruct azimuth
	// No need to protect with another if() because we touched Phi in the previous clause
	if ( !prtOut->IsPhiSmeared() ) {
	  cerr << "Momentum components are smeared, but neither phi nor px and py are." << endl;
	  cerr << "For legacy smear scripts, use det.SetLegacyMode ( true );" << endl;
	  throw std::runtime_error ("Failed consistency check in Detector::Smear()");	    
	} // 4 cases (4 choose 3), 29 running
	
	// Leaves (choose 2 out of P, pz, pt, theta) = 6 combinations, math checks out.
	// Using a complementary if clause instead of another else for readability
	if ( !prtOut->IsPxSmeared() && !prtOut->IsPySmeared() ) {
	  // NOTE: Currently, this is the only loop that should matter
	  // since we don't allow px and py smearing in the formulas,
	  // but there's no reason not to allow it in the future
	  if(prtOut->IsPSmeared() && prtOut->IsThetaSmeared() ) /*(1100)*/ { // P, theta
	    
	    prtOut->SetPt ( prtOut->GetP() * std::sin(prtOut->GetTheta()));
	    prtOut->SetPz ( prtOut->GetP() * std::cos(prtOut->GetTheta()));
	    
	  } else if( prtOut->IsPSmeared() && prtOut->IsPtSmeared() ) /*(1010)*/ { // P, pt
	    if ( prtOut->GetP() < prtOut->GetPt() ) prtOut->SetP(prtOut->GetPt(), false);
	    prtOut->SetPz( std::sqrt(std::pow(prtOut->GetP(), 2) - std::pow(prtOut->GetPt(), 2)));
	    prtOut->SetTheta ( std::atan2(prtOut->GetPt(),prtOut->GetPz()));
	    
	  } else if(prtOut->IsPzSmeared() && prtOut->IsPtSmeared() ) /*(0011)*/ { // pt, pz
	    
	    prtOut->SetP( std::sqrt(std::pow(prtOut->GetPt(), 2) + std::pow(prtOut->GetPz(), 2)));
	    prtOut->SetTheta ( std::atan2(prtOut->GetPt(),prtOut->GetPz()));
	    
	  } else if(prtOut->IsThetaSmeared() && prtOut->IsPtSmeared()) /*(0110)*/ { // pt, theta
	    // Note: Here (as in other cases), we rely on the device to deliver physically sound
	    // numbers (see and adapt HandleBogusValues). But for now we explicitly fault on division by zero
	    assert ( fabs(std::tan(prtOut->GetTheta())) > 1e-8 );
	    prtOut->SetPz( prtOut->GetPt() / std::tan(prtOut->GetTheta()) );
	    prtOut->SetP( std::sqrt(std::pow(prtOut->GetPt(), 2) + std::pow(prtOut->GetPz(), 2)));
	    
	  } else if(prtOut->IsPSmeared() && prtOut->IsPzSmeared()) /*(1001)*/ { // P, pz
	    if ( prtOut->GetP() < std::abs(prtOut->GetPz()) ) prtOut->SetP( std::abs(prtOut->GetPz()), false);
	    prtOut->SetPt( std::sqrt(std::pow(prtOut->GetP(), 2) - std::pow(prtOut->GetPz(), 2)));
	    prtOut->SetTheta( std::atan2(prtOut->GetPt() ,prtOut->GetPz() ) );		
	    
	  } else if(prtOut->IsThetaSmeared() && prtOut->IsPzSmeared())  /*(0101)*/ { // theta, pz
	    
	    prtOut->SetPt( prtOut->GetPz() * std::tan(prtOut->GetTheta()));
	    prtOut->SetP( std::sqrt( std::pow( prtOut->GetPt(),2) + std::pow( prtOut->GetPz(),2) ) );
	  }
	  
	  prtOut->SetPx (prtOut->GetP() * std::sin(prtOut->GetTheta()) * std::cos(prtOut->GetPhi()));
	  prtOut->SetPy (prtOut->GetP() * std::sin(prtOut->GetTheta()) * std::sin(prtOut->GetPhi()));
	  
	} // Know px and py are unsmeared, phi IS smeared
	
      } // case treatment for momentum components changed
      
    } // LegacyMode
  } // if devices not empty

  // Done.
  return prtOut;

}

std::vector<Smearer*> Detector::CopyDevices() const {
  std::vector<Smearer*> copies;
  for ( std::vector<Smearer*>::const_iterator it = Devices.begin();
	it!=Devices.end();
	++it ){
    copies.push_back ( (*it)->Clone(""));
  }
  return copies;
}

void Detector::Print(Option_t* o) const {
  for (unsigned i(0); i < GetNDevices(); ++i) {
    Devices.at(i)->Print(o);
  }  // for
}

  void Detector::SetLegacyMode( const bool mode ){
    LegacyMode = mode;
    if ( LegacyMode ){
      std::cout << "Warning: Turning on legacy mode, i.e. deactivating consistency checks and momentum regularization in Smear(). Use only for legacy smear scripts from earlier versions (<~1.0.4)" << endl;
    }
  }
  
  bool Detector::GetLegacyMode() const {
    return LegacyMode;
  }

}  // namespace Smear
