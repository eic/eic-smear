/**
 \file
 Implementation of class Smear::ParticleMCS.
 
 \author    Michael Savastio
 \date      2011-08-19
 \copyright 2011 Brookhaven National Lab
 */

#include "eicsmear/smear/ParticleMCS.h"

#include <iomanip>
#include <iostream>

#include <TMath.h>

namespace Smear {

ParticleMCS::ParticleMCS()
: status(0)
, id(0)
, px(0.)
, py(0.)
, pz(0.)
, E(0.)
, pt(0.)
, p(0.)
, theta(0.)
, phi(0.)
, numSigmaElectron( std::nan("") )
, numSigmaPion( std::nan("") )
, numSigmaProton( std::nan("") )
, numSigmaKaon( std::nan("") )
, numSigmaMuon( std::nan("") )
{
}

ParticleMCS::ParticleMCS(const TLorentzVector& ep, int pdg, int stat)
: status(stat)
, id(pdg)
, px(ep.Px())
, py(ep.Py())
, pz(ep.Pz())
, E(ep.E())
, pt(ep.Pt())
, p(ep.P())
, theta(ep.Theta())
, phi(ep.Phi())
, numSigmaElectron( std::nan("") )
, numSigmaPion( std::nan("") )
, numSigmaProton( std::nan("") )
, numSigmaKaon( std::nan("") )
, numSigmaMuon( std::nan("") )
{
}

ParticleMCS::~ParticleMCS() {
}

TLorentzVector ParticleMCS::Get4Vector() const {
  return TLorentzVector(px, py, pz, E);
}

void ParticleMCS::Print(Option_t* /* unused */) const {
  // Store initial flags for cout so we can restore them later.
  std::ios_base::fmtflags ff = std::cout.flags();
  // Configure flags for output
  std::cout.setf(std::ios_base::fixed, std::ios_base::basefield);
  std::cout.precision(4);
  // Output values.
  std::cout <<
  std::setw(3) << status <<
  std::setw(12) << id <<
  std::setw(10) << px <<
  std::setw(10) << py <<
  std::setw(10) << pz <<
  std::setw(10) << E << std::endl;
  // Restore initial cout flags.
  std::cout.flags(ff);
}

Double_t ParticleMCS::GetEta() const {
  // The default value of -19 is used in the main eRHIC code,
  // so use that for consistency.
  double eta(-19.);
  const double theta = GetTheta();
  if (theta > 0. && theta < TMath::Pi() && !TMath::IsNaN(theta)) {
    eta = -log(tan(theta / 2.));
  }  // if
  return eta;
}

Double_t ParticleMCS::GetRapidity() const {
  // The default value of -19 is used in the main eRHIC code,
  // so use that for consistency.
  double y(-19.);
  // Be careful here, as the energy and momentum can become
  // smeared such that the stored E < pz, giving NaN when
  // taking a log of a negative.
  // In that case, return the default value.
  // E > pz already takes care of the case if E is NaN, so
  // no need to check that again.
  if (E > pz && !TMath::IsNaN(pz)) {
    y = 0.5 * log((E + pz) / (E - pz));
  }  // if
  return y;
}

  // -----------------------------------------------------------
  
  Double_t ParticleMCS::GetPx() const {
    return p * sin(theta) * cos(phi);
  }

  Double_t ParticleMCS::GetPy() const {
    return p * sin(theta) * sin(phi);
  }

  Double_t ParticleMCS::GetPz() const {
    return pz;
  }

  Double_t ParticleMCS::GetE() const {
    return E;
  }

  Double_t ParticleMCS::GetM() const {
    return sqrt(pow(E, 2.) - pow(p, 2.));
  }

  Double_t ParticleMCS::GetPt() const {
    return pt;
  }

  TVector3 ParticleMCS::GetVertex() const {
    return TVector3();
  }

  Double_t ParticleMCS::GetP() const {
    return p;
  }

  Double_t ParticleMCS::GetTheta() const {
    return theta;
  }

  Double_t ParticleMCS::GetPhi() const {
    return phi;
  }

  UShort_t ParticleMCS::GetStatus() const {
    return status;
  }

  double ParticleMCS::GetNumSigmaElectron() const {
    return numSigmaElectron;
  }

  double ParticleMCS::GetNumSigmaPion() const {
    return numSigmaPion;
  }

  double ParticleMCS::GetNumSigmaProton() const {
    return numSigmaProton;
  }

  double ParticleMCS::GetNumSigmaKaon() const {
    return numSigmaKaon;
  }

  double ParticleMCS::GetNumSigmaMuon() const {
    return numSigmaMuon;
  }

  bool ParticleMCS::IsSmeared() const { return kParticleSmeared; }
  bool ParticleMCS::IsESmeared() const { return kESmeared; }
  bool ParticleMCS::IsPSmeared() const { return kPSmeared; }
  bool ParticleMCS::IsPtSmeared() const { return kPtSmeared; }
  bool ParticleMCS::IsPxSmeared() const { return kPxSmeared; }
  bool ParticleMCS::IsPySmeared() const { return kPySmeared; }
  bool ParticleMCS::IsPzSmeared() const { return kPzSmeared; }
  bool ParticleMCS::IsThetaSmeared() const { return kThetaSmeared; }
  bool ParticleMCS::IsPhiSmeared() const { return kPhiSmeared; }
  bool ParticleMCS::IsIdSmeared() const { return kIdSmeared; }
  bool ParticleMCS::IsNumSigmaSmeared() const { return kNumSigmaSmeared; }
  
  // -----------------------------------------------------------
  void ParticleMCS::SetVariable( const double z, const KinType kin) {
    switch (kin) {
    case kE:
      SetE(z); break;
    case kP:
      SetP(z); break;
    case kTheta:
      SetTheta(z); break;
    case kPhi:
      SetPhi(z); break;
    case kPz:
      SetPz(z); break;
    case kPt:
      SetPt(z); break;
    default:
      break;
    }  // switch
  }

  void ParticleMCS::SetE(const Double_t value, const bool CheckSetSmearFlag){
    if ( kESmeared && CheckSetSmearFlag ) throw std::runtime_error ("Attempting to smear E twice");
    else kESmeared = true;
    E = value;
  }

  void ParticleMCS::SetP(Double_t value, const bool CheckSetSmearFlag){
    if ( kPSmeared && CheckSetSmearFlag ) throw std::runtime_error ("Attempting to smear P twice");
    else kPSmeared = true;
    p = value;
  }

  void ParticleMCS::SetPt(Double_t value, const bool CheckSetSmearFlag){
    if ( kPtSmeared && CheckSetSmearFlag ) throw std::runtime_error ("Attempting to smear Pt twice");
    else kPtSmeared = true;
    pt = value;
  }

  void ParticleMCS::SetPx(Double_t value, const bool CheckSetSmearFlag){
    if ( kPxSmeared && CheckSetSmearFlag ) throw std::runtime_error ("Attempting to smear Px twice");
    else kPxSmeared = true;
    px = value;
  }

  void ParticleMCS::SetPy(Double_t value, const bool CheckSetSmearFlag){
    if ( kPySmeared && CheckSetSmearFlag ) throw std::runtime_error ("Attempting to smear Py twice");
    else kPySmeared = true;
    py = value;
  }

  void ParticleMCS::SetPz(Double_t value, const bool CheckSetSmearFlag){
    if ( kPzSmeared && CheckSetSmearFlag ) throw std::runtime_error ("Attempting to smear Pz twice");
    else kPzSmeared = true;
    pz = value;
  }

  void ParticleMCS::SetPhi(Double_t value, const bool CheckSetSmearFlag){
    if ( kPhiSmeared && CheckSetSmearFlag ) throw std::runtime_error ("Attempting to smear Phi twice");
    else kPhiSmeared = true;
    phi = value;
  }

  void ParticleMCS::SetTheta(Double_t value, const bool CheckSetSmearFlag){
    if ( kThetaSmeared && CheckSetSmearFlag ) throw std::runtime_error ("Attempting to smear Theta twice");
    else kThetaSmeared = true;
    theta = value;
  }

  void ParticleMCS::SetId(Int_t i, const bool CheckSetSmearFlag){
    if ( kIdSmeared && CheckSetSmearFlag ) throw std::runtime_error ("Attempting to smear Id twice");
    else kIdSmeared = true;	
    id = i;
  }

  void ParticleMCS::SetStatus(Int_t i) {
    status = i;
  }

  void ParticleMCS::SetNumSigmaElectron( const double d, const bool CheckSetSmearFlag){
    if ( kNumSigmaSmeared && CheckSetSmearFlag ) throw std::runtime_error ("Attempting to smear numSigma twice");
    else kNumSigmaSmeared = true;
    numSigmaElectron = d;
  }

  void ParticleMCS::SetNumSigmaPion( const double d, const bool CheckSetSmearFlag){
    if ( kNumSigmaSmeared && CheckSetSmearFlag ) throw std::runtime_error ("Attempting to smear numSigma twice");
    else kNumSigmaSmeared = true;
    numSigmaPion = d;
  }

  void ParticleMCS::SetNumSigmaProton( const double d, const bool CheckSetSmearFlag){
    if ( kNumSigmaSmeared && CheckSetSmearFlag ) throw std::runtime_error ("Attempting to smear numSigma twice");
    else kNumSigmaSmeared = true;
    numSigmaProton = d;
  }

  void ParticleMCS::SetNumSigmaKaon( const double d, const bool CheckSetSmearFlag){
    if ( kNumSigmaSmeared && CheckSetSmearFlag ) throw std::runtime_error ("Attempting to smear numSigma twice");
    else kNumSigmaSmeared = true;
    numSigmaKaon = d;
  }

  void ParticleMCS::SetNumSigmaMuon( const double d, const bool CheckSetSmearFlag){
    if ( kNumSigmaSmeared && CheckSetSmearFlag ) throw std::runtime_error ("Attempting to smear numSigma twice");
    else kNumSigmaSmeared = true;
    numSigmaMuon = d;
  }
  
  erhic::Pid ParticleMCS::Id() const {
    return ::erhic::Pid(id);
  }

  void ParticleMCS::HandleBogusValues( KinType kin ) {
    double fault(0.);
    if (kE == kin && GetE() < 0.) {
      SetE(fault, false);
    } else if (kP == kin && GetP() < 0.) {
      SetP(fault, false);
    } else if (kPt == kin && GetPt() < 0.) {
      SetPt(fault, false);
    }
  }  
  
  void ParticleMCS::SetSmeared( bool flag) {kParticleSmeared=flag;}
  void ParticleMCS::SetESmeared( bool flag) {kESmeared=flag;}
  void ParticleMCS::SetPSmeared( bool flag) {kPSmeared=flag;}
  void ParticleMCS::SetPtSmeared( bool flag) {kPtSmeared=flag;}
  void ParticleMCS::SetPxSmeared( bool flag) {kPxSmeared=flag;}
  void ParticleMCS::SetPySmeared( bool flag) {kPySmeared=flag;}
  void ParticleMCS::SetPzSmeared( bool flag) {kPzSmeared=flag;}
  void ParticleMCS::SetThetaSmeared( bool flag) {kThetaSmeared=flag;}
  void ParticleMCS::SetPhiSmeared( bool flag) {kPhiSmeared=flag;}
  void ParticleMCS::SetIdSmeared( bool flag) {kIdSmeared=flag;}
  void ParticleMCS::SetNumSigmaSmeared( bool flag) {kNumSigmaSmeared=flag;}

}  // namespace Smear
