#include "eicsmear/smear/NumSigmaPid.h"

namespace Smear {

  // -----------------------------------------------------------
  NumSigmaPid::NumSigmaPid(){}
  
  // -----------------------------------------------------------
  void NumSigmaPid::Smear(const erhic::VirtualParticle& prt,
			  ParticleMCS& prtOut) {
    auto p = prt.GetP();
    auto eta = prt.GetEta();
    auto pdgtruth = prt.Id().Code();
    PID::Species pidtruth = PID::kElectron;
    switch ( std::abs(pdgtruth) ){
    case 11:   pidtruth = PID::kElectron; break;
    case 211:  pidtruth = PID::kPion;     break;
    case 2212: pidtruth = PID::kProton;   break;
    case 321:  pidtruth = PID::kKaon;     break;
    case 13:   pidtruth = PID::kMuon;     break;
    default : // not a recognized species. Do nohing.
      // Warn?
      std::cerr << "Warning: Not a recognized PID species " << pdgtruth << std::endl;
      return;
    }

    if ( ThePidObject->valid(eta,p) ){
      prtOut.SetNumSigmaElectron ( ThePidObject->numSigma(eta, p, pidtruth, PID::kElectron), false );
      prtOut.SetNumSigmaPion     ( ThePidObject->numSigma(eta, p, pidtruth, PID::kPion), false ); // only check/set the smearing flag the first time
      prtOut.SetNumSigmaProton   ( ThePidObject->numSigma(eta, p, pidtruth, PID::kProton), false );
      prtOut.SetNumSigmaKaon     ( ThePidObject->numSigma(eta, p, pidtruth, PID::kKaon), false );
      prtOut.SetNumSigmaMuon     ( ThePidObject->numSigma(eta, p, pidtruth, PID::kMuon), false );      
      
    }
  }
  
  // -----------------------------------------------------------
  NumSigmaPid* NumSigmaPid::Clone(const char*) const {
    // TODO: Probably should add a proper copy ctor to hand over type etc.
    return new NumSigmaPid(*this);
  }

  // -----------------------------------------------------------
  bool NumSigmaPid::valid (double eta, double p ) const {
    return ThePidObject->valid( eta, p);
  }

  double NumSigmaPid::maxP (double eta, double numSigma, PID::type PID) const{
    return ThePidObject->maxP (eta, numSigma, PID );
  }

  double NumSigmaPid::minP (double eta, double numSigma, PID::type PID) const{
    return ThePidObject->minP (eta, numSigma, PID );
  }

  std::string NumSigmaPid::name() const {
    return ThePidObject->name();
  }

  void NumSigmaPid::description() const{
    ThePidObject->description();
  }
  // -----------------------------------------------------------

  
}
