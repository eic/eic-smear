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
    if ( ThePidObject->valid(eta,p) ){
      prtOut.SetNumSigmaElectron ( ThePidObject->numSigma(eta, p, pdgtruth, PID::kElectron) );
      prtOut.SetNumSigmaPion     ( ThePidObject->numSigma(eta, p, pdgtruth, PID::kPion), false ); // only check/set the smearing flag the first time
      prtOut.SetNumSigmaProton   ( ThePidObject->numSigma(eta, p, pdgtruth, PID::kProton), false );
      prtOut.SetNumSigmaKaon     ( ThePidObject->numSigma(eta, p, pdgtruth, PID::kKaon), false );
      prtOut.SetNumSigmaMuon     ( ThePidObject->numSigma(eta, p, pdgtruth, PID::kMuon), false );
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
