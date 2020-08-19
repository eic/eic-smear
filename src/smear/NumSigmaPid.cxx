#include "eicsmear/smear/NumSigmaPid.h"

namespace Smear {

  // -----------------------------------------------------------
  NumSigmaPid::NumSigmaPid(){}
  
  // -----------------------------------------------------------
  void NumSigmaPid::Smear(const erhic::VirtualParticle& prt,
			  ParticleMCS& prtOut) {
    auto p = prt.GetP();
    auto eta = prt.GetEta();
    // std::cout << "prt is at eta = " << eta << " p= " << p << std::endl;
    if ( ThePidObject->valid(eta,p) ){
      // std::cout << "--> prt is valid at eta = " << eta << " p= " << p << std::endl;
      prtOut.SetNumSigma ( ThePidObject->numSigma(eta, p, EnumType) );
      // std::cout << EnumType << "  " << ThePidObject->numSigma(eta, p, EnumType) << "  " << prtOut.GetNumSigma() << std::endl;
    }
  }
  
  // -----------------------------------------------------------
  NumSigmaPid* NumSigmaPid::Clone(const char*) const {
    // TODO: Probably should add a proper copy ctor to hand over type etc.
    return new NumSigmaPid(*this);
  }

  // -----------------------------------------------------------
  void NumSigmaPid::SetNumSigmaType( const int i ){    
    switch ( i ){
    case 0 :
      NumSigmaType = 0;
      EnumType = PID::pi_k;
      break;
    case 1 :
      NumSigmaType = 1;
      EnumType = PID::k_p;
      break;
    default :
      std::cerr << "Unrecognized NumSigmaType " << i  << std::endl;
      throw;
      break;
    }
    
  }
  
  // -----------------------------------------------------------
  int NumSigmaPid::GetNumSigmaType( ) const {
    return NumSigmaType;
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
