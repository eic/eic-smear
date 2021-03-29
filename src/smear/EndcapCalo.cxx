/**
 \file
 Implementation of class Smear::EndcapCalo.
 
 \author    Kolja Kauder
 \date      2020-11-05
 \copyright 2020 Brookhaven National Lab
 */

#include "eicsmear/smear/EndcapCalo.h"

#include <algorithm>
#include <functional>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

#include <TUUID.h>
#include <TDatabasePDG.h>

#include "eicsmear/smear/FormulaString.h"
#include "eicsmear/smear/ParticleMCS.h"

using std::cout;
using std::cerr;
using std::endl;

namespace Smear {

EndcapCalo::EndcapCalo(const TString SpatialFormula, const EGenre genre,
		       const double Zposition,
		       const bool Projective, const double ProjPhiAngle, const double ProjThetaAngle,
		       const TString EResFormula)
  : mZposition(Zposition)
  , mProjective(Projective)
  , mProjPhiAngle(ProjPhiAngle)
  , mProjThetaAngle(ProjThetaAngle)
{  
  mSmeared = kE;
  mFormula=0;
  if (EResFormula !=""){
    mFormula = new FormulaString(EResFormula.Data());
  }
  mSpatialFormula = new FormulaString(SpatialFormula.Data());
  
  Accept.SetGenre( genre );
}

EndcapCalo::EndcapCalo(const EndcapCalo& that)
: Device(that)
, mZposition(that.mZposition)
, mProjective(that.mProjective)
, mProjPhiAngle(that.mProjPhiAngle)
, mProjThetaAngle(that.mProjThetaAngle)
{
  if (that.mSpatialFormula) {
    mSpatialFormula = static_cast<FormulaString*>(that.mSpatialFormula->Clone());
  }  // if
}

  

EndcapCalo::~EndcapCalo() {
  if (mFormula) {
    delete mFormula;
    mFormula = NULL;
  }  // if
  if (mSpatialFormula) {
    delete mSpatialFormula;
    mSpatialFormula = NULL;
  }  // if
}

void EndcapCalo::Smear(const erhic::VirtualParticle &prt, ParticleMCS &out) {
  // Test for acceptance and do nothing if it fails.
  if (!Accept.Is(prt)) {
    return;
  }  // if

  // if we have an E-Smearing, apply it now via the baseclass
  if (mFormula) {
    Device::Smear( prt, out);
  }

  // Get the spatial resolution (reusing code from the base class)
  const std::vector<KinType> vars = mSpatialFormula->Variables();
  std::vector<double> args;
  for ( unsigned int i=0; i < vars.size(); ++i) {
    args.push_back(GetVariable(prt, vars.at(i)));
  }  // for
  double spresolution = mSpatialFormula->Eval(args);
  double spdeltaphi   = mDistribution.Generate(0, spresolution);
  double spdeltatheta = mDistribution.Generate(0, spresolution);

  // // adapt for (imperfect) projectivity
  // if ( mProjective ){
  //   if ( mProjPhiAngle!=0 ) spdeltaphi = std::sqrt ( pow( spdeltaphi,2 ) + pow ( mRadLength*std::sin( mProjPhiAngle ),2) );
  //   if ( mProjThetaAngle!=0 ) spdeltatheta = std::sqrt ( pow( spdeltatheta,2 ) + pow ( mRadLength*std::sin( mProjThetaAngle ),2) );
  // }

  // truth values
  auto phitruth = GetVariable(prt, kPhi);
  auto thetatruth = GetVariable(prt, kTheta);

  // translate to angles
  double deltaphi = spdeltaphi / (mZposition * std::tan (thetatruth ));
  double deltatheta = pow ( std::cos(thetatruth),2)  * spdeltatheta  / mZposition;

  out.SetVariable(phitruth + deltaphi, kPhi);
  out.SetPhi ( FixPhi(out.GetPhi() ), false);
  
  out.SetVariable(thetatruth + deltatheta, kTheta);
  out.SetTheta ( FixTheta(out.GetTheta() ), false);

}

EndcapCalo* EndcapCalo::Clone(const char* /** Unused */) const {
  return new EndcapCalo(*this);
}

// void EndcapCalo::Print(Option_t* /* option */) const {
//   const std::string name = FormulaString::GetKinName(mSmeared);
//   std::cout << "EndcapCalo smearing " << name << " with sigma(" << name <<
//   ") = " << mFormula->GetInputString() << std::endl;
// }

}  // namespace Smear
