/**
 \file
 Implementation of class Smear::BarrelCalo.
 
 \author    Kolja Kauder
 \date      2020-11-05
 \copyright 2020 Brookhaven National Lab
 */

#include "eicsmear/smear/BarrelCalo.h"

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

BarrelCalo::BarrelCalo(const TString SpatialFormula, const EGenre genre,
		       const double InnerRadius, const double OuterRadius,	   
		       const bool Projective, const double ProjPhiAngle, const double ProjThetaAngle,
		       const TString EResFormula)
  : mInnerRadius(InnerRadius)
  , mOuterRadius(OuterRadius)
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
  // cout << SpatialFormula.Data() << endl;
  cout << mSpatialFormula->GetString() << endl;
  
  Accept.SetGenre( genre );
}

BarrelCalo::BarrelCalo(const BarrelCalo& that)
: Device(that)
, mInnerRadius(that.mInnerRadius)
, mOuterRadius(that.mOuterRadius)
, mProjective(that.mProjective)
, mProjPhiAngle(that.mProjPhiAngle)
, mProjThetaAngle(that.mProjThetaAngle)
{
  if (that.mSpatialFormula) {
    mSpatialFormula = static_cast<FormulaString*>(that.mSpatialFormula->Clone());
  }  // if
}

  

BarrelCalo::~BarrelCalo() {
  if (mFormula) {
    delete mFormula;
    mFormula = NULL;
  }  // if
  if (mSpatialFormula) {
    delete mSpatialFormula;
    mSpatialFormula = NULL;
  }  // if
}

void BarrelCalo::Smear(const erhic::VirtualParticle &prt, ParticleMCS &out) {
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

  // adapt for (imperfect) projectivity
  if ( mProjective ){
    if ( mProjPhiAngle!=0 ) spdeltaphi = std::sqrt ( pow( spdeltaphi,2 ) + pow ( mRadLength*std::sin( mProjPhiAngle ),2) );
    if ( mProjThetaAngle!=0 ) spdeltatheta = std::sqrt ( pow( spdeltatheta,2 ) + pow ( mRadLength*std::sin( mProjThetaAngle ),2) );
  }

  // truth values
  auto phitruth = GetVariable(prt, kPhi);
  auto thetatruth = GetVariable(prt, kTheta);

  // // translate to angles
  double deltaphi = spdeltaphi / mInnerRadius;
  double deltatheta = pow ( std::sin(thetatruth),2)  * spdeltatheta  / mInnerRadius;
  if ( mProjective ){
    deltatheta = std::sin(thetatruth)  * spdeltatheta  / mInnerRadius;
  }

  out.SetVariable(phitruth + deltaphi, kPhi);
  out.SetPhi ( FixPhi(out.GetPhi() ), false);
  
  out.SetVariable(thetatruth + deltatheta, kTheta);
  out.SetTheta ( FixTheta(out.GetTheta() ), false);

}

BarrelCalo* BarrelCalo::Clone(const char* /** Unused */) const {
  return new BarrelCalo(*this);
}

// void BarrelCalo::Print(Option_t* /* option */) const {
//   const std::string name = FormulaString::GetKinName(mSmeared);
//   std::cout << "BarrelCalo smearing " << name << " with sigma(" << name <<
//   ") = " << mFormula->GetInputString() << std::endl;
// }

}  // namespace Smear
