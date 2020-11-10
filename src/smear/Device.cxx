/**
 \file
 Implementation of class Smear::Device.
 
 \author    Michael Savastio
 \date      2011-08-19
 \copyright 2011 Brookhaven National Lab
 */

#include "eicsmear/smear/Device.h"

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


bool Device::Init( const TString& resolutionFunction, int genre) {
  Accept.SetGenre(genre);
  // Set the resolution function.
  mFormula = new FormulaString(resolutionFunction.Data());
  return mFormula;
}

Device::Device(KinType type, const TString& formula, EGenre genre)
: mSmeared(type)
, mFormula(NULL) {
  Accept.SetGenre(genre);
  Init(formula, genre);
}

Device::Device(const TString& variable, const TString& resolution,
               EGenre genre)
: mSmeared(kInvalidKinType)
, mFormula(NULL) {
  // Parse the "variable" string. It should have exactly one recognized KinType
  KinType kindummy;
  // need to make a copy since the string is destroyed
  TString vtmp = variable;
  int d = ParseInputFunction( vtmp, mSmeared, kindummy);
  if (d!=1) {
    std::cerr<< "Bad use of Device(const TString& variable, ...) constructor. Cannot parse variable" << endl;
    throw;
  }
    
  Init(resolution, genre);
}

Device::Device(const Device& that)
: Smearer(that)
, mSmeared(that.mSmeared)
, mFormula(NULL)
, mDimensions(that.mDimensions) {
  if (that.mFormula) {
    mFormula = static_cast<FormulaString*>(that.mFormula->Clone());
  }  // if
}

Device::~Device() {
  if (mFormula) {
    delete mFormula;
    mFormula = NULL;
  }  // if
}

void Device::Smear(const erhic::VirtualParticle &prt, ParticleMCS &out) {
  // Test for acceptance and do nothing if it fails.
  if (!Accept.Is(prt)) {
    return;
  }  // if
  // Get each argument for the resolution function from the particle.
  const std::vector<KinType> vars = mFormula->Variables();
  std::vector<double> args;
  for (unsigned i(0); i < vars.size(); ++i) {
    args.push_back(GetVariable(prt, vars.at(i)));
  }  // for
  // Evaluate the quantity to smear and the resolution, then throw
  // a random smeared value.
  double unsmeared = GetVariable(prt, mSmeared);
  double resolution = mFormula->Eval(args);
  double smeared = mDistribution.Generate(unsmeared, resolution);
  // mDistribution.Print();
  if ( false && abs(prt.Id())==11){
    // if ( mSmeared != kE ){std::cerr << " ======================  " << mSmeared << std::endl; throw(-1);}
    std::cout << " ======================  mSmeared = " << mSmeared << std::endl;
    // prt.Print();
    std::cout << " orig eta = " << prt.GetEta()<< std::endl;
    std::cout << "unsmeared " << unsmeared << std::endl;
    std::cout << "resolution " << resolution << std::endl;
    std::cout << "smeared " << smeared << std::endl;
    std::cout << "mSmeared " << mSmeared << std::endl;
  }

  out.SetVariable(smeared, mSmeared);
  // Fix angles to the correct ranges.
  if (kTheta == mSmeared) {
    out.SetTheta ( FixTheta(out.GetTheta() ), false);
  } else if (kPhi == mSmeared) {
    out.SetPhi ( FixPhi(out.GetPhi() ), false);
  }  // else if
  // Ensure E, p are positive definite
  out.HandleBogusValues(mSmeared);
  if ( std::isnan( GetVariable (out, mSmeared ) ) ) cout << "Problem. smeared is " << smeared << " but propagated nan in " << mSmeared<< endl;
  // std::cout << "Bye Smear" << std::endl << std::endl;
}

Device* Device::Clone(const char* /** Unused */) const {
  return new Device(*this);
}

void Device::Print(Option_t* /* option */) const {
  const std::string name = FormulaString::GetKinName(mSmeared);
  std::cout << "Device smearing " << name << " with sigma(" << name <<
  ") = " << mFormula->GetInputString() << std::endl;
}

}  // namespace Smear
