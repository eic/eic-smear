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


  /// TODO: KK It seems the use of kinematicFunction is unnecessary
  /// All it does is
  // - create a string that identifies the smeared observable (say, phi)
  // - create a string ("phi"), translate it to a TF1 ("x")
  // - relate the smeared value to the observable using GetX()
  //But this TF1 is always identity, modulo range restrictions that cause obscure errors
bool Device::Init(const TString& kinematicFunction,
                  const TString& resolutionFunction, int genre) {
  Accept.SetGenre(genre);
  // Use FormulaString to parse the kinematic function,
  // though we can't use a FormulaString for the kinematic
  // function as we need TF1 functionality.
  FormulaString f(kinematicFunction.Data());
  // The expression has to have exactly one variable.
  mSmeared = f.Variables().front();
  // Use UUID for ROOT function name to avoid instances clashing
  mKinematicFunction = new TF1(TUUID().AsString(),
                               f.GetString().c_str(), -1e15, 1.e16);
  // Set the resolution function.
  mFormula = new FormulaString(resolutionFunction.Data());
  // cout << mFormula->GetString() << endl;
  return mFormula && mKinematicFunction;
}

Device::Device(KinType type, const TString& formula, EGenre genre)
: mSmeared(type)
, mKinematicFunction(NULL)
, mFormula(NULL) {
  Accept.SetGenre(genre);
  Init(FormulaString::GetKinName(type), formula, genre);
}

Device::Device(const TString& variable, const TString& resolution,
               EGenre genre)
: mSmeared(kInvalidKinType)
, mKinematicFunction(NULL)
, mFormula(NULL) {
  Init(variable, resolution, genre);
}

Device::Device(const Device& that)
: Smearer(that)
, mSmeared(that.mSmeared)
, mKinematicFunction(NULL)
, mFormula(NULL)
, mDimensions(that.mDimensions) {
  if (that.mKinematicFunction) {
    mKinematicFunction = static_cast<TF1*>(
        that.mKinematicFunction->Clone(TUUID().AsString()));
  }  // if
  if (that.mFormula) {
    mFormula = static_cast<FormulaString*>(that.mFormula->Clone());
  }  // if
}

Device::~Device() {
  if (mFormula) {
    delete mFormula;
    mFormula = NULL;
  }  // if
  if (mKinematicFunction) {
    delete mKinematicFunction;
    mKinematicFunction = NULL;
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
  double unsmeared = mKinematicFunction->Eval(GetVariable(prt, mSmeared));
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
  // mKinematicFunction->Print();
  if ( mKinematicFunction->GetFormula()->GetExpFormula() != "x" ){
    std::cerr << "Formula = " << mKinematicFunction->GetFormula()->GetExpFormula() << std::endl;
    throw -1;
  }
  // if ( smeared < 0 || fabs(mKinematicFunction->GetX(smeared)-smeared)>1e-5 ) cout << mKinematicFunction->GetX(smeared) << " " << smeared << endl;
  out.SetVariable(mKinematicFunction->GetX(smeared), mSmeared);
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
