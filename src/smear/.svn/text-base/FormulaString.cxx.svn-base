/**
 \file
 Implementation of class Smear::FormulaString.
 
 \author    Michael Savastio
 \date      2011-08-19
 \copyright 2011 Brookhaven National Lab
 */

#include "eicsmear/smear/FormulaString.h"

#include <algorithm>
#include <iterator>
#include <list>
#include <map>
#include <string>
#include <vector>

#include <TFormula.h>
#include <TString.h>
#include <TUUID.h>

namespace {

// Valid patterns for the names of smearing function variables.
// These are replaced with "x", "y", "z", "t" when the functions are parsed.
std::vector<std::string> functionPatterns;

// Returns a reference to the function pattern vector, filling
// it first if needed.
const std::vector<std::string>& patterns() {
  if (functionPatterns.empty()) {
    const std::string array[6] = {"E", "P", "theta", "phi", "pT", "pZ"};
    functionPatterns.assign(array, array + 6);
  }  // if
  return functionPatterns;
}

// Variables recognised by ROOT::TFormula.
// These replace the functionPatterns above so ROOT::TFormula can
// understand the formula string.
const std::string tformulaVariables[4] = {"x", "y", "z", "t"};

// KinType keyed by corresponding string name.
std::map<std::string, Smear::KinType> kinTypes;

// Returns a reference to the KinType table, populating it first if needed.
const std::map<std::string, Smear::KinType>& kinTypeTable() {
  if (kinTypes.empty()) {
    kinTypes.insert(std::make_pair("E", Smear::kE));
    kinTypes.insert(std::make_pair("P", Smear::kP));
    kinTypes.insert(std::make_pair("theta", Smear::kTheta));
    kinTypes.insert(std::make_pair("phi", Smear::kPhi));
    kinTypes.insert(std::make_pair("pT", Smear::kPt));
    kinTypes.insert(std::make_pair("pZ", Smear::kPz));
  }  // if
  return kinTypes;
}

}  // anonymous namespace

namespace Smear {

// This is basically a wrapper around ROOT::TFormula, with functionality
// to parse functions containing "P", "theta" etc into "x, y, z, t"
// as accepted by TFormula.

FormulaString::~FormulaString() {
  if (mFormula) {
    delete mFormula;
    mFormula = NULL;
  }  // if
}

FormulaString::FormulaString()
: mFormula(NULL) {
}

FormulaString::FormulaString(const std::string& formula)
: mFormula(NULL)
, mInput(formula) {
  std::string f = Parse(formula);
  mFormula = new TFormula(TUUID().AsString(), f.c_str());
}

double FormulaString::Eval(const std::vector<double>& args) const {
  if (args.size() != mVariables.size()) {
    std::cerr << "FormulaString::Eval() got " << args.size() <<
    " arguments, expected " << mVariables.size() << std::endl;
  }  // if
  // TFormula accepts up to four arguments.
  // Use default zeros for absent arguments.
  std::vector<double> a(args);
  a.resize(4, 0.);
  return mFormula->Eval(a.at(0), a.at(1), a.at(2), a.at(3));
}

std::vector<KinType> FormulaString::Variables() const {
  return mVariables;
}

std::string FormulaString::Parse(const std::string& formula) {
  using std::string;
  using std::vector;
  mVariables.clear();
  // Get all the patterns that appear in the input string
  // We want to preserve the order in which the arguments appear
  // so store the position of each argument with its name.
  std::map<int, string> args;
  typedef vector<string>::const_iterator StrIter;
  for (StrIter i = patterns().begin(); i != patterns().end(); ++i) {
    size_t position = formula.find(*i);
    if (position != string::npos) {
      args.insert(std::make_pair(position, *i));
    }  // if
  }  // for
  // Substitute each pattern (E, P etc) with
  // variables recognised by ROOT::TFormula (x, y, z, t)
  // Also, accumulate the Smear::KinType corresponding to each pattern
  // in the vector function argument.
  std::list<string> substitutions(tformulaVariables,
                                  tformulaVariables + 4);
  TString s(formula);
  for (std::map<int, string>::const_iterator i = args.begin();
      i != args.end(); ++i) {
    s.ReplaceAll(i->second, substitutions.front());
    substitutions.pop_front();
    mVariables.push_back(GetKinType(i->second));
  }  // for
  return string(s.Data());
}

std::string FormulaString::GetString() const {
  std::string str;
  if (mFormula) {
    str = mFormula->GetTitle();
  }  // if
  return str;
}

std::string FormulaString::GetInputString() const {
  return mInput;
}

KinType FormulaString::GetKinType(const std::string& name) {
  KinType type(kInvalidKinType);
  if (kinTypeTable().find(name) != kinTypeTable().end()) {
    type = kinTypes[name];
  }  // if
  return type;
}

std::string FormulaString::GetKinName(KinType type) {
  std::string name;
  std::map<std::string, KinType>::const_iterator i;
  for (i = kinTypeTable().begin(); i != kinTypeTable().end(); ++i) {
    if (i->second == type) {
      name = i->first;
    }  // if
  }  // for
  return name;
}

}  // namespace Smear
