//
// Device.cxx
//
// Created by TB on 8/15/11.
// Copyright 2011 BNL. All rights reserved.
//

#include "eicsmear/smear/Device.h"

#include <algorithm>
#include <functional>
#include <iterator>
#include <stdexcept>

#include <TUUID.h>

#include "eicsmear/smear/ParticleMCS.h"

namespace {
   // Valid patterns for the names of smearing function variables.
   // These are replaced with "x", "y", "z", "t" when the functions are parsed.
   TString patterns[6] = { "E", "P", "theta", "phi", "pT", "pZ" };

   // Variables recognised by ROOT::TFormula.
   // These replace the patterns above.
   TString tformulaVariables[4] = { "x", "y", "z", "t" };

   // Helper functor searching for a substring in a ROOT::TString.
   struct TStringContains
   : public std::binary_function<TString, TString, bool> {
      // Returns true if at least one instance of "pattern" is present in "str".
      // Returns false if not.
      bool operator()(const TString& str, const TString& pattern) const {
         return str.Contains(pattern);
      }
   };
} // namespace <anonymous>

namespace Smear {

   //
   // class Device
   //

   // static collection of valid patterns for parameterisations
   const std::vector<TString> Device::smPatterns(patterns, patterns + 6);

   // static lookup table of Smear::KinType keyed by text pattern
   std::map<TString, KinType> Device::smKinTypes;

   KinType Device::GetKinType(const TString& name) {
      FillKinTypeTable();
      KinType type(kInvalidKinType);
      if(smKinTypes.find(name) not_eq smKinTypes.end()) {
         type = smKinTypes[name];
      } // if
      return type;
   }

   TString Device::GetKinName(KinType type) {
      FillKinTypeTable();
      TString name;
      std::map<TString, KinType>::const_iterator iter;
      for(iter = smKinTypes.begin(); iter not_eq smKinTypes.end(); ++iter) {
         if(iter->second == type) {
            name = iter->first;
         } // if
      } // for
      return name;
   }

   bool Device::Init(const TString& kinematicFunction,
                     const TString& resolutionFunction, int genre) {
                     std::cout << "Initialising for " << kinematicFunction<<std::endl;
      Accept.SetGenre(genre);
      // Parse the kinematic function. This has to have exactly one variable.
      TString formula(kinematicFunction);
      std::vector<KinType> types;
      if(Parse(formula, types) not_eq 1) {
         return false;
      } // if
      mSmeared = types.front();
      // Use UUID for ROOT function name to avoid instances clashing
      mKinematicFunction = new TF1(TUUID().AsString(), formula, 0., 1.e16);
      // Parse the resolution function.
      // This can have between 0 and 4 variables.
      formula = resolutionFunction;
      int nVariables = Parse(formula, mDimensions);
      if(nVariables > 4 or nVariables < 0) {
         return false;
      } // if
      // Use UUID for ROOT formula name to avoid instances clashing
      mFormula = new TFormula(TUUID().AsString(), formula);
      return mFormula and mKinematicFunction;
   }
   
   Device::Device(KinType type, const TString& formula, EGenre genre)
   : mSmeared(type)
   , mKinematicFunction(NULL)
   , mFormula(NULL) {
      Accept.SetGenre(genre);
      FillKinTypeTable();
      Init(GetKinName(type), formula, genre);
   }
   Device::Device(const TString& variable, const TString& resolution, EGenre genre)
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
      //      mFormula = new TFormula(TUUID().AsString(), that.mFormula->GetTitle());
      if(that.mKinematicFunction) {
         mKinematicFunction = dynamic_cast<TF1*>(
            that.mKinematicFunction->Clone(TUUID().AsString()));
      } // if
      if(that.mFormula) {
         mFormula = dynamic_cast<TFormula*>(
            that.mFormula->Clone(TUUID().AsString()));
      } // if
      FillKinTypeTable();
   }

   Device::~Device() {
      if(mFormula) {
         delete mFormula;
         mFormula = NULL;
      } // if
      if(mKinematicFunction) {
         delete mKinematicFunction;
         mKinematicFunction = NULL;
      } // if
   }

   void Device::FillKinTypeTable() {
      if(smKinTypes.empty()) {
         smKinTypes.insert(std::make_pair("E", kE));
         smKinTypes.insert(std::make_pair("P", kP));
         smKinTypes.insert(std::make_pair("theta", kTheta));
         smKinTypes.insert(std::make_pair("phi", kPhi));
         smKinTypes.insert(std::make_pair("pT", kPt));
         smKinTypes.insert(std::make_pair("pZ", kPz));
      } // if
   }

   int Device::Parse(TString& s, std::vector<KinType>& types) {
      typedef std::list<TString> SList;
      // Get all the patterns that appear in the input string
      SList args;
      std::remove_copy_if(smPatterns.begin(),
                          smPatterns.end(),
                          std::back_inserter(args),
                          std::not1(std::bind1st(TStringContains(), s)));
      // TFormula supports at most four variables.
      // Zero variables (constant formula) is OK though.
      if(args.size() > 4) {
         return -1;
      } // if
      // Substitute each pattern (E, P etc) with
      // variables recognised by ROOT::TFormula (x, y, z, t)
      // Also, accumulate the Smear::KinType corresponding to each pattern
      // in the vector function argument.
      SList substitutions(tformulaVariables,
                          tformulaVariables + 4);
      
      for(SList::const_iterator iter = args.begin();
          iter not_eq args.end() and not substitutions.empty();
          ++iter) {
         s.ReplaceAll(*iter, substitutions.front());
         substitutions.pop_front();
         types.push_back(smKinTypes[*iter]);
      } // for
      return args.size();
   }

   void Device::Smear(const erhic::VirtualParticle &prt,
                           ParticleMCS &prtOut) {
      if(not Accept.Is(prt)) {
         return;
      } // if
      // Evaluate each argument for the resolution function.
      // Get the x, y, z, t to pass to TFormula.
      std::vector<double> args(4, 0.);
      for(unsigned i(0); i < mDimensions.size(); ++i) {
         args.at(i) = SwitchKinGetFromParticle(prt, mDimensions.at(i));
      } // for
      double unsmeared = mKinematicFunction->Eval(
                            SwitchKinGetFromParticle(prt, mSmeared));
      double resolution = mFormula->Eval(args[0], args[1], args[2], args[3]);
      double smeared = mDistribution.Generate(unsmeared, resolution);
      SwitchKinStoreToParticle(prtOut, mKinematicFunction->GetX(smeared),
                                  mSmeared);
      if(kTheta == mSmeared) {
         prtOut.theta = FixTopologyTheta(prtOut.theta);
      } // if
      else if(kPhi == mSmeared) {
         prtOut.phi = FixTopologyPhi(prtOut.phi);
      } // else if
      // Ensure E, p are positive definite
      HandleBogusValues(prtOut, mSmeared);
   }
   
   Device* Device::Clone(const char* /** Unused */) const {
      return new Device(*this);
   }
} // namespace Smear
