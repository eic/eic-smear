//
// Device.cxx
//
// Created by TB on 8/15/11.
// Copyright 2011 BNL. All rights reserved.
//

#include "Device.h"

#include <algorithm>
#include <functional>
#include <iterator>
#include <stdexcept>

#include <TRandom.h>
#include <TUUID.h>

namespace {
   // Valid patterns for the names of smearing function dimensions.
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
   // class ParamSimple
   //

   // static collection of valid patterns for parameterisations
   const std::vector<TString> ParamSimple::smPatterns(patterns, patterns + 6);

   // static lookup table of KinType keyed by text pattern
   std::map<TString, KinType> ParamSimple::smKinTypes;

   ParamSimple::ParamSimple(KinType type, const TString& formula, EGenre genre)
   : mSmeared(type)
   , mFormula(NULL) {
      Accept.SetGenre(genre);
      TString form(formula);
      std::vector<KinType> types;
      int nVariables = Parse(form, mDimensions);
      if(nVariables > 4 or nVariables < 0) {
         throw std::invalid_argument("ParamSimple limited to 4 dimesions");
      } // if
      // Use UUID for ROOT formula name to avoid instances clashing
      mFormula = new TFormula(TUUID().AsString(), form);
   }

   ParamSimple::ParamSimple(const ParamSimple& other)
   : Smearer(other)
   , mSmeared(other.mSmeared)
   , mFormula(NULL)
   , mDimensions(other.mDimensions) {
      mFormula = new TFormula(TUUID().AsString(), other.mFormula->GetTitle());
   }

   ParamSimple::~ParamSimple() {
      if(mFormula) {
         delete mFormula;
         mFormula = NULL;
      } // if
   }

   void ParamSimple::FillKinTypeTable() {
      if(smKinTypes.empty()) {
         smKinTypes.insert(std::make_pair("E", kE));
         smKinTypes.insert(std::make_pair("P", kP));
         smKinTypes.insert(std::make_pair("theta", kTheta));
         smKinTypes.insert(std::make_pair("phi", kPhi));
         smKinTypes.insert(std::make_pair("pT", kPt));
         smKinTypes.insert(std::make_pair("pZ", kPz));
      } // if
   }      

   int ParamSimple::Parse(TString& s, std::vector<KinType>& types) {
      using namespace std;
      FillKinTypeTable();
      // Get the patterns that appear in the input string
      list<TString> args;
      remove_copy_if(smPatterns.begin(), smPatterns.end(), back_inserter(args),
                     not1(bind1st(TStringContains(), s)));
      if(args.size() > 4) {
         return -1;
      } // if
      // Substitute each pattern (E, P, ...) with
      // variables recognised by ROOT::TFormula (x, y, z, t)
      list<TString> substitutions(tformulaVariables, tformulaVariables + 4);
      list<TString>::const_iterator iter = args.begin();
      for(; iter not_eq args.end() and not substitutions.empty(); ++iter) {
         s.ReplaceAll(*iter, substitutions.front());
         substitutions.pop_front();
         types.push_back(smKinTypes[*iter]);
      } // for
      return args.size();
   }

   void ParamSimple::DevSmear(const erhic::ParticleMC &prt,
                              ParticleS &prtOut) {
      if(Accept.Is(prt)) {
         // Evaluate each argument for the resolution function
         std::vector<double> args(4, 0.); // x, y, z, t to pass to TFormula
         for(unsigned i(0); i < mDimensions.size(); ++i) {
            args[i] = SwitchKinGetFromParticle(prt, mDimensions[i]);
         } // for
         Double_t unsmeared = SwitchKinGetFromParticle(prt, mSmeared);
         double resolution = mFormula->Eval(args[0], args[1], args[2], args[3]);
         double smeared = gRandom->Gaus(unsmeared, resolution);
         SwitchKinStoreToParticle(prtOut, smeared, mSmeared);
         if(kTheta == mSmeared) {
            prtOut.theta = FixTopologyTheta(prtOut.theta);
         } // if
         else if(kPhi == mSmeared) {
            prtOut.phi = FixTopologyPhi(prtOut.phi);
         } // else if
         // Ensure E, p are positive definite
         HandleBogusValues(prtOut, mSmeared);
      } //if
   }
   
   ParamSimple* ParamSimple::Clone(const char* /** Unused */) const {
      return new ParamSimple(*this);
   }

   //
   // class Device
   //

   Device::Device(KinType kin, TString parameterisation, int genre)
   : InKin1(kE)
   , InKin2(kTheta)
   , OutKin(kin)
   , ParDim(1)
   , ParMin1(0.)
   , ParMax1(1.e6)
   , ParMin2(0.)
   , ParMax2(1.e6)
   , name("Generic Device") {
      SetParametrization(parameterisation);
      SetGenre(genre);  
   }
   void Device::SetGenre(int n) {
      Accept.SetGenre(n);
   }
/*   
   void Device::SetGenre(const TString& genre) {
      Accept.SetGenre(genre);
   }
   */
   void Device::SetParDim(int n) {
      if(n==1 or n==2) {
         ParDim = n;	
      } // if
   }
   
   void Device::SetInOutKinematics(KinType In1, KinType In2, KinType Out) {
//      InKin1 = In1;
//      InKin2 = In2;
//      OutKin = Out;
//      ParDim = 2;
      SetInputKinematics(In1, In2);
      SetSmearedKinematics(Out);
   }
   
   void Device::SetInOutKinematics(KinType In, KinType Out) {
//      InKin1 = In;
//      OutKin = Out;
//      ParDim = 1;
      SetInputKinematics(In);
      SetSmearedKinematics(Out);
   }
   
   void Device::SetInputKinematics(KinType In1, KinType In2) {
      InKin1 = In1;
      InKin2 = In2;
      ParDim = 2;
   }
   
   void Device::SetInputKinematics(KinType In) { 
      InKin1 = In;
      ParDim = 1;
   }
   
   void Device::SetSmearedKinematics(KinType Out) {
      OutKin = Out;
   }
   
   void Device::SetDistribution(TString& s) {
      Distribution.SetDistribution(s);
   }
   
   void Device::SetDistribution(Distributor::Function f, int npars) {
      Distribution.SetDistribution(f,npars);
   }
   
   void Device::SetDistributionRange(double min, double max) {
      Distribution.SetRange(min,max);
   }
   
   void Device::SetMoveableRange(double plus, double minus) {
      Distribution.SetMoveableRange(plus,minus);
   }
   
   void Device::SetMoveable(bool b) {
      Distribution.SetMoveable(b);
   }
   
   void Device::SetDistributionBias(double bias) {
      Distribution.SetBias(bias);
   }
   
   KinType Device::GetTargetVariable() {
      return OutKin;
   }
   
   const TString& Device::GetDeviceType() {
      return name;
   }
   
   void Device::SetParametrization(TString par, bool parse) {
      if(parse) {
         ParDim = ParseInputFunction(par,InKin1,InKin2);
      } // if
      
      if(ParDim==2) {
         Param2 = TF2("f2",par,ParMin1,ParMin2,ParMax1,ParMax2);
      } // if
      else {
         Param1 = TF1("f1",par,ParMin1,ParMax1);
      } // else
   }
   
   void Device::SetParametrization(void *f, int npars) {
      if(ParDim==2) {
         Param2 = TF2("f2",f,ParMin1,ParMin2,ParMax1,ParMax2,npars);
      } // if
      else {
         Param1 = TF1("f1",f,ParMin1,ParMax1,npars);
      } // else
   }
   
   void Device::GetAcceptanceFromDevice(Device dev) {
      Accept = dev.Accept;
   }
   
   void Device::SetParRange(double min, double max) {
      Param2.SetRange(min,max);
      Param1.SetRange(min,max);
      ParMin1=min; ParMax1=max;
   }
   
   void Device::SetParRange(double min1, double min2,
                            double max1, double max2) {
      Param2.SetRange(min1,min2,max1,max2);
      Param1.SetRange(min1,max1);
      ParMin1=min1; ParMin2=min2; ParMax1=max1; ParMax2=max2;
   }
   
   void Device::SetParams(double xi0, double xi1, double xi2,
                          double xi3, double xi4, double xi5) {
      Param1.SetParameters(xi0,xi1,xi2,xi3,xi4,xi5);
      Param2.SetParameters(xi0,xi1,xi2,xi3,xi4,xi5);
   }
   
   void Device::SetParams(int n, double xi) {
      Param1.SetParameter(n,xi);
      Param2.SetParameter(n,xi);
   }
   
   TRandom3& Device::GetRandomGenerator() {
      return Distribution.Ran;
   }
   
   void Device::SetRanSeed(int n) {
      Distribution.SetSeed(n);
   }
   
   Device* Device::Clone(const char* /* Unused */) const {
      return new Device(*this);
   }
   
   double Device::EvaluateRes(double x) {
      return Param1.Eval(x);
   }
   
   double Device::EvaluateRes(double x, double y) {
      return Param2.Eval(x,y);
   }
   
   // TODO just rename to Smear
   void Device::DevSmear(const erhic::ParticleMC &prt, ParticleS &prtOut) {
//      bool isHadron = prt.KS==1 and abs(prt.id)>110 and Accept.Genre == 2;
//      if(isHadron)std::cout<<"Hadron accepted = "<<Accept.Is(prt)<<std::endl;
//      std::cout << "\t" << ClassName() << "::DevSmear" << std::endl;
      if(Accept.Is(prt)) {
         Double_t x1 = SwitchKinGetFromParticle(prt,InKin1);
         Double_t x2 = SwitchKinGetFromParticle(prt,InKin2);
         
         Double_t y = SwitchKinGetFromParticle(prt,OutKin);
//         if(isHadron)std::cout<<"x1="<<x1<<" x2="<<x2<<" y="<<y<<std::endl;
         if(ParDim==2) {
//            std::cout << InKin1 << " " << x1 << " " << x2 << " " << EvaluateRes(x1,x2) << std::endl;
            y = Distribution.Generate(y,EvaluateRes(x1,x2));
         } // if
         else {
//            std::cout << InKin1 << " " << x1 << " " << EvaluateRes(x1) << std::endl;
            y = Distribution.Generate(y,EvaluateRes(x1));
         } // else
//         std::cout << "Storing " << y << std::endl;
         SwitchKinStoreToParticle(prtOut,y,OutKin);
         
         //make sure angular coordinates live in S^2
         if(OutKin==kTheta) {
            prtOut.theta = FixTopologyTheta(prtOut.theta);
         } // if
         
         if(OutKin==kPhi) {
            prtOut.phi = FixTopologyPhi(prtOut.phi);
         } // if
         
         //make sure E, p are positive definite
         HandleBogusValues(prtOut,OutKin);
//         std::cout << "done" << std::endl;
      } //if
   }
} // namespace Smear
