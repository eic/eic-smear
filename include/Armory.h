#ifndef _ERHIC_BUILDTREE_ARMORY_
#define _ERHIC_BUILDTREE_ARMORY_

// This file contains a list of specialized Device classes.

#include <TString.h>
#include <TRandom3.h>
#include <TF1.h>
#include <fstream>
#include <TF2.h>
#include <TLorentzVector.h>
#include <iostream>
#include <sstream>
#include <cmath>

#include "VirtualEvent.h"
#include "VirtualParticle.h"
#include "VirtualEvent.h"
#include "Kinematics.h"
#include "ParticleIdentifier.h"

#include "Smear.h"
#include "Acceptance.h"

#include "Device.h"

namespace Smear {
	
  /**
	 Specialized Device class for simulating tracking detectors.  Contains a special parametrization suited for tracking
	 detectors and methods for setting the various input parameters (such as magnetic field and track length).
	 If you want to simulate a tracker by entering a parametrization direclty, you will have to use the genereic Device class.
	 */
  struct Tracker: public Smear::Device {
	
		Tracker() {
			B = 2.;
			Nrl = 0.03;
			SigmaRPhi = 0.001;
			n = 25.;
			R1=0.1;
			R2=1.;
			l=6.;
			fault = -999.;
			
			bMoreCrit = false;
			
			P = NULL;
			
			SetInputKinematics(kP);
			SetSmearedKinematics(kP);
			
			FixAcceptance();
			
			name = "Tracker";
		}
		
		Tracker(double inner, double outer, double L, double Bb, double nrl, double sigmaRPhi, double N) {
			Tracker();
			SetVals(inner,outer,L,Bb,nrl,sigmaRPhi,N);
		}
		
		Tracker(double inner, double outer, double L) {
			Tracker();
			SetDimensions(inner,outer,L);
		}
		
		double B;
		double Nrl;
		double SigmaRPhi;
		double n;
		double R1;
		double R2;
		double l;
		double fault;
		
		bool bMoreCrit;
		
		const Particle *P;
		
		double MultipleScatteringContribution() {
			if (!P) return fault;
			double c1=0.0136/0.3;
			return c1*P->p*sqrt(Nrl)/(L()*beta()*B*cosSquaredgamma());
		}
		
		double IntrinsicContribution() {
			if (!P) return fault;
			double c1=sqrt(720.)/0.3;
			double over = c1*P->p*P->p*SigmaRPhi;
			return over/(B*pow(LPrime(),2)*sqrt(n+4.));
		}
		
		double EvaluateRes() {
			if (P) {
				//std::cout << "Mult, Intrins " << MultipleScatteringContribution() << " " << IntrinsicContribution() << std::endl;
				//std::cout << ">>>>>>>>>>>>>>>>> " << P->p << " " << P->theta << " " << std::endl;
				//std::cout << "     ------  -- " << L() << " " << LPrime() << " " << beta() << " " << cosSquaredgamma() << std::endl;
				//std::cout << "RES " << MultipleScatteringContribution()+IntrinsicContribution() << std::endl;
			  return MultipleScatteringContribution()+IntrinsicContribution();
			} else {
				return fault;
			}
		}
		
		/**
		 Set the magnetic field in the tracker in teslas.
		 */
		void SetMagneticField(double Bb) {
			B = Bb;
		}
		
		/**
		 Set the inner and outer radius of your cylindrical tracker in meters.  Automatically calls FixAcceptance,
		 so if you want a custom acceptance you should call it AFTER calling this.
		 */
		void SetRadii(double inner, double outer) {
			if (inner < outer) {
			  R1 = inner;  R2 = outer;
			  FixAcceptance();
			} else {
				std::cerr << "ERROR! Outer radius must exceed inner.\n";
				return;
			}
		}
		
		/**
		 Set the total length of your tracker (automatically centered at IP) in meters.  Automatically calls FixAcceptance,
		 so if you want a custom acceptance you should call it AFTER calling this.
		 */
		void SetLength(double L) {
			l = L;
			FixAcceptance();
		}
		
		/**
		 Set the radii and length in meters.  Automatically calls FixAcceptance, so if you want a custom acceptance, you should
		 set it AFTER calling this.
		 */
		void SetDimensions(double inner, double outer, double L) {
			SetRadii(inner,outer);
			SetLength(L);
		}
		
		double cosSquaredgamma() {
			if (!P) return fault;
			return cos(pi/2. - P->theta)*cos(pi/2. - P->theta);
		}
		
		double L() {
			if (!P) return fault;
			if (bMoreCrit) {
				return LPrime()/sin(P->theta);
			} else {
				return sqrt(pow(LPrime(),2)+pow((l/2.)-(R1/tanTheta()),2));
			}
		}
		
		double LPrime() {
			if (!P) return fault; 
			if (bMoreCrit) {
				return R2-R1;
			} else {
				return (l/2.)*tanTheta()-R1;
			}
		}
		
		double tanTheta() {
			if (!P) return fault;
			return fabs(tan(P->theta));
		}
		
		double ThetaCrit() {
			return atan(2.*R2/l);
		}
		
		/**
		 Set the number of radiation lengths in the tracker.
		 */
		void SetNumberOfRadiationLengths(double nrl) {
			Nrl = nrl;
		}
		
		double beta() {
			if (!P) return fault;
			return (P->p)/(P->E);
		}
		
		/**
		 Set the position resolution in meters.
		 */
		void SetPositionResolution(double sigmaRPhi) {
			SigmaRPhi = sigmaRPhi;
		}
		
		/**
		 Set the number of measurements.
		 */
		void SetNumberOfMeasurements(double N) {
			n = N;
		}
											
		/**
		 Returns the minimum theta of particle accepted by a tracker with the length and radii you specified.
		 If you must set a different acceptance in theta, it is recomended that you allow it no closer to the beam
		 line than this.
		 */
		double GetThetaMin() {
			return atan(2.*R1/l);
		}
		
		/**
		 Automatically set acceptance to that of a cylinder with the dimensions you set.
		 */
		void FixAcceptance() {
			Accept.SetTheta(GetThetaMin(),pi-GetThetaMin());									
		}
		
		//note: every inherited device should contain this function to work properly with Detector
		Tracker *Clone() {
			Tracker *dev = new Tracker();
			*dev = *this;
			return dev;
		}
		
		void SetParticle(const Particle &prt) {
			P = &prt;
			if (P->theta >= ThetaCrit() && P->theta < pi-ThetaCrit()) {
				bMoreCrit = true;
			} else {
				bMoreCrit = false;
			}
		}
		
		void SetVals(double inner, double outer, double L, double Bb, double nrl, double sigmaRPhi, double N) {
			SetDimensions(inner,outer,L);  SetMagneticField(Bb);  SetNumberOfRadiationLengths(nrl);
			SetPositionResolution(sigmaRPhi);  SetNumberOfMeasurements(N);
		}
		
		void DevSmear(const Particle &prt, ParticleS &prtOut) {
         
			Double32_t x1; Double32_t x2;
			Double32_t y;
			
			if (Accept.Is(prt)) {
				
				SetParticle(prt);
				
				x1 = SwitchKinGetFromParticle(prt,InKin1);
				x2 = SwitchKinGetFromParticle(prt,InKin2);
				
				y = SwitchKinGetFromParticle(prt,OutKin);
				
				y = Distribution.Generate(y,EvaluateRes());
				
				SwitchKinStoreToParticle(prtOut,y,OutKin);
				
				//make sure angular coordinates live in S^2
				if (OutKin==kTheta) prtOut.theta = FixTopologyTheta(prtOut.theta);
				if (OutKin==kPhi) prtOut.phi = FixTopologyPhi(prtOut.phi);
				
				//make sure E, p are positive definite
				HandleBogusValues(prtOut,OutKin);
				
			} //if
		}
	
		ClassDef(Tracker,1)
	};
	
	/**
	 Specialized EM calorimeter device.
	 */
	struct EMCalorimeter: public Smear::Device {
		
		EMCalorimeter() {
		  c1 = 0.3;
			c2 = 0.;
			
			Accept.SetGenre(1);
			SetInputKinematics(kE);
			SetSmearedKinematics(kE);
			
			name = "EM Calorimeter";
		}
		
		EMCalorimeter(double C1, double C2=0.) {
			EMCalorimeter();
			SetVals(C1,C2);
		}
		
		double c1;
		double c2;
		
		double EvaluateRes(double E) {
			return c1*sqrt(E)+c2*E;
		}
		
		EMCalorimeter *Clone() {
			EMCalorimeter *dev = new EMCalorimeter();
			*dev = *this;
			return dev;
		}
		
		void SetCoefficient(double C1, double C2=0.) {
			c1 = C1;  c2 = C2;
		}
		
		void SetVals(double C1, double C2) {
			SetCoefficient(C1,C2);
		}
		
		ClassDef(EMCalorimeter,1)
	};
	
	
	/**
	 Specialized hadronic calorimeter device.
	 */
	struct HCalorimeter: public Smear::Device {
		
		HCalorimeter() {
			c1 = 0.3;
			c2 = 0.;
			
			Accept.SetGenre(2);
			SetInputKinematics(kE);
			SetSmearedKinematics(kE);
			
			name = "Hadronic Calorimeter";
		}
		
		HCalorimeter(double C1, double C2=0.) {
			HCalorimeter();
			SetVals(C1,C2);
		}
		
		double c1;
		double c2;
		
		double EvaluateRes(double E) {
			return c1*sqrt(E)+c2*E;
		}
		
		HCalorimeter *Clone() {
			HCalorimeter *dev = new HCalorimeter();
			*dev = *this;
			return dev;
		}
		
		void SetCoefficient(double C1, double C2=0.) {
			c1 = C1;  c2 = C2;
		}
		
		void SetVals(double C1, double C2) {
			SetCoefficient(C1,C2);
		}
		
		ClassDef(HCalorimeter,1)
	};
	
	
	/**
	 Device class for smearing arbitrary functions of input kinematics.  By default, if set to smear kinematic variable X,
	 this device will smear 1/X and then store the smeared value back into Particle.X.  If provided with some function via
	 the constructor, or with SetInputFunction, if the input function of f(X) (you can specify X using E,P,theta etc as with
	 parametrizations), then the Devious will smear to y=f(x) and store x=f^-1(y) in Particle.X.
	 */
	struct Devious: public Smear::Device {
	
		Devious() {
			
			Param1 = TF1("f1","0.",1.e6);
			Param2 = TF2("f2","0.",1.e6);
			
			InputFunc = TF1("in","x",0.,1.e6);
			
			ParDim = 1;
			
			ParMin1 = 0.; ParMax1 = 1.e6;
			ParMin2 = 0.; ParMax2 = 1.e6;
			
			InFuncMin = 0.; InFuncMax = 1.e6;
			
			bInverse = true;
			
			InKin1 = kE;
			InKin2 = kTheta;
			OutKin = kE;
			
			name = "Arbitrary Input Device";
		}
		
		/**
		 Executes default constructor, then does SetInputFunction(f).
		 */
		Devious(TString f) {
			Devious();
			SetInputFunction(f);
		}
		
		/**
		 Executes default constructor, then sets input function to f, and parametrization to param.
		 */
		Devious(TString f, TString param) {
			Devious();
			SetInputFunction(f);
			SetParametrization(param);
		}
		
		TF1 InputFunc;
		
		double InFuncMin; double InFuncMax;
		bool bInverse;
		
		
		/**
		 Set the function f(X) of a kinematic variable X to be smeared.  X will still be stored in the output (not f(X)).
		 If your string contains "INVERSE" it will set f(X)=1/X and ensure fast smearing (i.e. it will not use root finding).
		 If you used the default constructor, this is already set to "INVERSE".
		 */
		void SetInputFunction(TString f) {
			KinType dumb;
			if (f.Contains("INVERSE") || f.Contains("inverse")) {
				bInverse = true;
				return;
			}
			bInverse = false;
			int d = ParseInputFunction(f,OutKin,dumb);
			if (d!=1) {
				std::cerr << "ERROR! Input function does not have dimension 1.  Function not applied.\n";
				return;
			}
			SetupRange(OutKin,InFuncMin,InFuncMax);
			InputFunc = TF1("in",f,InFuncMin,InFuncMax);
		}
		
		/**
		 Set the range of the input function.  If for some reason you must use this, you should call it after you call
		 SetInputFunction as a call to that automatically sets the range to appropriate values.
		 */
		void SetInputFunctionRange(double min, double max) {
			InFuncMin = min;  InFuncMax = max;
			InputFunc.SetRange(min,max);
		}
		
		/**
		 Set the kinematic variable of which a function will be smeared.  Alternatively, if you input a function string,
		 it can contain one of E,P,theat,phi,pT,pZ to automatically set this.
		 */
		void SetSmearedKinematics(KinType kin) {
			OutKin = kin;
			SetupRange(kin,InFuncMin,InFuncMax);
		}
		
		Devious *Clone() {
			Devious *dev = new Devious();
			*dev = *this;
			return dev;
		}
		
		/**
		 Device level smearing adapted for smearing functions of particle kinematics.  See documentation in Device.
		 */
		void DevSmear(Particle &prt, ParticleS &prtOut) {
			
			Double32_t x1; Double32_t x2;
			Double32_t y;
			
			if (Accept.Is(prt)) {
				
				x1 = SwitchKinGetFromParticle(prt,InKin1);
				x2 = SwitchKinGetFromParticle(prt,InKin2);
				
				if (bInverse) {
					y = 1./SwitchKinGetFromParticle(prt,OutKin);
				} else {
				  y = InputFunc.Eval(SwitchKinGetFromParticle(prt,OutKin));
				}
				
				if (ParDim==2) {
					y = Distribution.Generate(y,EvaluateRes(x1,x2));
				}
				else {
					y = Distribution.Generate(y,EvaluateRes(x1));
				}
				
				if (bInverse) {
				  SwitchKinStoreToParticle(prtOut,1./y,OutKin);
				} else {
					SwitchKinStoreToParticle(prtOut,InputFunc.GetX(y),OutKin);
				}
					
				//make sure angular coordinates live in S^2
				if (OutKin==kTheta) prtOut.theta = FixTopologyTheta(prtOut.theta);
				if (OutKin==kPhi) prtOut.phi = FixTopologyPhi(prtOut.phi);
				
				//make sure E, p are positive definite
				HandleBogusValues(prtOut,OutKin);
				
			} //if
		}
		
		ClassDef(Devious,1)
	};

	
}

// 2011/10/18: Changed SetParticle and DevSmear arguments to const Particle&
// to match Device.

#endif