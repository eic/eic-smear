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

#include "eicsmear/erhic/VirtualEvent.h"
#include "eicsmear/erhic/VirtualParticle.h"
#include "eicsmear/erhic/VirtualEvent.h"
#include "eicsmear/erhic/Kinematics.h"
#include "eicsmear/erhic/ParticleIdentifier.h"
#include "eicsmear/smear/Acceptance.h"
#include "eicsmear/smear/Device.h"
#include "eicsmear/smear/ParticleMCS.h"
#include "eicsmear/smear/Smear.h"
#include "eicsmear/smear/Smearer.h"

namespace erhic {
   class VirtualParticle;
} // namespace erhic

namespace Smear {
	
   /**
	 Specialized Device class for simulating tracking detectors.
    Contains a special parametrization suited for tracking
	 detectors and methods for setting the various input parameters
    (such as magnetic field and track length).
	 If you want to simulate a tracker by entering a parametrization direclty,
    you will have to use the genereic Device class.
    \todo Remove SetParticle(), P and bMoreCrit - can pass Particle as function
          argument when required, as have a ComputePathLength(Particle) method
          serving bMoreCrit's function.
	 */
   struct Tracker: public Smearer {

      /** Default constructor */
      Tracker();

      /**
       Constructor.
       Initialise with all detector characteristics.
       */
		Tracker(double inner, double outer, double L, double Bb,
              double nrl, double sigmaRPhi, double N);		

      /**
       Constructor.
       Initialise with the detector dimensions and default values for
       other characteristics.
       */
		Tracker(double inner, double outer, double L);

      /** Destructor */
      virtual ~Tracker();

      /**
       Returns a new copy of this Tracker.
       The argument is not used.
      */
		virtual Tracker* Clone(const char*) const;

      /**
       Smear the kinematics of prt and store the result in prtOut.
      */
		void DevSmear(const erhic::VirtualParticle& prt, ParticleMCS& prtOut);

		/**
		 Set the magnetic field in the tracker in teslas.
		 */
		void SetMagneticField(double Bb);

		/**
		 Set the inner and outer radius of your cylindrical tracker in meters.
       Automatically calls FixAcceptance, so if you want a custom acceptance
       you should call it AFTER calling this.
		 */
//		void SetRadii(double inner, double outer);
		
		/**
		 Set the total length of your tracker (centered at IP) in meters.
       Automatically calls FixAcceptance, so if you want a custom acceptance
       you should call it AFTER calling this.
		 */
		void SetLength(double L);

		/**
		 Set the radii and length in meters.
       Automatically calls FixAcceptance, so if you want a custom acceptance,
       you should set it AFTER calling this.
		 */
		void SetDimensions(double inner, double outer, double L);
      
		/**
		 Set the number of radiation lengths in the tracker.
		 */
		void SetNumberOfRadiationLengths(double);

		/**
		 Set the position resolution in meters.
		 */
		void SetPositionResolution(double);

		/**
		 Set the number of measurements.
		 */
		void SetNumberOfMeasurements(double);

		double beta();

		double cosSquaredgamma();

		double L();

      double LPrime();

		double tanTheta();

		double ThetaCrit();

		/**
		 Returns the minimum theta of particle accepted by a tracker with
       the length and radii you specified.
		 If you must set a different acceptance in theta, it is recomended
       that you allow it no closer to the beam line than this.
		 */
		double GetThetaMin();

		void SetParticle(const erhic::VirtualParticle& prt);

		void SetVals(double inner, double outer, double L, double Bb,
                   double nrl, double sigmaRPhi, double N);

   protected:

		double MultipleScatteringContribution();
      
		double IntrinsicContribution();
      
		double EvaluateRes();

		/**
		 Automatically set acceptance to that of a cylinder with the dimensions
       you set.
		 */
		void FixAcceptance();

		double B; ///< Magnetic field strength in Tesla
		double Nrl; ///< Number of radiation lengths (dimensionless)
		double SigmaRPhi; ///< Point resolution (
		double n; ///< Number of fit points
		double R1; ///< Inner radius (cm)
		double R2; ///< Outer radius (cm)
		double l; ///< Total tracker length (cm)
		double fault;
		bool bMoreCrit;
      Distributor Distribution; ///< Random distribution
		const erhic::VirtualParticle *P;
      KinType OutKin;

		ClassDef(Tracker,1)
	};

   inline void Tracker::SetMagneticField(double Bb) {
      B = Bb;
   }

   inline void Tracker::SetLength(double L) {
      l = L;
      FixAcceptance();
   }

//   inline void Tracker::SetDimensions(double inner, double outer, double L) {
//      SetRadii(inner, outer); // Calls FixAcceptance
//      SetLength(L); // Calls FixAcceptance again
//   }

   inline double Tracker::cosSquaredgamma() {
      if(not P) {
         return fault;
      } // if
      return cos(pi / 2. - P->GetTheta()) * cos(pi / 2. - P->GetTheta());
   }

   inline double Tracker::tanTheta() {
      if(not P) {
         return fault;
      } // if
      return fabs(tan(P->GetTheta()));
   }
   
   inline double Tracker::ThetaCrit() {
      return atan(2. * R2 / l);
   }

   inline void Tracker::SetNumberOfRadiationLengths(double nrl) {
      Nrl = nrl;
   }

   inline double Tracker::beta() {
      if(not P) {
         return fault;
      } // if
      return P->GetP() / P->GetE();
   }

   inline void Tracker::SetPositionResolution(double sigmaRPhi) {
      SigmaRPhi = sigmaRPhi;
   }

   inline void Tracker::SetNumberOfMeasurements(double N) {
      n = N;
   }

   inline double Tracker::GetThetaMin() {
      return atan(2. * R1 / l);
   }

//   inline void Tracker::FixAcceptance() {
//      Accept.SetTheta(GetThetaMin(), pi - GetThetaMin());
//   }

   inline Tracker* Tracker::Clone(const char*) const {
      return new Tracker(*this);
   }

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
		void DevSmear(Particle &prt, ParticleMCS &prtOut) {
			
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