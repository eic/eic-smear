#ifndef _ERHIC_BUILDTREE_SMEAR_
#define _ERHIC_BUILDTREE_SMEAR_

#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TRandom3.h>
#include <TF1.h>
#include <fstream>
#include <TF2.h>
#include <TLorentzVector.h>
#include <iostream>
#include <sstream>
#include <cmath>
#include <TMath.h>
#include <TVector2.h>
#include "VirtualParticle.h"
#include "Kinematics.h"
#include "ParticleIdentifier.h"

#include <memory>

#include "Particle.h"
#include "EventS.h"

namespace Smear {
	
   //	const double pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170680;
	const double pi = TMath::Pi();
   
	/**
	 Enumerator listing particle wise kinematic variables.  Naming is self explanitory.
	 These will be used when specifying arguments and outputs for device parametrizations.
	 */
	enum KinType { 
      kE, kP, kTheta, kPhi, kPz, kPt	
	};

   /** Classes of particles */
   enum EGenre {
      kAll = 0,
      kElectromagnetic = 1,
      kHadronic = 2
   };
   
	/**
	 This determines whether the particle in the argument is
	 
	 0: Not 1 or 2
	 1: Stable (i.e. in the Pythia sense) and a photon or lepton
	 2: Stable and a hadron.
	 */
	inline int PGenre(Particle prt) {
		int o;
      if (prt.GetStatus()==1 && abs(prt.Id())>10 && abs(prt.Id())<23 && abs(prt.Id())!=21) {
			o=1;
		}
		else if (prt.GetStatus() == 1 && abs(prt.Id())>110) {
			o=2;
		}
		else {
			o=0;
		}
		return o;
	}
   
	/**
	 Fix a polar angle so that it lies within [0,pi].
    TODO Nothing Smear-specific here - move to general functions file.
	 */
	inline double FixTopologyTheta(double theta) {
		while (theta < 0. || theta > pi) {
			if (theta<0.) {
				theta = -theta;
			}
			if (theta>pi) {
				theta = 2.*pi - theta;
			} 
		}
		return theta;
	}
	
	/**
	 Fix an azimuthal angle so that it lies within [0,2*pi).
    TODO Nothing Smear-specific here - move to general functions file.
	 */
	inline double FixTopologyPhi(double phi) {
      /*
       while (phi < 0. || phi >= 2.*pi) {
       if (phi < 0.) {
       phi = -phi;
       }
       if (phi >= 2.*pi) {
       phi = phi - 2.*pi;
       }
       }
       return phi;
       */
      return TVector2::Phi_0_2pi(phi);
	}
	
	/**
	 Returns the kinematic variable associated with kin from the input particle.
	 */
	inline double SwitchKinGetFromParticle(Particle prt, KinType kin) {
		double z;
		switch (kin) {
			case kE:
				z = prt.GetE(); break;
			case kP:
				z = prt.GetP(); break;
			case kTheta:
				z = prt.GetTheta(); break;
			case kPhi:
				z = prt.GetPhi(); break;
			case kPz:
				z = prt.GetPz(); break;
			case kPt:
				z = prt.GetPt(); break;
			default:
				z = prt.GetE(); break;
		}	
		return z;
	}
	
	/**
	 Stores z in the ParticleS.K where K is the kinematic variable associated with kin.
	 */
	inline void SwitchKinStoreToParticle(ParticleS &prt, double z, KinType kin) {
		switch (kin) {
			case kE:
				prt.SetE(z); break;
			case kP:
				prt.SetP(z); break;
			case kTheta:
				prt.SetTheta(z); break;
			case kPhi:
				prt.SetPhi(z); break;
			case kPz:
				prt.SetPz(z); break;
			case kPt:
				prt.SetPt(z); break;
			default:
            // TODO remove this default action!
				prt.SetE(z); break;
		}
	}
	
	/**
	 This dictates how the namespace deals with positive definite variables
    which have been smeared to negative values.
    Currently they are set to -999.
    TODO Provide static 'bogus value' for the user to specify.
	 */
	inline void HandleBogusValues(ParticleS &prt, KinType kin) {
		double fault = NAN;//-999.; 
      // TODO replace with switch statement
		if (prt.GetE() < 0. && kin==kE) {
			prt.SetE(fault);
		}
		if (prt.GetP() < 0. && kin==kP) {
			prt.SetP(fault);
		}
		if (prt.GetPt() < 0. && kin==kPt) {
			prt.SetPt(fault);
		}
		if (prt.GetPz() < 0. && kin==kPz) {
			prt.SetPz(fault);
		}	
	}
	
   /*
    TODO Make the "bogus value" static and settable
    */
	inline void HandleBogusValues(ParticleS& prt) {
		double fault = NAN;//-999.;
		if (prt.GetE() < 0.) {
			prt.SetE(fault);
		}
		if (prt.GetP() < 0.) {
			prt.SetP(fault);
		}
		if (prt.GetPt() < 0.) { 
			prt.SetPt(fault);
		}
		if (prt.GetPz() < 0.) {
			prt.SetPz(fault);
		}
	}
	
	//deteremines if a vector contains the obj object, returns position.  if not found, returns -1.
	template <class Tem> 
	int Vcontains(Tem obj,std::vector<Tem> v) {
		for (int i=0; i<v.size(); i++) {
			if (v.at(i)==obj) return i;
		}
		return -1;
	}
	
	inline bool IsCoreType(KinType kin) {
		if (kin==kE || kin==kP || kin==kTheta || kin==kPhi) return true;
		return false;
	}
	
	/**
	 Spawns a test Particle with some generic properties (electron, KS=1).  Useful for troubleshooting.
	 */
	inline Particle SpawnTestParticle() {
//		Particle p;
//		p.SetE(1.);
//		p.SetP(1.);
//		p.SetTheta(0.);
//		p.SetPhi(0.);
//		p.SetStatus(1);
//		p.SetId(11);
//		return p;
      // I KS id orig daughter ldaughter px py pz m E xv yv zv
      Particle p("0 0 0 0 0 0 0.0 0.0 1.0 0.0 1.0 0.0 0.0 0.0");
      p.ComputeDerivedQuantities();
      return p;
	}
	
	inline ParticleS SpawnTestParticleS() {
		ParticleS p;
		p.E = 1.;
		p.p = 1.;
		p.theta = 0.;
		p.phi = 0.;
		p.id = 11;
		return p;
	}
	
   /**
    
    */
	inline void SetupRange(KinType kin, double &min, double &max) {
		if (kin==kTheta) {
			min = 0.; max = pi;
		} else if (kin==kPhi) {
			min = 0.; max = 2.*pi;
		} else {
			min = 0.; max = 1.e6;
		}
	}
	
   int ParseInputFunction(TString &s, KinType &kin1, KinType &kin2);
   
	/**
	 Used by devices to generate smearing on some distribution.
    By default, smearing is generated on a Gaussian
	 of which the Monte Carlo value is the mean, and the standard deviation
    is given by the device parametrization.
	 It is intended that you access this object via the methods provided
    in the Device class.
    
    \todo This could do with some reworking, to make it compatible with
     an arbitrary function via a functor (which I think is the best way
     to go). Move to its own files.
	 */
	class Distributor {
      
   public:
      
		Distributor();
      
      virtual ~Distributor() { }
		
		double min;
		double max;
		double bias;
      
		double plus;
		double minus;
		
//		bool CustomDist;
		bool bMoveable;		
		
		TF1* Distrib;
		TRandom3 Ran;
		
		void SetDistribution(const TString& s);
		
      typedef double (*Function)(double*, double*);
//      typedef void* Function;
      void SetDistribution(Function, int npars);
         
		void SetBias(double b);
		
		void SetSeed(int n);
		
      /**
       Restricts range to [min, max].
       */
      void SetRange(double min, double max);
		
      /**
       Restricts range to [mean - minus, mean + plus].
       i.e. A variable range depending on the mean passed to Generate().
       plus and minus should both be greater than zero.
       Pass plus and/or minus <= 0 to revert to normal range.
       */
      void SetMoveableRange(double plus, double minus);
		
		void SetMoveable(bool b);
		
      double Generate(double mean, double sigma);
		
		ClassDef(Distributor, 1)
	};
   
}

#endif
