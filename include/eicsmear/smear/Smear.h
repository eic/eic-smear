#ifndef _EICSMEAR_SMEAR_SMEAR_H_
#define _EICSMEAR_SMEAR_SMEAR_H_

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
#include "eicsmear/erhic/VirtualParticle.h"
#include "eicsmear/erhic/Kinematics.h"
#include "eicsmear/erhic/ParticleIdentifier.h"

#include <memory>

#include "eicsmear/erhic/Particle.h"
#include "eicsmear/smear/ParticleMCS.h"

namespace Smear {

	/**
	 Enumerator listing particle wise kinematic variables.
    Naming is self explanitory.
	 These will be used when specifying arguments and outputs
    for device parametrizations.
	 */
	enum KinType { 
      kE, kP, kTheta, kPhi, kPz, kPt, kInvalidKinType
	};

   /** Classes of particles */
   enum EGenre {
      kAll = 0, kElectromagnetic = 1, kHadronic = 2
   };
   
   /** Particle charged **/
   enum ECharge {
      kNeutral, kCharged, kAllCharges
   };
   
	/**
    Determine particle "genre".

	 This refers to whether the particle in the argument is
	 "electromagnetic" or "hadronic" from the perspective of calorimetry.

    Returns one of the following (see enum EGenre):
	  kElectromagnetic: photon or electron/positron.
	  kHadronic:        stable hadron.
	  kAll:             neither of the above.
	 */
	inline int PGenre(const erhic::VirtualParticle& prt) {
		int genre(kAll);
      const int id = abs(prt.Id()); // Sign doesn't matter
      if(1 == prt.GetStatus()) { // Only check stable particles.
         if(id == 11 or id == 22) {
			   genre = kElectromagnetic;
		   } // if
		   else if(id >110) {
			   genre = kHadronic;
		   } // else if
      } // if
		return genre;
	}
   
	/**
	 Fix a polar angle so that it lies within [0,pi].
    TODO Nothing Smear-specific here - move to general functions file.
	 */
	inline double FixTheta(double theta) {
		while (theta < 0. || theta > TMath::Pi()) {
			if (theta<0.) {
				theta = -theta;
			}
			if (theta > TMath::Pi()) {
				theta = TMath::TwoPi() - theta;
			} 
		}
		return theta;
	}
	
	/**
	 Fix an azimuthal angle so that it lies within [0,2*pi).
    TODO Nothing Smear-specific here - move to general functions file.
	 */
	inline double FixPhi(double phi) {
      return TVector2::Phi_0_2pi(phi);
	}
	
	/**
	 Returns the kinematic variable associated with kin from the input particle.
	 */
	inline double GetVariable(const erhic::VirtualParticle& prt, KinType kin) {
		double z(0.);
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
            break;
		}	
		return z;
	}
	
	/**
	 Stores z in the ParticleS.K where K is the kinematic variable associated with kin.
	 */
	inline void SetVariable(ParticleMCS &prt, double z, KinType kin) {
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
            break;
		}
	}
	
	/**
	 This dictates how the namespace deals with positive definite variables
    which have been smeared to negative values.
	 */
	inline void HandleBogusValues(ParticleMCS &prt, KinType kin) {
		double fault(0.);
		if(kE == kin and prt.GetE() < 0.) {
			prt.SetE(fault);
		} // if
		else if(kP == kin and prt.GetP() < 0.) {
			prt.SetP(fault);
		} // else if
		else if(kPt == kin and prt.GetPt() < 0.) {
			prt.SetPt(fault);
		} // else if
		else if(kPz == kin and prt.GetPz() < 0.) {
			prt.SetPz(fault);
      } // else if
	}
	
   /*
    */
	inline void HandleBogusValues(ParticleMCS& prt) {
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

	inline bool IsCoreType(KinType kin) {
		if (kin==kE || kin==kP || kin==kTheta || kin==kPhi) return true;
		return false;
	}

   int ParseInputFunction(TString &s, KinType &kin1, KinType &kin2);
}

#endif // _EICSMEAR_SMEAR_SMEAR_H_
