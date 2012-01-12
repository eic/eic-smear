#ifndef _ERHIC_BUILDTREE_BREMSSTRAHLUNG_
#define _ERHIC_BUILDTREE_BREMSSTRAHLUNG_

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

#include "VirtualParticle.h"
#include "Kinematics.h"
#include "ParticleIdentifier.h"

#include "Smear.h"
#include "Acceptance.h"

#include "Device.h"

namespace Smear {

  struct Dubious: public Device {
		
		Dubious() {
			epsilon = 0.01;
			Traversed = 10.;
			RadLength = 47.1;
			
			Accept.AddParticle(11);
			Accept.AddParticle(-11);
		}
		
		Particle *P;
		
		TF1 PDF;
		
		double kMin;
		double kMax;
		double epsilon;
		double Traversed;
		double RadLength;
		
		double dSigmadK(double *x, double *par) {
			double k = x[0];
			par = NULL;
			double ret = 4./3.;
			ret += -4.*k/(3.*P->E);
			ret += pow(k/P->E,2);
			ret /= k;
			return ret;
		}
		
		void SetupPDF() {
			PDF = TF1("PDF",this,&Smear::Dubious::dSigmadK,kMin,kMax,0,"Smear::Dubious");
		}
		
		int NGamma() {
			double ret = 4.*log(kMax/kMin)/3.;
			ret += -4.*(kMax-kMin)/(3.*P->E);
			ret += 0.5*pow((kMax-kMin)/P->E,2);
			ret *= Traversed/RadLength;
			int n = (int) ret;
			if (fabs(ret-n)<fabs(ret-n-1)) {
				return n;
			} else {
				return n+1;
			}
		}
		
		void SetParticle(Particle &prt) {
			P = &prt;
			kMin = epsilon;
			kMax = P->E-epsilon;
			SetupPDF();
		}
		
		void FixParticleKinematics(ParticleS& prt) {
			prt.p = sqrt(prt.GetE()*prt.GetE() - prt.M()*prt.M());
			prt.pt = prt.p*sin(prt.theta);
			prt.pz = prt.p*cos(prt.theta);
		}
		
		Dubious *Clone() {
			Dubious *dev = new Dubious();
			*dev = *this;
			return dev;
		}
		
     void DevSmear(Particle &prt, ParticleS& prtOut) {
			
			//double before;  double after;
			
			//before = prt.E;
			//std::cout << "before " << before << std::endl;
			
			SetParticle(prt);
			
			for (int i=0; i<NGamma(); i++) {
				prt.E = prt.E - PDF.GetRandom();
			}

			// Is this right?
         // Should it not be prtOut that is modified?
//			FixParticleKinematics(prt);
        FixParticleKinematics(prtOut);
//			HandleBogusValues(prt);
        HandleBogusValues(prtOut);
			
			//after = prt.E;
			//std::cout << "after " << after << std::endl;
			
			//if (after>before) std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!WTF!!!!!!!!!!!!!!!!!!!!!!!!\n";
		}
		
		ClassDef(Dubious,1)
	};

}

#endif