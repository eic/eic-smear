//
// Smear.cxx
//
// Created by TB on 8/16/11.
// Copyright 2011 BNL. All rights reserved.
//

#include "Smear.h"

#include <cmath>
#include <algorithm>
#include <functional>
#include <list>

namespace Smear {

   Distributor::Distributor()
   : min(-1.e6)
   , max( 1.e6)
   , bias(0.)
   , plus(NAN)
   , minus(NAN)
   , bMoveable(false)
   , Distrib(NULL) {
      Ran.SetSeed(0);
   }

   void Distributor::SetDistribution(const TString& s) {
      if(not s.Contains("GAUSS")) {
         Distrib = new TF1("f", s, min, max);
      } // if
   }

   void Distributor::SetDistribution(Function f, int npars) {
      Distrib = new TF1("f",f,min,max,npars);
   }

   void Distributor::SetBias(double b) {
      if(b > 0.) {
         bias = b;
      } // if
      else {
         bias = 0.;
      } // else
   }

   void Distributor::SetSeed(int n) {
      Ran.SetSeed(n);
   }

   void Distributor::SetRange(double Min, double Max) {
      min = Min;
      max = Max;
      Distrib->SetRange(min,max);
   }

   void Distributor::SetMoveableRange(double aplus, double aminus) {
      if(aplus > 0. and aminus > 0.) {
         plus = aplus;
         minus = aminus;
         bMoveable = true;
      } // if
      else {
         plus = 0.;
         minus = 0.;
         bMoveable = false;
      } // else
   }

   void Distributor::SetMoveable(bool b) {
      bMoveable = b;
   }

   double Distributor::Generate(double mean, double sigma) {
      double y;
      if(Distrib) {
         Distrib->SetParameters(mean + bias, sigma);
         if(bMoveable) {
            y = Distrib->GetRandom(mean - minus, mean + plus);
         } // if
         else {
            y = Distrib->GetRandom();
         } // else
      } // if
      else {
         y = Ran.Gaus(mean + bias, sigma);
      } // else
      return y;
   }

	int ParseInputFunction(TString &s, KinType &kin1, KinType &kin2) {
		int d=0;
		if (s.Contains("E")) {
			s.ReplaceAll("E","x");
			d++;
			kin1 = kE;
		}
		if (s.Contains("P")) {
			d++;
			if (d==2) {
				s.ReplaceAll("P","y");
				kin2 = kP;
			} else {
				s.ReplaceAll("P","x");
				kin1 = kP;
			}
		}
		if (s.Contains("theta") || s.Contains("Theta")) {
			d++;
			if (d==2) {
				s.ReplaceAll("theta","y"); s.ReplaceAll("Theta","y");
				kin2 = kTheta;
			} else {
				s.ReplaceAll("theta","x"); s.ReplaceAll("Theta","x");
				kin1 = kTheta;
			}
		}
		if (s.Contains("phi")) {
			d++;
			if (d==2) {
				s.ReplaceAll("phi","y"); 
				kin2 = kPhi;
			} else {
				s.ReplaceAll("phi","x"); 
				kin1 = kPhi;
			}
		}
		if (s.Contains("pT")) {
			d++;
			if (d==2) {
				s.ReplaceAll("pT","y"); 
				kin2 = kPt;
			} else {
				s.ReplaceAll("pT","x"); 
				kin1 = kPt;
			}
		}
		if (s.Contains("pZ")) {
			d++;
			if (d==2) {
				s.ReplaceAll("pZ","y"); 
				kin2 = kPz;
			} else {
				s.ReplaceAll("pZ","x"); 
				kin1 = kPz;
			}
		}
		if (d>2 || d<0) {
			std::cout << "!WARNING! Not enough, or too many parameters detected in parametrization (d=";
			std::cout << d << ").\n";
		}
		if (d==0) d++;
		return d;
	}
} // namespace Smear
