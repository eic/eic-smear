//
// Smear.cxx
//
// Created by TB on 8/16/11.
// Copyright 2011 BNL. All rights reserved.
//

#include <cmath>

#include "Smear.h"

namespace Smear {
   
   Distributor::Distributor()
   : min(-1.e6)
   , max( 1.e6)
   , bias(0.)
   , plus(NAN)
   , minus(NAN)
   , bMoveable(false) {
/*      
      min = -1.e6;
      max = 1.e6;
      bias = 0.;
      
      plus = NAN;
      minus = NAN;
      
//      CustomDist = false;
      bMoveable = false;
      */
      Ran.SetSeed(0);
   }
   
/*   Distributor::~Distributor() {
      if(Distrib) {
         
      } // if
   }*/
   
   void Distributor::SetDistribution(const TString& s) {
      if (s.Contains("GAUSS")) {
//         CustomDist = false;
      } else {
//         if(Distrib) delete Distrib;
//         Distrib = new TF1("f",s,min,max);
         Distrib.reset(new TF1("f", s, min, max));
//         CustomDist = true;
      }
   }
   
   void Distributor::SetDistribution(Function f, int npars) {
//      if(Distrib) delete Distrib;
//      Distrib = new TF1("f",f,min,max,npars);
      Distrib.reset(new TF1("f", f, min, max, npars));
//      CustomDist = true;
   }
	
   void Distributor::SetBias(double b) {
      if(bias > 0.) {
         bias = b;
      } // if...
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
//      if (CustomDist && bMoveable) {
      if(Distrib.get()) {// and bMoveable) {
         Distrib->SetParameters(mean + bias, sigma);
         if(bMoveable) {
            y = Distrib->GetRandom(mean - minus, mean + plus);
         } // if...
//      }
//      else if (CustomDist) {
//      else if(Distrib) {
//         Distrib->SetParameters(mean+bias,sigma);
         else {
            y = Distrib->GetRandom();
         } // else
      }
      else {
         y = Ran.Gaus(mean + bias, sigma);
      }
      return y;
   }
   
} // namespace Smear
