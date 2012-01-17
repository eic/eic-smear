/**
 \file
 Implementation of class Bremsstrahlung.
 
 \author Michael Savastio
 */

#include <cassert>
#include <cmath>

#include "Bremsstrahlung.h"

namespace Smear {
   
   
   Bremsstrahlung::Bremsstrahlung(double epsilon,
                                  double traversed,
                                  double radLength)
   : mParticle(NULL)
   , mKMin(0.)
   , mKMax(0.)
   , mEpsilon(epsilon)
   , mTraversed(traversed)
   , mRadLength(radLength)
   , mPdf(NULL)
   {
      Accept.AddParticle(11);
      Accept.AddParticle(-11);
   }
   
   
   double Bremsstrahlung::dSigmadK(double *x, double*) {
      double k = x[0];
      double ret = 4. / 3.;
      ret += -4. * k / (3. * mParticle->E);
      ret += pow(k / mParticle->E, 2.);
      ret /= k;
      return ret;
   }
   
   
   bool Bremsstrahlung::SetupPDF() {
      
      double lower = mEpsilon;
      double upper = mParticle->E - mEpsilon;
      
      if(upper < lower or isnan(upper) or isnan(lower)) return false;
      
      mKMin = lower;
      mKMax = upper;
      mPdf->SetRange(mKMin, mKMax);
      
      return true;
   }
   
   
   int Bremsstrahlung::NGamma() {
      double ret = 4. * log(mKMax / mKMin) / 3.;
      ret += -4. * (mKMax - mKMin) / (3. * mParticle->E);
      ret += 0.5* pow((mKMax - mKMin) / mParticle->E, 2.);
      ret *= mTraversed / mRadLength;
      int n = (int)ret;
      if(fabs(ret - n) < fabs(ret - n - 1)) {
         return n;
      } // if
      else {
         return n + 1;
      } // else
   }
   
   
   void Bremsstrahlung::SetParticle(const Particle& prt) {
      mParticle.reset((Particle*)prt.Clone());
      if(not mPdf) {
         mPdf = new TF1("Smear_Bremsstrahlung_PDF",
                        this,
                        &Smear::Bremsstrahlung::dSigmadK,
                        mKMin, mKMax, 0);
      } // if
      SetupPDF();
   }
   
   
   void Bremsstrahlung::FixParticleKinematics(ParticleS& prt) {
      prt.p = sqrt(prt.GetE() * prt.GetE() - prt.GetM() * prt.GetM());
      if(prt.p < 0. or isnan(prt.p)) prt.p = 0.;
      prt.pt = prt.p * sin(prt.theta);
      prt.pz = prt.p * cos(prt.theta);
   }
   
   
   Bremsstrahlung* Bremsstrahlung::Clone() {
      return new Bremsstrahlung(*this);
   }
   
   
   void Bremsstrahlung::DevSmear(const Particle &prt, ParticleS& prtOut) {
      
      SetParticle(prt);
      
      const int nGamma = NGamma();
      
      for(int i = 0; i < nGamma; i++) {
         if(not SetupPDF()) break;
         double loss = mPdf->GetRandom();
         mParticle->E -= loss;
      } // for
      
      prtOut.E = mParticle->GetE();
      
      FixParticleKinematics(prtOut);
      HandleBogusValues(prtOut);
      
      mParticle.reset(NULL);
   }

} // namespace Smear
