/**
 \file
 Implementation of class Bremsstrahlung.
 
 \author Michael Savastio
 */

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
   , mPdf("Smear_Bremsstrahlung_PDF",
         this,
         &Smear::Bremsstrahlung::dSigmadK,
         mKMin, mKMax, 0)
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
   
   
   void Bremsstrahlung::SetupPDF() {
      mPdf.SetRange(mKMin, mKMax);
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
   
   
   void Bremsstrahlung::SetParticle(Particle &prt) {
      mParticle = &prt;
      mKMin = mEpsilon;
      mKMax = mParticle->E - mEpsilon;
      SetupPDF();
   }
   
   
   void Bremsstrahlung::FixParticleKinematics(ParticleS& prt) {
      prt.p = sqrt(prt.GetE() * prt.GetE() - prt.M() * prt.M());
      prt.pt = prt.p * sin(prt.theta);
      prt.pz = prt.p * cos(prt.theta);
   }
   
   
   Bremsstrahlung* Bremsstrahlung::Clone() {
      return new Bremsstrahlung(*this);
   }
   
   
   void Bremsstrahlung::DevSmear(Particle &prt, ParticleS& prtOut) {
      
      SetParticle(prt);
      
      for(int i = 0; i < NGamma(); i++) {
         prt.E = prt.E - mPdf.GetRandom();
      } // for
      
      FixParticleKinematics(prtOut);
      HandleBogusValues(prtOut);
      
      mParticle = NULL;
   }

} // namespace Smear
