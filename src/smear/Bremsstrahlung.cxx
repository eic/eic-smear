/**
 \file
 Implementation of class Smear::Bremsstrahlung.
 
 \author    Michael Savastio
 \date      2011-08-19
 \copyright 2011 Brookhaven National Lab
 */

#include "eicsmear/smear/Bremsstrahlung.h"

#include <cassert>
#include <cmath>

#include "eicsmear/erhic/VirtualParticle.h"

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
, mPdf(NULL) {
  Accept.AddParticle(11);
  Accept.AddParticle(-11);
}

Bremsstrahlung::Bremsstrahlung(const Bremsstrahlung& other)
    : Device(other), mParticle(NULL), mKMin(other.mKMin), mKMax(other.mKMax),
      mEpsilon(other.mEpsilon), mTraversed(other.mTraversed),
      mRadLength(other.mRadLength), mPdf(NULL) {
  // Duplicate cached particle if there is one
  // SetParticle() duplicates the particle and sets up the PDF
  //if (other.mParticle.get()) {
  //  SetParticle(*mParticle);
  //}  // if
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
  if (upper < lower || isnan(upper) || isnan(lower)) {
    return false;
  }  // if
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
  int n = static_cast<int>(ret);
  if (fabs(ret - n) < fabs(ret - n - 1)) {
    return n;
  } else {
    return n + 1;
  }  // if
}

void Bremsstrahlung::SetParticle(const erhic::VirtualParticle& prt) {
  mParticle.reset(static_cast<erhic::ParticleMC*>(prt.Clone()));
  if (!mPdf) {
    mPdf = new TF1("Smear_Bremsstrahlung_PDF",
                   this,
                   &Smear::Bremsstrahlung::dSigmadK,
                   mKMin, mKMax, 0);
  }  // if
  SetupPDF();
}

void Bremsstrahlung::FixParticleKinematics(ParticleMCS& prt) {
  prt.p = sqrt(prt.GetE() * prt.GetE() - prt.GetM() * prt.GetM());
  if (prt.p < 0. || isnan(prt.p)) prt.p = 0.;
  prt.pt = prt.p * sin(prt.theta);
  prt.pz = prt.p * cos(prt.theta);
}

Bremsstrahlung* Bremsstrahlung::Clone(Option_t* /* not used */) const {
  return new Bremsstrahlung(*this);
}

void Bremsstrahlung::Smear(const erhic::VirtualParticle& prt,
                           ParticleMCS& prtOut) {
  SetParticle(prt);
  const int nGamma = NGamma();
  for (int i = 0; i < nGamma; i++) {
    if (!SetupPDF()) break;
    double loss = mPdf->GetRandom();
    mParticle->E -= loss;
  }  // for
  prtOut.E = mParticle->GetE();
  FixParticleKinematics(prtOut);
  HandleBogusValues(prtOut);
  mParticle.reset(NULL);
}

}  // namespace Smear
