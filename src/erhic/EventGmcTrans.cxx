/**
 \file
 Implementation of class erhic::EventGmcTrans.
 
 \author    Thomas Burton
 \date      2012-09-20
 \copyright 2012 Brookhaven National Lab
 */

#include "eicsmear/erhic/EventGmcTrans.h"

#include <cmath>
#include <sstream>
#include <string>

namespace erhic {

EventGmcTrans::EventGmcTrans(const std::string& s)
// Initialise all floats to NAN for consistency with EventDis
: mStruckQuark(0)
, mQSquared(NAN)
, mBjorkenX(NAN)
, mInelasticity(NAN)
, mWSquared(NAN)
, mNu(NAN)
, mS(NAN)
, mZ(NAN)
, mHadronPt(NAN)
, mLeptonTheta(NAN)
, mLeptonPhi(NAN)
, mPhiSpin(NAN)
, mPhiHadron(NAN)
, mF1(NAN)
, mG1(NAN)
, mH1(NAN)
, mD1(NAN)
, mF1TPerp(NAN)
, mF1TPerp1(NAN)
, mF1TPerp12(NAN)
, mH1Perp(NAN)
, mH1Perp1(NAN)
, mH1Perp12(NAN)
, mAutSiv(NAN)
, mAutWtSiv(NAN)
, mAutSivAllQ(NAN)
, mAutWtSivAllQ(NAN)
, mAutSivPiDiff(NAN)
, mAutWtSivPiDiff(NAN)
, mAutCol(NAN)
, mAutWtCol(NAN)
, mAutTw3Col(NAN)
, mAutWtTw3Col(NAN)
, mAutColAllQ(NAN)
, mAutWtColAllQ(NAN)
, mXUnpolarised(NAN)
, mXSivers(NAN)
, mXCollins(NAN) {
  // Initialise from a string if provided.
  if (!s.empty()) {
    Parse(s);
  }  // if
}

bool EventGmcTrans::Parse(const std::string& line) {
  // Save ourselves the overhead of a new stringstream with each event read.
  static std::stringstream stream;
  // Clear the stream contents and flags from any previous use.
  stream.str("");
  stream.clear();
  // Read values from the input line.
  stream << line;
  stream
  // The first integer is always 0 (indicating start of event record)
  // so skip that and read the next int, which is the struck quark.
  >> mStruckQuark >> mStruckQuark
  >> mBjorkenX
  >> mQSquared
  >> mNu
  >> mInelasticity
  >> mWSquared
  >> mZ
  >> mHadronPt
  >> mLeptonTheta
  >> mLeptonPhi
  >> mPhiSpin
  >> mPhiHadron
  >> mF1 >> mG1 >> mH1 >> mD1
  >> mF1TPerp >> mF1TPerp1 >> mF1TPerp12
  >> mH1Perp >> mH1Perp1 >> mH1Perp12
  >> mAutSiv >> mAutWtSiv >> mAutSivAllQ >> mAutWtSivAllQ
  >> mAutSivPiDiff >> mAutWtSivPiDiff
  >> mAutCol >> mAutWtCol >> mAutTw3Col >> mAutWtTw3Col
  >> mAutColAllQ >> mAutWtColAllQ
  >> mXUnpolarised >> mXSivers >> mXCollins;
  // Inherits process from EventMC. Always set this equal to 99 for DIS
  process = 99;
  // The stream state should still be good if processing
  // the string went OK.
  return !stream.fail();
}

}  // namespace erhic
