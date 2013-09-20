/**
 \file
 Declaration of class erhic::EventGmcTrans.
 
 \author    Thomas Burton
 \date      2012-09-20
 \copyright 2012 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_ERHIC_EVENTGMCTRANS_H_
#define INCLUDE_EICSMEAR_ERHIC_EVENTGMCTRANS_H_

#include <string>

#include <TLorentzVector.h>

#include "eicsmear/erhic/EventMC.h"

namespace erhic {
/**
 Describes an event from the gmc_trans Monte Carlo.
 */
class EventGmcTrans : public EventMC {
 public:
  /**
   Constructor.

   Optionally initialise from a string. See Parse() for the string format.
   */
  explicit EventGmcTrans(const std::string& = "");

  /**
   Destructor.
   */
  virtual ~EventGmcTrans() { }

  /**
   Set the event properties from a string.
   The string format should be as follows:
   \verbatim
   quark, x, Q2, nu, y, W2, z, hadronPt, leptonTheta, leptonPhi, phiSpin,
   phiHadron, f1, g1, h1, D1, f1Tperp, f1Tperp1, f1Tperp12 H1perp, H1perp1,
   H1perp12, AUT_Siv, AUT_WtSiv, AUT_Siv_allq, AUT_WtSiv_allq,
   AUT_Siv_piDiff, AUT_WtSiv_piDiff, AUT_Col, AUT_WtCol, AUT_Col_allq,
   AUT_WtCol_allq, AUT_Tw3_Col, AUT_WtTw3_Col, X_UU, X_SIV, X_COL
   \endverbatim
   Returns true if the string is processed successfully, false if not.
   */
  virtual bool Parse(const std::string&);

  /**
   Azimuthal angle of spin vector in gamma*-proton centre of mass frame.
   */
  Double_t GetPhiSpin() const { return mPhiSpin; }

  /**
   Azimuthal angle of hadron in gamma*-proton centre of mass frame.
   */
  Double_t GetPhiHadron() const { return mPhiHadron; }

  /**
   \remark Not called z() to avoid hiding Event::z().
   */
  Double_t GetHadronZ() const { return mZ; }

  Double_t GetHadronPt() const { return mHadronPt; }

  Double_t GetF1() const { return mF1; }

  Double_t GetG1() const { return mG1; }

  Double_t GetH1() const { return mH1; }

  Double_t GetD1() const { return mD1; }

  Double_t GetF1TPerp() const { return mF1TPerp; }

  Double_t GetF1TPerp1() const { return mF1TPerp1; }

  Double_t GetH1Perp() const { return mH1Perp; }

  Double_t GetH1Perp1() const { return mH1Perp1; }

  Double_t GetH1Perp12() const { return mH1Perp12; }

  Double_t GetSivers() const { return mAutSivAllQ; }

  Double_t GetSiversWeight() const { return mAutWtSivAllQ; }

  Double_t GetSiversStruckQuark() const { return mAutSiv; }

  Double_t GetSiversStruckQuarkWeight() const { return mAutWtSiv; }

  Double_t GetSiversPiDifference() const { return mAutSivPiDiff; }

  Double_t GetSiversPiDifferenceWeight() const { return mAutWtSivPiDiff; }

  Double_t GetCollins() const { return mAutColAllQ; }

  Double_t GetCollinsWeight() const { return mAutWtColAllQ; }

  Double_t GetCollinsStruckQuark() const { return mAutCol; }

  Double_t GetCollinsStructQuarkWeight() const { return mAutWtCol; }

  Double_t GetCollinsTwist3() const { return mAutTw3Col; }

  Double_t GetCollinsTwist3Weight() const { return mAutWtTw3Col; }

  Double_t GetCrossSectionUnpolarised() const { return mXUnpolarised; }

  Double_t GetCrossSectionSivers() const { return mXSivers; }

  Double_t GetCrossSectionCollins() const { return mXCollins; }

  /**
   Returns the hadron beam spin vector in the current frame.
   If the spin vector in the rest frame of the hadron (of mass m)
   is n and the 4-momentum of the hadron in the current frame
   is p = (E,P), then the hadron spin 4-vector in the current frame
   is (Bjorken & Drell, 1965):
   \verbatim
   s = (s0, S)
   \endverbatim
   where
   \verbatim
   s0 = n.P / m
   S  = n + (E/m - 1)(P.n) P / |P|^2
   \endverbatim
   */
  TLorentzVector GetHadronPolarisation() const { return TLorentzVector(); }

  /**
   Azimuthal angle of the produced hadron around the spin direction
   in the proton rest frame, as defined by the HERMES experiment.
   phi_s = [(qxk.S)/|qxk.S|]arccos[(qxk.qxS)/(|qxk||qxS|)],
   0 <= arccos <= pi.
   Gives result in range [0,2pi].
   */
  virtual double GetHermesPhiS() const { return 0.; }

 protected:
  Int_t mStruckQuark;  ///< Flavour of struck quark
  Double32_t mQSquared;  ///< Negative squared four-momentum of
                         ///< the exchanged boson
  Double32_t mBjorkenX;  ///< Longitudinal momentum fraction in the
                         ///< infinite momentum frame
  Double32_t mInelasticity;  // or mY?
  Double32_t mWSquared;  ///< Invariant mass of the hadronic final state
  Double32_t mNu;  ///< Energy of the exchanged boson in the hadron rest frame
  Double32_t mS;  ///< Square of the centre of mass energy
  Double32_t mZ;  ///< z of the produced hadron
  Double32_t mHadronPt;  ///< pT of the produced hadron
  Double32_t mLeptonTheta;  ///< Polar angle of the scattered lepton
  Double32_t mLeptonPhi;  ///< Azimuthal angle of the scattered lepton
  Double32_t mPhiSpin;  ///< Azimuthal angle of spin vector in gamma*-proton
                        ///< centre-of-mass frame.
  Double32_t mPhiHadron;
  Double32_t mF1;
  Double32_t mG1;
  Double32_t mH1;
  Double32_t mD1;
  Double32_t mF1TPerp;
  Double32_t mF1TPerp1;
  Double32_t mF1TPerp12;
  Double32_t mH1Perp;
  Double32_t mH1Perp1;
  Double32_t mH1Perp12;
  Double32_t mAutSiv;
  Double32_t mAutWtSiv;
  Double32_t mAutSivAllQ;
  Double32_t mAutWtSivAllQ;
  Double32_t mAutSivPiDiff;
  Double32_t mAutWtSivPiDiff;
  Double32_t mAutCol;
  Double32_t mAutWtCol;
  Double32_t mAutTw3Col;
  Double32_t mAutWtTw3Col;
  Double32_t mAutColAllQ;
  Double32_t mAutWtColAllQ;
  Double32_t mXUnpolarised;
  Double32_t mXSivers;
  Double32_t mXCollins;  ///< Cross sections

  ClassDef(erhic::EventGmcTrans, 1)
};

}  // namespace erhic

#endif  // INCLUDE_EICSMEAR_ERHIC_EVENTGMCTRANS_H_
