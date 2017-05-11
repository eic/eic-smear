/**
 \file
 Declaration of class erhic::EventRapgap.
 
 \author    Thomas Burton
 \date      2011-08-31
 \copyright 2011 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_ERHIC_EVENTPYTHIA_H_
#define INCLUDE_EICSMEAR_ERHIC_EVENTPYTHIA_H_

#include <string>

#include <Rtypes.h>

#include "eicsmear/erhic/EventMC.h"

namespace erhic {

/**
 \brief Pythia 6 DIS event
 
 Describes an event from the EIC PYTHIA6 implementation.
 
 \todo Change sHat/t_hat/u_hat naming to be consistent
 */
class EventPythia : public EventMC {
 public:
  /**
   Constructor.
   @param [in] str A text string setting event-wise quantities. The string
   format should be (no newlines):
   "I ievent genevent subprocess nucleon targetparton xtargparton
   beamparton xbeamparton thetabeamprtn truey trueQ2 truex trueW2 trueNu
   leptonphi s_hat t_hat u_hat pt2_hat Q2_hat F2 F1 R sigma_rad SigRadCor
   EBrems photonflux nrTracks"
   */
  explicit EventPythia(const std::string& str = "");

  /**
   Destructor
   */
  virtual ~EventPythia();

  /**
   Parses the event information from a text string.
   See the constructor for the string format.
   Returns true in the event of a successful read operation,
   false in case of an error.
   */
  virtual bool Parse(const std::string&);

  /**
   Sets the nucleon species.
   @param n [in] PDG code of the hadron beam, see MSTI(12)
   */
  virtual void SetNucleon(int n);

  /**
   Sets the target parton species.
   @param n [in] PDG code of the struck parton in the hadron beam,
   see MSTI(16)
   */
  virtual void SetTargetParton(int n);

  /**
   Sets the beam parton species.
   @param n [in] PDG code of the parton interacting with the hadron beam in
   the case of resolved photon processes and soft VMD, see MSTI(15)
   */
  virtual void SetBeamParton(int n);

  /**
   Sets the number of trials required to generate this event
   @param n [in] The number of trials
   */
  virtual void SetGenEvent(int n);

  /**
   Sets the x of the target parton.
   */
  virtual void SetTargetPartonX(double xB);

  /**
   Sets the x of the beam parton.
   */
  virtual void SetBeamPartonX(double xB);

  /**
   Sets the polar angle of the beam parton in the cm frame
   \todo enforce range [0, pi] when setting
   */
  virtual void SetBeamPartonTheta(double radians);

  /**
   Azimuthal angle of the scattered lepton in the cm frame
   \todo enforce range [0, 2pi) when setting
   */
  virtual void SetLeptonPhi(double radians);

  /**
   Used for radiative corrections.
   */
  virtual void SetF1(double f1);

  /**
   Used for radiative corrections.
   */
  virtual void SetF2(double f2);

  /**
   Used for radiative corrections.
   */
  virtual void SetSigmaRad(double sr);

  /**
   Sets the Mandelstamm s of the hard interaction.
   */
  virtual void SetHardS(double s);

  /**
   Sets the Mandelstamm t of the hard interaction.
   */
  virtual void SetHardT(double t);

  /**
   Sets the Mandelstamm u of the hard interaction.
   */
  virtual void SetHardU(double u);

  /**
   Sets the Q<sup>2</sup> of the hard interaction.
   */
  virtual void SetHardQ2(double Q2);

  /**
   Sets the squared p<sub>T</sub> of the hard interaction.
   */
  virtual void SetHardPt2(double pt2);

  /**
   Used for radiative corrections.
   */
  virtual void SetSigRadCor(double s);

  /**
   Sets the energy per radiative photon in the nuclear rest frame.
   */
  virtual void SetEBrems(double e);

  /**
   Flux factor, see VINT(319) in PYTHIA 6.
   */
  virtual void SetPhotonFlux(double f);

  /**
   Sets the true (not reconstructed) value for inelasticity.
   */
  virtual void SetTrueY(double inelasticity);

  /**
   Sets the true (not reconstructed) value for Q<sup>2</sup>.
   */
  virtual void SetTrueQ2(double Q2);

  /**
   Sets the true (not reconstructed) value for x.
   */
  virtual void SetTrueX(double x);

  /**
   Sets the true (not reconstructed) value for W<sup>2</sup>.
   */
  virtual void SetTrueW2(double W2);

  /**
   Sets the true (not reconstructed) value for nu.
   */
  virtual void SetTrueNu(double Nu);

  /**
   Used for radiative corrections.
   */
  virtual void SetR(double r);

  /**
   Returns the number of trials required to generate this event.
   */
  virtual int GetGenEvent() const;

  /**
   Returns the x of the target parton.
   */
  virtual double GetTargetPartonX() const;

  /**
   Returns the x of the beam parton.
   */
  virtual double GetBeamPartonX() const;

  /**
   Returns the polar angle of the beam parton in the cm frame,
   in radians in the range [0, pi]
   */
  virtual double GetBeamPartonTheta() const;

  /**
   Returns the azimuthal angle of the scattered lepton.
   
   Angle is given in the centre-of-mass frame, in radians in the range [0, 2pi).
   */
  virtual double GetLeptonPhi() const;

  /**
   Used for radiative corrections.
   */
  virtual double GetF1() const;

  /**
   Used for radiative corrections.
   */
  virtual double GetF2() const;

  /**
   Used for radiative corrections.
   */
  virtual double GetSigmaRad() const;

  /**
   Returns the Mandelstamm s of the hard interaction.
   */
  virtual double GetHardS() const;

  /**
   Returns the Mandelstamm t of the hard interaction.
   */
  virtual double GetHardT() const;

  /**
   Returns the Mandelstamm u of the hard interaction.
   */
  virtual double GetHardU() const;

  /**
   Returns the Q<sup>2</sup> of the hard interaction.
   */
  virtual double GetHardQ2() const;

  /**
   Returns the squared p<sub>T</sub> of the hard interaction.
   */
  virtual double GetHardPt2() const;

  /**
   Used for radiative corrections.
   */
  virtual double GetSigRadCor() const;

  /**
   Returnss the energy per radiative photon in the nuclear rest frame.
   */
  virtual double GetEBrems() const;

  /**
   Returns the flux factor, see VINT(319) in PYTHIA 6.
   */
  virtual double GetPhotonFlux() const;

  /**
   Sets the true (not reconstructed) value for inelasticity.
   */
  virtual double GetTrueY() const;

  /**
   Sets the true (not reconstructed) value for Q<sup>2</sup>.
   */
  virtual double GetTrueQ2() const;

  /**
   Sets the true (not reconstructed) value for x.
   */
  virtual double GetTrueX() const;

  /**
   Sets the true (not reconstructed) value for W<sup>2</sup>.
   */
  virtual double GetTrueW2() const;

  /**
   Sets the true (not reconstructed) value for &nu;.
   */
  virtual double GetTrueNu() const;

  /**
   Used for radiative corrections.
   */
  virtual double GetR() const;

  /**
   Returns the scattered lepton.
   This is the first final state particle with the same species
   as the beam lepton and parent index equal to three
   (counting index from 1).
   */
  virtual const ParticleMC* ScatteredLepton() const;

  // Let them all be public; this access method dances for POD does not make sense;
  //protected:
  // Inline comments after field names will appear in ROOT
  // when EventPythia::Dump() is called.
  Int_t       nucleon;          ///< PDG code of the hadron beam,
                                ///< see MSTI(12)
  Int_t       tgtparton;        ///< PDG code of the struck parton
                                ///< in the hadron beam, see MSTI(16)
  Int_t       beamparton;       ///< Parton interacting with the hadron
                                ///< beam in the case of resolved photon
                                ///< processes and soft VMD, see MSTI(15)
  Int_t       genevent;         ///< Trials required for this event
  Double32_t  xtgtparton;       ///< Momentum fraction taken by the
                                ///< target parton, see PARI(34)
  Double32_t  xbeamparton;      ///< Momentum fraction taken by the
                                ///< beam parton, see PARI(33)
  Double32_t  thetabeamparton;  ///< Polar angle of the beam parton in
                                ///< the cm frame, between 0 and pi
                                ///< radians, see PARI(53)
  Double32_t  leptonphi;        ///< Azimuthal angle of the scattered
                                ///< lepton in the cm frame
  Double32_t  F1;               ///< Value used for radiative corrections
  Double32_t  sigma_rad;        ///< Value used for radiative corrections
  Double32_t  t_hat;            ///< Mandelstam t of the hard subprocess,
                                ///< see PARI(15)
  Double32_t  u_hat;            ///< Mandelstam u of the hard subprocess,
                                ///< see PARI(16)
  Double32_t  Q2_hat;           ///< Q<sup>2</sup> of the hard subprocess,
                                ///< see PARI(22)
  Double32_t  SigRadCor;        ///< Value used for radiative corrections
  Double32_t  EBrems;           ///< Energy per radiative photon in the
                                ///< nuclear rest frame.
  Double32_t  photonflux;       ///< Flux factor, see VINT(319)
  Double32_t  trueY;            ///< Generated y of the event,
                                ///< see VINT(309)
  Double32_t  trueQ2;           ///< Generated Q<sup>2</sup> of the event,
                                ///< see VINT(307)
  Double32_t  trueX;            ///< Generated x of the event
  Double32_t  trueW2;           ///< Generated W<sup>2</sup> of the event
  Double32_t  trueNu;           ///< Generated nu of the event
  Double32_t  F2;               ///< Value used for radiative corrections
  Double32_t  R;                ///< Value used for radiative corrections
  Double32_t  pt2_hat;          ///< Squared p<sub>T</sub> of the hard
                                ///< subprocess, see PARI(18)
  Double32_t  sHat;             ///< Mandelstam s of the hard subprocess,
                                ///< see PARI(14)
  ClassDef(erhic::EventPythia, 2)
};

inline void EventPythia::SetNucleon(int n) {
  nucleon = n;
}

inline void EventPythia::SetTargetParton(int n) {
  tgtparton = n;
}

inline void EventPythia::SetBeamParton(int n) {
  beamparton = n;
}

inline void EventPythia::SetGenEvent(int n) {
  genevent = n;
}

inline void EventPythia::SetTargetPartonX(double xB) {
  xtgtparton = xB;
}

inline void EventPythia::SetBeamPartonX(double xB) {
  xbeamparton = xB;
}

inline void EventPythia::SetBeamPartonTheta(double radians) {
  thetabeamparton = radians;
}

inline void EventPythia::SetLeptonPhi(double radians) {
  leptonphi = radians;
}

inline void EventPythia::SetF1(double f1) {
  F1 = f1;
}

inline void EventPythia::SetF2(double f2) {
  F2 = f2;
}

inline void EventPythia::SetSigmaRad(double sr) {
  sigma_rad = sr;
}

inline void EventPythia::SetHardS(double s) {
  sHat = s;
}

inline void EventPythia::SetHardT(double t) {
  t_hat = t;
}

inline void EventPythia::SetHardU(double u) {
  u_hat = u;
}

inline void EventPythia::SetHardQ2(double Q2) {
  Q2_hat = Q2;
}

inline void EventPythia::SetHardPt2(double pt2) {
  pt2_hat = pt2;
}

inline void EventPythia::SetSigRadCor(double s) {
  SigRadCor = s;
}

inline void EventPythia::SetEBrems(double e) {
  EBrems = e;
}

inline void EventPythia::SetPhotonFlux(double f) {
  photonflux = f;
}

inline void EventPythia::SetTrueY(double inelasticity) {
  trueY = inelasticity;
}

inline void EventPythia::SetTrueQ2(double Q2) {
  trueQ2 = Q2;
}

inline void EventPythia::SetTrueX(double xB) {
  trueX = xB;
}

inline void EventPythia::SetTrueW2(double W2) {
  trueW2 = W2;
}

inline void EventPythia::SetTrueNu(double Nu) {
  trueNu = Nu;
}

inline void EventPythia::SetR(double r) {
  R = r;
}

inline int EventPythia::GetGenEvent() const {
  return genevent;
}

inline double EventPythia::GetTargetPartonX() const {
  return xtgtparton;
}

inline double EventPythia::GetBeamPartonX() const {
  return xbeamparton;
}

inline double EventPythia::GetBeamPartonTheta() const {
  return thetabeamparton;
}

inline double EventPythia::GetLeptonPhi() const {
  return leptonphi;
}

inline double EventPythia::GetF1() const {
  return F1;
}

inline double EventPythia::GetF2() const {
  return F2;
}

inline double EventPythia::GetSigmaRad() const {
  return sigma_rad;
}

inline double EventPythia::GetHardS() const {
  return sHat;
}

inline double EventPythia::GetHardT() const {
  return t_hat;
}

inline double EventPythia::GetHardU() const {
  return u_hat;
}

inline double EventPythia::GetHardQ2() const {
  return Q2_hat;
}

inline double EventPythia::GetHardPt2() const {
  return pt2_hat;
}

inline double EventPythia::GetSigRadCor() const {
  return SigRadCor;
}

inline double EventPythia::GetEBrems() const {
  return EBrems;
}

inline double EventPythia::GetPhotonFlux() const {
  return photonflux;
}

inline double EventPythia::GetTrueY() const {
  return trueY;
}

inline double EventPythia::GetTrueQ2() const {
  return trueQ2;
}

inline double EventPythia::GetTrueX() const {
  return trueX;
}

inline double EventPythia::GetTrueW2() const {
  return trueW2;
}

inline double EventPythia::GetTrueNu() const {
  return trueNu;
}

inline double EventPythia::GetR() const {
  return R;
}


class EventBeagle : public EventPythia {
 public:
  explicit EventBeagle(const std::string& str = "");

  ~EventBeagle();

  bool RequiresEaParticleFields() { return true; };

  bool Parse(const std::string&);

  ///=================additional variables for BeAGLE==================
	//put in the public region to be easily accessed
  Int_t lepton;  ///< lepton beam ID
  Int_t Atarg;	 ///< mass number of target beam
  Int_t Ztarg;	 ///< charge number of target beam
  Double32_t pzlep;		///< lepton beam momentum
  Double32_t pztarg;		///< target beam momentum
  Double32_t pznucl;		///< target nucleon momentum
	Double32_t crang;			///< crossing angle (mr). crang=1000*atan(px/pz), one beam px=py=0, the other py=0
	Double32_t crori;			///< crossing angle orientation, +-1 lepton beam along +-z, +-2 hadron beam along +-z, 0 meas no crossing angle
  Double32_t b;		///< impact parameter value
  Double32_t Phib;		///< phi of impact parameter vector
	Double32_t Thickness; ///< T(b) in nucleons/fm^2
	Double32_t ThickScl; ///< T(b)/rho0 in fm
	Int_t Ncollt; ///< Number of collisions in target
	Int_t Ncolli; ///< Number of inelastic collisions in target
	Int_t Nwound; ///< Number of wounded nucleon including those in INC
	Int_t Nwdch; ///< Number of wounded proton including those in INC
	Int_t Nnevap; ///< Number of neutrons from evaporation
	Int_t Npevap; ///< Number of protons from evaporation
	Int_t Aremn; ///< A of the nuclear remnant after evaporation and breakup
  Int_t NINC; ///< Number of stable hadrons from intranuclear cascade
	Int_t NINCch; ///< Number of charged stable hadrons from intranuclear cascade
	Double32_t d1st; ///< density-weighted distance from first collision to the edge of the nucleus (amount of material traversed / rho0)
	Double32_t davg; ///< Average density-weighted distance from all inelastic collisions to the edge of the nucleus
	Double32_t pxf; ///< Sum fermi momentum of all inelastic participant in target rest frame z along gamma*
	Double32_t pyf; ///< Sum fermi momentum of all inelastic participant in target rest frame z along gamma*
	Double32_t pzf; ///< Sum fermi momentum of all inelastic participant in target rest frame z along gamma*
	Double32_t Eexc; ///< Excitation energy in the nuclear remnant before evaporation and breakup 
	Double32_t RAevt; ///< Nuclear PDF ratio for the up sea for the given event kinematics 
	Double32_t User1; ///< User variables to prevent/delay future format changes 
	Double32_t User2; ///< User variables to prevent/delay future format changes 
	Double32_t User3; ///< User variables to prevent/delay future format changes 

  ClassDef(erhic::EventBeagle, 1)
};

}  // namespace erhic

#endif  // INCLUDE_EICSMEAR_ERHIC_EVENTPYTHIA_H_
