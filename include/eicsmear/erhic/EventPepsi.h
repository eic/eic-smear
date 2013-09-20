/**
 \file
 Declaration of class erhic::EventPepsi.
 
 \author    Thomas Burton
 \date      2011-07-07
 \copyright 2011 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_ERHIC_EVENTPEPSI_H_
#define INCLUDE_EICSMEAR_ERHIC_EVENTPEPSI_H_

#include <string>

#include <Rtypes.h>

#include "eicsmear/erhic/EventMC.h"

namespace erhic {

/**
 Describes an event from the generator PEPSI.

 \todo Add accessor and setter methods
 */
class EventPepsi : public EventMC {
 public:
  /**
   Parses the event information from a text string with the following format
   (no newlines):
   \verbatime
   "I ievent genevent process subprocess nucleon struckparton,
   partontrck trueY trueQ2 trueX trueW2 trueNu FixedWeight,
   weight dxsec Extraweight dilut F1 F2 A1 A2 R Depol d,
   eta eps chi gendilut genF1 genF2 genA1 genA2 genR genDepol,
   gend geneta geneps genchi Sigcor radgamEnucl nrTracks"
   \endverbatim
   Returns true in the event of a successful read operation,
   false in case of an error.
   */
  virtual bool Parse(const std::string&);

  /**
   Returns a pointer to the exchange boson, or NULL if it cannot be found.
   */
  virtual const ParticleMC* ExchangeBoson() const;

  /**
   Returns a pointer to the scattered lepton, or NULL if it cannot be found.
   */
  virtual const ParticleMC* ScatteredLepton() const;

  Int_t nucleon;  ///< PDG code of the hadron beam
  Int_t struckparton;  ///< Parton hit in the target LST(25)
  Int_t partontrck;  ///< Number of parton track LST(26)
  Int_t genevent;  ///< Trials required for this event
  Int_t subprocess;  ///< PEPSI subprocess LST(23)
  Double32_t trueY;  ///< Generated y of the event
  Double32_t trueQ2;  ///< Generated Q<sup>2</sup> of the event
  Double32_t trueX;  ///< Generated x of the event
  Double32_t trueW2;  ///< Generated W<sup>2</sup> of the event
  Double32_t trueNu;  ///< Generated nu of the event
  Double32_t FixedWeight;  ///< Weight calculated from generation limits
  Double32_t Weight;  ///< Total weight including everything
  Double32_t dxsec;  ///< Cross section included in the weight
  Double32_t ExtraWeight;  ///< PEPSI total cross section in pb from
                           ///< numerical integration PARL(23)
  Double32_t Dilute;  ///< True variables needed to calculate g<sub>1</sub>
  Double32_t F1;  ///< True variables needed to calculate g<sub>1</sub>
  Double32_t F2;  ///< True variables needed to calculate g<sub>1</sub>
  Double32_t A1;  ///< True variables needed to calculate g<sub>1</sub>
  Double32_t A2;  ///< True variables needed to calculate g<sub>1</sub>
  Double32_t R;  ///< True variables needed to calculate g<sub>1</sub>
  Double32_t DePol;  ///< True variables needed to calculate g<sub>1</sub>
  Double32_t D;  ///< True variables needed to calculate g<sub>1</sub>
  Double32_t Eta;  ///< True variables needed to calculate g<sub>1</sub>
  Double32_t Eps;  ///< True variables needed to calculate g<sub>1</sub>
  Double32_t Chi;  ///< True variables needed to calculate g<sub>1</sub>
  Double32_t gendilut;  ///< Needed to calculate g1
  Double32_t genF1;  ///< Needed to calculate g1
  Double32_t genF2;  ///< Needed to calculate g1
  Double32_t genA1;  ///< Needed to calculate g1
  Double32_t genA2;  ///< Needed to calculate g1
  Double32_t genR;  ///< Needed to calculate g1
  Double32_t genDepol;  ///< Needed to calculate g1
  Double32_t gend;  ///< Needed to calculate g1
  Double32_t geneta;  ///< Needed to calculate g1
  Double32_t geneps;  ///< Needed to calculate g1
  Double32_t genchi;  ///< Needed to calculate g1
  Double32_t SigCorr;  ///< Needed in the radiative correction code
  Double32_t radgamEnucl;

  ClassDef(erhic::EventPepsi, 1)
};

// DJANGOH gives particle output according to the LEPTO convention, whereby
// the exchange boson comes before the scattered lepton. This is different
// to the PYTHIA convention (lepton then boson), which is the default from
// EventMC.
inline const ParticleMC* EventPepsi::ExchangeBoson() const {
  return GetTrack(2);
}

inline const ParticleMC* EventPepsi::ScatteredLepton() const {
  return GetTrack(3);
}

}  // namespace erhic

#endif  // INCLUDE_EICSMEAR_ERHIC_EVENTPEPSI_H_
