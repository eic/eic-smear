/**
 \file
 Declaration of class erhic::EventDjangoh.
 
 \author    Thomas Burton
 \date      2011-07-07
 \copyright 2011 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_ERHIC_EVENTDJANGOH_H_
#define INCLUDE_EICSMEAR_ERHIC_EVENTDJANGOH_H_

#include <string>

#include <Rtypes.h>

#include "eicsmear/erhic/EventMC.h"

namespace erhic {

class ParticleMC;

/**
 Describes an event from the generator DJANGOH.
 */
class EventDjangoh : public EventMC {
 public:
  /**
   Constructor.
   */
  EventDjangoh() { }
  /**
   Parses the event information from a text string with the following format
   (no newlines):
   \verbatime
   "0 eventNum channel process subprocess nucleon parton partonTrack y Q2 x W2
   nu trueY tueQ2 trueX trueW2 trueNu crossSection crossSectionError
   depolarisation F1NC F3NC G1NC G3NC A1NC F1CC F3CC G1CC G5CC numTracks"
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

  Int_t nucleon;
  Int_t IChannel;
  Int_t dprocess;
  Int_t dstruckparton;
  Int_t dpartontrck;
  Double32_t dY;
  Double32_t dQ2;
  Double32_t dX;
  Double32_t dW2;
  Double32_t dNu;
  Double32_t dtrueY;
  Double32_t dtrueQ2;
  Double32_t dtrueX;
  Double32_t dtrueW2;
  Double32_t dtrueNu;
  Double32_t sigTot;
  Double32_t sigTotErr;
  Double32_t D;
  Double32_t F1NC;
  Double32_t F3NC;
  Double32_t G1NC;
  Double32_t G3NC;
  Double32_t A1NC;
  Double32_t F1CC;
  Double32_t F3CC;
  Double32_t G1CC;
  Double32_t G5CC;

  ClassDef(erhic::EventDjangoh, 2)
};

// DJANGOH gives particle output according to the LEPTO convention, whereby
// the exchange boson comes before the scattered lepton. This is different
// to the PYTHIA convention (lepton then boson), which is the default from
// EventMC.
inline const ParticleMC* EventDjangoh::ExchangeBoson() const {
  return GetTrack(2);
}

inline const ParticleMC* EventDjangoh::ScatteredLepton() const {
  return GetTrack(3);
}

}  // namespace erhic

#endif  // INCLUDE_EICSMEAR_ERHIC_EVENTDJANGOH_H_
