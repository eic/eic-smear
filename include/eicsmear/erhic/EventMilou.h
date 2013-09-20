/**
 \file
 Declaration of class erhic::EventMilou.
 
 \author    Thomas Burton
 \date      2011-07-07
 \copyright 2011 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_ERHIC_EVENTMILOU_H_
#define INCLUDE_EICSMEAR_ERHIC_EVENTMILOU_H_

#include <string>

#include <Rtypes.h>

#include "eicsmear/erhic/EventMC.h"

namespace erhic {

/**
 Describes an event from the generator MILOU.
 */
class EventMilou : public EventMC {
 public:
  /**
   Constructor.
   */
  EventMilou();

  /**
   Parses the event information from a text string.
   
   The string must have the following format (no newlines):
   \verbatim
   "0 eventnumber numTracks weight processId radiativeCorrectionFlag
   trueX trueQ2 trueY trueT truePhi phi phiResolution reconstructedPhi"
   \endverbatim
   Returns true in the event of a successful read operation,
   false in case of an error.
   */
  virtual bool Parse(const std::string&);

  /**
   Azimuthal angle between the production and the scattering plane.
   */
  Double_t GetPhiBelGen() const;

  /**
   Resolution in azimuthal angle.
   */
  Double_t GetPhiBelRes() const;

  /**
   Reconstructed azimuthal angle.
   */
  Double_t GetPhiBelRec() const;

  Bool_t radcorr;
  Double32_t weight;
  Double32_t trueX;
  Double32_t trueQ2;
  Double32_t trueY;
  Double32_t trueT;
  Double32_t truePhi;
  Double32_t phibelgen;  ///> the azimuthal angle between the production
                         ///> and the scattering plane
  Double32_t phibelres;  ///> the resolution of the previous angle
                         ///> according to H1
  Double32_t phibelrec;  ///> the reconstructed phi angle

  ClassDef(erhic::EventMilou, 1)
};

inline Double_t EventMilou::GetPhiBelGen() const {
  return phibelgen;
}

inline Double_t EventMilou::GetPhiBelRes() const {
  return phibelres;
}

inline Double_t EventMilou::GetPhiBelRec() const {
  return phibelrec;
}

}  // namespace erhic

#endif  // INCLUDE_EICSMEAR_ERHIC_EVENTMILOU_H_
