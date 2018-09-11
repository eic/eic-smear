/**
 \file
 Declaration of class erhic::EventDis.
 
 \author    Thomas Burton 
 \date      2012-05-01
 \copyright 2012 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_ERHIC_EVENTDIS_H_
#define INCLUDE_EICSMEAR_ERHIC_EVENTDIS_H_

#include <cmath>  // For pow

#include <Rtypes.h>  // For ClassDef

#include "eicsmear/erhic/VirtualEvent.h"

namespace erhic {

struct DisKinematics;
class VirtualParticle;

/**
 A deeply inelastic scattering event.
 Stores kinematics computed by different methods:
 <ul>
 <li>Using the scattered electron</li>
 <li>Jacquet-Blondel method (uses hadrons)</li>
 <li>Double-angle method (uses both hadrons and the scattered electron)</li>
 </ul>
 This is an abstract class, as it does not implement track methods
 inherited from VirtualEvent.
 The user must implement the appropriate methods for their track type
 in an inheriting class.
 */
class EventDis : public VirtualEvent {
 public:
  /**
   Destructor.
   */
  virtual ~EventDis();

  /**
   Default constructor.
   */
  EventDis();

  /**
   Constructor copying another event's kinematics.
   */
  EventDis(const EventDis&);

  /**
   Assign another event's kinematics to this EventDis.
   */
  EventDis& operator=(const EventDis&);

  /**
   Returns Bjorken-x of the event.
   x<sub>B</sub> = Q<sup>2</sup>/(2p.q)
   */
  virtual Double_t GetX() const;

  /**
   Returns the four-momentum transfer (exchange boson mass) Q<sup>2</sup>.
   Q<sup>2</sup> = 2EE`(1+cos(theta)) = (e-e`)<sup>2</sup>
   */

  virtual Double_t GetQ2() const;
  /**
   Returns the event inelasticity.
   y = (p.q)/(p.e)
   */

  virtual Double_t GetY() const;
  /**
   Returns Y+ = y<sup>2</sup> / (1 + (1-y)<sup>2</sup>)
   */

  virtual Double_t GetYPlus() const;
  /**
   Returns the invariant mass of the hadronic final state.
   W<sup>2</sup> = M<sup>2</sup> + Q<sup>2</sup>(1-x)/x
   */

  virtual Double_t GetW2() const;
  /**
   Returns the exchange boson energy in the beam hadron rest frame.
   nu = q.p/M
   */

  virtual Double_t GetNu() const;
  /**
   Returns Bjorken x computed via the double-angle method.
   */
  virtual double GetXDoubleAngle() const;

  /**
   Returns Q-squared computed via the double-angle method.
   */
  virtual double GetQ2DoubleAngle() const;

  /**
   Returns inelasticity computed via the double-angle method.
   */
  virtual double GetYDoubleAngle() const;

  /**
   Returns W-squared computed via the double-angle method.
   */
  virtual double GetW2DoubleAngle() const;

  /**
   Returns Bjorken x computed via the Jacquet-Blondel method.
   */
  virtual double GetXJacquetBlondel() const;
  /**
   Returns Q-squared computed via the Jacquet-Blondel method.
    */
  virtual double GetQ2JacquetBlondel() const;

  /**
   Returns inelasticity computed via the Jacquet-Blondel method.
   */
  virtual double GetYJacquetBlondel() const;

  /**
   Returns W-squared computed via the Jacquet-Blondel method.
   */
  virtual double GetW2JacquetBlondel() const;

  /**
   Set the kinematics computed from the scattered lepton.
   */
  virtual void SetLeptonKinematics(const DisKinematics&);

  /**
   Set the kinematics computed from the Jacquet-Blondel method.
   */
  virtual void SetJacquetBlondelKinematics(const DisKinematics&);

  /**
   Set the kinematics computed from the double-angle method.
   */
  virtual void SetDoubleAngleKinematics(const DisKinematics&);

  /**
   Returns a pointer to the incident lepton beam particle.
   Returns NULL if the particle cannot be located in the event.
   IMPORTANT - DO NOT DELETE THE POINTER OR BAD THINGS WILL HAPPEN!
   */
  virtual const VirtualParticle* BeamLepton() const = 0;

  /**
   Returns a pointer to the incident hadron beam particle.
   Returns NULL if the particle cannot be located in the event.
   IMPORTANT - DO NOT DELETE THE POINTER OR BAD THINGS WILL HAPPEN!
   */
  virtual const VirtualParticle* BeamHadron() const = 0;

  /**
   Returns a pointer to the exchanged boson.
   Returns NULL if the particle cannot be located in the event.
   IMPORTANT - DO NOT DELETE THE POINTER OR BAD THINGS WILL HAPPEN!
   */
  virtual const VirtualParticle* ExchangeBoson() const = 0;

  /**
   Returns a pointer to the lepton beam particle after scattering.
   Returns NULL if the particle cannot be located in the event.
   IMPORTANT - DO NOT DELETE THE POINTER OR BAD THINGS WILL HAPPEN!
   */
  virtual const VirtualParticle* ScatteredLepton() const = 0;

// protected:
  /**
   Copy the kinematics from another event to this event.
   */
  virtual void CopyKinematics(const EventDis&);
  Double32_t x;  ///< Bjorken scaling variable
  Double32_t QSquared;  ///< Q<sup>2</sup> calculated from scattered electron
  Double32_t y;  ///< Inelasticity
  Double32_t WSquared;  ///< Invariant mass of the hadronic system
  Double32_t nu;  ///< Energy transfer from the electron
  Double32_t yJB;  ///< y calculated via the Jacquet-Blondel method
  Double32_t QSquaredJB;  ///< Q2 calculated via the Jacquet-Blondel method
  Double32_t xJB;  ///< x calculated via the Jacquet-Blondel method
  Double32_t WSquaredJB;  ///< W2 calculated via the Jacquet-Blondel method
  Double32_t yDA;  ///< y calculated via the double-angle method
  Double32_t QSquaredDA;  ///< Q2 calculated via the double-angle method
  Double32_t xDA;  ///< x calculated via the double-angle method
  Double32_t WSquaredDA;  ///< W2 calculated via the double-angle method

  ClassDef(erhic::EventDis, 1)
};

inline Double_t EventDis::GetX() const {
  return x;
}

inline Double_t EventDis::GetNu() const {
  return nu;
}

inline Double_t EventDis::GetQ2() const {
  return QSquared;
}

inline Double_t EventDis::GetW2() const {
  return WSquared;
}

inline Double_t EventDis::GetY() const {
  return y;
}

inline Double_t EventDis::GetYPlus() const {
  return pow(GetY(), 2.) / (1. + pow(1. - GetY(), 2.));
}

inline double EventDis::GetXDoubleAngle() const {
  return xDA;
}

inline double EventDis::GetQ2DoubleAngle() const {
  return QSquaredDA;
}

inline double EventDis::GetYDoubleAngle() const {
  return yDA;
}

inline double EventDis::GetW2DoubleAngle() const {
  return WSquaredDA;
}

inline double EventDis::GetXJacquetBlondel() const {
  return xJB;
}

inline double EventDis::GetQ2JacquetBlondel() const {
  return QSquaredJB;
}

inline double EventDis::GetYJacquetBlondel() const {
  return yJB;
}

inline double EventDis::GetW2JacquetBlondel() const {
  return WSquaredJB;
}

}  // namespace erhic

#endif  // INCLUDE_EICSMEAR_ERHIC_EVENTDIS_H_
