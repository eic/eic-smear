/**
 \file
 Implementation of class erhic::EventDis.
 
 \author    Thomas Burton 
 \date      2012-05-01
 \copyright 2012 Brookhaven National Lab
 */

#include "eicsmear/erhic/EventDis.h"

#include "eicsmear/erhic/Kinematics.h"

namespace erhic {

EventDis::~EventDis() { }

/**
 \todo Consider initialising float fields to zero, not NAN.
 e.g. in smeared output, particles outside acceptance retain default values,
 therefore in e.g. jacquet-blondel calculations, those not detected screw
 up the calculation because they have E, p = NAN, not zero!
 */
EventDis::EventDis()
: x(NAN)
, QSquared(NAN)
, y(NAN)
, WSquared(NAN)
, nu(NAN)
, yJB(NAN)
, QSquaredJB(NAN)
, xJB(NAN)
, WSquaredJB(NAN)
, yDA(NAN)
, QSquaredDA(NAN)
, xDA(NAN)
, WSquaredDA(NAN) {
}

EventDis::EventDis(const EventDis& that)
: VirtualEvent(that) {
  CopyKinematics(that);
}

EventDis& EventDis::operator=(const EventDis& that) {
  if (this != &that) {  // Protect against self-assignment
    CopyKinematics(that);
  }  // if
  return *this;
}

void EventDis::CopyKinematics(const EventDis& that) {
  SetLeptonKinematics(DisKinematics(that.x, that.y, that.nu,
                                    that.QSquared, that.WSquared));
  SetJacquetBlondelKinematics(DisKinematics(that.xJB, that.yJB, 0.,
                                            that.QSquaredJB, that.WSquaredJB));
  SetDoubleAngleKinematics(DisKinematics(that.xDA, that.yDA, 0.,
                                         that.QSquaredDA, that.WSquaredDA));
}

void EventDis::SetLeptonKinematics(const DisKinematics& kin) {
  x = kin.mX;
  QSquared = kin.mQ2;
  WSquared = kin.mW2;
  nu = kin.mNu;
  y = kin.mY;
}

void EventDis::SetJacquetBlondelKinematics(const DisKinematics& kin) {
  xJB = kin.mX;
  QSquaredJB = kin.mQ2;
  WSquaredJB = kin.mW2;
  yJB = kin.mY;
}

void EventDis::SetDoubleAngleKinematics(const DisKinematics& kin) {
  xDA = kin.mX;
  QSquaredDA = kin.mQ2;
  WSquaredDA = kin.mW2;
  yDA = kin.mY;
}

}  // namespace erhic
