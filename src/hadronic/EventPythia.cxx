/**
 \file
 Implementation of class erhic::hadronic::EventPythiaPP.
 
 \author    Thomas Burton
 \date      2012-05-02
 \copyright 2012 Brookhaven National Lab
 */

#include <eicsmear/hadronic/EventPythia.h>

namespace erhic {
namespace hadronic {

EventPythiaPP::~EventPythiaPP() {
}

EventPythiaPP::EventPythiaPP()
: QSquared(0.)
, x1(0.)
, x2(0.) {
}

EventPythiaPP::EventPythiaPP(double Q2, double xa, double xb)
: QSquared(Q2)
, x1(xa)
, x2(xb) {
}

EventPythiaPP::EventPythiaPP(const EventPythiaPP& that)
: EventMC(that)
, QSquared(that.QSquared)
, x1(that.x1)
, x2(that.x2) {
}

EventPythiaPP& EventPythiaPP::operator=(const EventPythiaPP& that) {
  if (this != &that) {
    EventMC::operator=(that);
    QSquared = that.QSquared;
    x1 = that.x1;
    x2 = that.x2;
  }  // if
  return *this;
}

Double_t EventPythiaPP::GetQ2() const {
  return QSquared;
}

Double_t EventPythiaPP::GetX1() const {
  return x1;
}

Double_t EventPythiaPP::GetX2() const {
  return x2;
}

}  // namespace hadronic
}  // namespace erhic
