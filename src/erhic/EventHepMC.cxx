/**
 \file
 Implementation of class erhic::EventMilou.
 
 \author    Thomas Burton
 \date      2011-07-07
 \copyright 2011 Brookhaven National Lab
 */

#include "eicsmear/erhic/EventHepMC.h"

#include <cmath>
#include <sstream>
#include <string>
#include <iostream>

namespace erhic {

EventHepMC::EventHepMC()
: radcorr(-1)
, weight(NAN)
, trueX(NAN)
, trueQ2(NAN)
, trueY(NAN)
, trueT(NAN)
, truePhi(NAN)
, phibelgen(NAN)
, phibelres(NAN)
, phibelrec(NAN) {
}

}  // namespace erhic
