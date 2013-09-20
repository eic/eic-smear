/**
 \file
 For backward compatibility with earlier code, where
 the event base class was named EventBase.
 
 \author    Thomas Burton
 \date      2011-07-07
 \copyright 2011 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_ERHIC_EVENTBASE_H_
#define INCLUDE_EICSMEAR_ERHIC_EVENTBASE_H_

#include "eicsmear/erhic/EventMC.h"

// For backward compatibility with earlier code, where
// the event base class was named EventBase.
// EventBase is now a typedef to erhic::EventMC.
// Remove this typedef if the global definition of EventBase clashes
// with something and you want just to use the namespace version.

typedef erhic::EventMC EventBase;

#endif  // INCLUDE_EICSMEAR_ERHIC_EVENTBASE_H_
