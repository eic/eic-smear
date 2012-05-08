#ifndef _EICSMEAR_EVENTBASE_
#define _EICSMEAR_EVENTBASE_

#include "eicsmear/erhic/EventMC.h"

// For backward compatibility with earlier code, where
// the event base class was named EventBase.
// EventBase is now a typedef to erhic::EventMC.
// Remove this typedef if the global definition ofEventBase clashes
// with something and you want just to use the namespace version.

typedef erhic::EventMC EventBase;

#endif
