#ifndef _ERHIC_BUILDTREE_EVENTBASE_
#define _ERHIC_BUILDTREE_EVENTBASE_

#include "EventMC.h"

// For backward compatibility with earlier code.
// EventBase is a typedef to erhic::EventMC.
// Remove this typedef if the global EventBase clashes with something and you
// want just to use the namespace version.

typedef erhic::EventMC EventBase;

#endif
