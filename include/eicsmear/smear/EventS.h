#ifndef _EICSMEAR_SMEAR_EVENTS_H_
#define _EICSMEAR_SMEAR_EVENTS_H_

// For backward compatibility with version 1
// of the code, where the EventS and ParticleS
// classes were both declared in EventS.h

#include "eicsmear/erhic/ParticleMC.h"
#include "eicsmear/smear/EventSmear.h"

//typedef Smear::Event::TrackType ParticleS;
typedef Smear::ParticleMCS ParticleS;

#endif // _EICSMEAR_SMEAR_EVENTS_H_
