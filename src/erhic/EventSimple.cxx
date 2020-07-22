/**
 \file Implementation of class erhic::EventSimple. 
 \author    Barak Schmookler
 \date      2020-03-19
 */

#include "eicsmear/erhic/EventSimple.h"

#include <cmath>
#include <sstream>
#include <string>

namespace erhic {

EventSimple::EventSimple()
: numParticles(NAN){
}

bool EventSimple::Parse(const std::string& line) {
  static std::stringstream ss;
  ss.str("");
  ss.clear();
  ss << line;
  ss >>
  number >> number >>  // Skip first int in the line
    numParticles;
  // Protect against errors in the input file or the stream
  return !ss.fail();
}

// Look for the scattered lepton in the event record.
// This is the first (only?) particle that matches the following:
//  1) pdg code equals that of incident lepton beam.
//  2) status code is 1 i.e. it's a stable/final-state particle.
const ParticleMC* EventSimple::ScatteredLepton() const {
  // Look for the lepton beam to get the species.
  // If we don't get it we can't find the scattered
  // lepton so return NULL.
  const VirtualParticle* beam = BeamLepton();
  if (!beam) {
    return nullptr;
  }  // if
  const int species = beam->Id().Code();
  // Get the final state particles and search them for
  // the scattered lepton.
  std::vector<const VirtualParticle*> final;
  FinalState(final);
  // std::vector<const VirtualParticle*>::const_iterator iter;
  for (auto & iter : final ) {
    // We already know the particle is final state
    // could check for its parent but we don't need to
    if ( iter->Id().Code() == species) {
      // Found it, cast to required particle type and return.
      return static_cast<const ParticleMC*>(iter);
    }
  }

  // No luck, couldn't find the scattered lepton.
  return nullptr;
}

}  // namespace erhic
