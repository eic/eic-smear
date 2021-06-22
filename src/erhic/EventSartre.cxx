/**
 \file Implementation of class erhic::EventSartre. 
 \author    Barak Schmookler
 \date      2020-04-08
 */

#include "eicsmear/erhic/EventSartre.h"

#include <cmath>
#include <sstream>
#include <string>

namespace erhic {

EventSartre::EventSartre()
: genevent(-1)
, trueT(NAN)
, trueQ2(NAN)
, trueX(NAN)
, trueY(NAN)
, trueW2(NAN)
, trueNu(NAN)
, trueXpom(NAN)
, s_cm(NAN)
, pol(-1)
, dmode(-1)
, bup(-1){
}

bool EventSartre::Parse(const std::string& line) {
  static std::stringstream ss;
  ss.str("");
  ss.clear();
  ss << line;
  ss >>
  number >> number >>  // Skip first int in the line
  genevent >>
  trueT >> trueQ2 >> trueX >> trueY >> trueW2 >> trueNu >>
  trueXpom >> s_cm >> pol >> dmode >> bup;
  // Protect against errors in the input file or the stream
  return !ss.fail();
}

  // Look for the scattered lepton in the event record.
  // This is the first (only?) particle that matches the following:
  //  1) pdg code equals that of incident lepton beam.
  //  2) status code is 1 i.e. it's a stable/final-state particle.
  //  3) the parent is track 1 or 2
  const ParticleMC* EventSartre::ScatteredLepton() const {
    // Look for the lepton beam to get the species.
    // If we don't get it we can't find the scattered
    // lepton so return NULL.
    const VirtualParticle* beam = BeamLepton();
    if (!beam) {
      return NULL;
    }  // if
    const int species = beam->Id().Code();
    // Get the final state particles and search them for
    // the scattered lepton.
    std::vector<const VirtualParticle*> final;
    FinalState(final);
    std::vector<const VirtualParticle*>::const_iterator iter;
    for (iter = final.begin(); iter != final.end(); ++iter) {
      // We already know the particle is final state, so
      // check its species and parent index.
      if ( (*iter)->Id().Code() == species && 
	   ( (*iter)->GetParentIndex() == 1 || (*iter)->GetParentIndex() == 2 )
	   ) {
	// Found it, cast to required particle type and return.
	return static_cast<const ParticleMC*>(*iter);
      }  // if
    }  // for
    // No luck, couldn't find the scattered lepton.
    return nullptr;
  }


  // Look for the exchange boson in the event record.
  // It would probably be the third track, but we'll go with the first status=21 boson
  //    that has particle 1 or 2 as parent
  const ParticleMC* EventSartre::ExchangeBoson() const {
    for ( auto particle : EventMC::GetTracks() ){
      if ( particle->GetStatus() != 21 ) continue;
      if ( particle->GetParentIndex() != 1 &&  particle->GetParentIndex() == 2 ) continue;
      if ( particle->Id().Code() == 22 ||  particle->Id().Code() == 23 ||  std::abs(particle->Id().Code()) == 24 ){
	// Found it, cast to required particle type and return.
	return static_cast<const ParticleMC*>(particle);
      }
    }  // for
    // No luck, couldn't find the exchange boson.
    return nullptr;
  }

}  // namespace erhic



