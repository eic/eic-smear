#ifndef INCLUDE_EICSMEAR_ERHIC_EVENTFACTORYHEPMC_H_
#define INCLUDE_EICSMEAR_ERHIC_EVENTFACTORYHEPMC_H_

#include "eicsmear/erhic/EventFactory.h"

#include <HepMC3/GenEvent.h>
#include <HepMC3/GenVertex.h>
#include <HepMC3/GenParticle.h>

#include<map>

namespace erhic {

  /**
     Helper. This clumsy construction would be much better handled with private members, but
     it doesn't work like that if we want to specialize from EventFromAsciiFactory while keeping some functions.
   */
  void HandleHepmcParticle( const HepMC3::GenParticlePtr& p, std::map < HepMC3::GenParticlePtr, int >& hepmcp_index, int& particleindex, std::unique_ptr<erhic::EventHepMC>& mEvent );

}  // namespace erhic

#endif  // INCLUDE_EICSMEAR_ERHIC_EVENTFACTORY_H_
