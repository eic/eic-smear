/**
 \file
 Declaration of class Smear::AcceSmearerptance.
 
 \author    Michael Savastio
 \date      2011-08-19
 \copyright 2011 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_SMEAR_SMEARER_H_
#define INCLUDE_EICSMEAR_SMEAR_SMEARER_H_

#include <TObject.h>

#include "eicsmear/smear/Acceptance.h"

namespace erhic {

class VirtualParticle;

}  // namespace erhic

namespace Smear {

class ParticleMCS;

/**
 Abstract base class for objects performing smearing.
 Has no constructors or assignemt operator, as it is an abstract class.
 \todo Data hiding for Acceptance, but need to address more general data
 hiding issues in smearing code first.
 */
class Smearer : public TObject {
 public:
  /**
   Destructor.
   */
  virtual ~Smearer() { }

  /**
   Inherited from TObject.
   */
  virtual Smearer* Clone(const char* = "") const = 0;

  /**
   Smears the input ParticleMC and stores the result(s) in the ParticleMCS.
   */
  virtual void Smear(const erhic::VirtualParticle&, ParticleMCS&) = 0;

  Acceptance Accept;

  ClassDef(Smear::Smearer, 1)
};

}  // namespace Smear

#endif  // INCLUDE_EICSMEAR_SMEAR_SMEARER_H_
