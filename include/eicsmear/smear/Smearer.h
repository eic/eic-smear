/**
 Smearer.h

 \file
 Declaration of class Smearer.

 \author Thomas Burton 
 \date 3/28/12
 \copyright 2012 BNL. All rights reserved.
*/

#ifndef _EICSMEAR_SMEAR_SMEARER_H_
#define _EICSMEAR_SMEAR_SMEARER_H_

#include <TObject.h>

#include "eicsmear/smear/Acceptance.h"

namespace erhic {
   class VirtualParticle;
} // namespace erhic

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

      /** Destructor */
      virtual ~Smearer() { }

      /** Inherited from TObject */
      virtual Smearer* Clone(const char* = "") const = 0;

      /**
       Smears the input ParticleMC and stores the result(s) in the ParticleMCS
       */
      virtual void Smear(const erhic::VirtualParticle&, ParticleMCS&) = 0;
      
      Acceptance Accept;

      ClassDef(Smear::Smearer, 1)
   };

} // namespace Smear

#endif // _EICSMEAR_SMEAR_SMEARER_H_
