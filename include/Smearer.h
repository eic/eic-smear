/**
 Smearer.h

 \file
 Declaration of class Smearer.

 \author Thomas Burton 
 \date 3/28/12
 \copyright 2012 BNL. All rights reserved.
*/

#ifndef _Smearer_H_
#define _Smearer_H_

#include <TObject.h>

#include "Acceptance.h"

namespace erhic {
   class ParticleMC;
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
      virtual void DevSmear(const erhic::ParticleMC&, ParticleMCS&) = 0;
      
      Acceptance Accept;

      ClassDef(Smear::Smearer, 1)
   };

} // namespace Smear

#endif
