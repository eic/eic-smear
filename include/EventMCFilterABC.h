/**
 EventMCFilterABC.h

 \file
 Declaration of class EventMCFilterABC.

 Created by TB on 1/18/12.
 Copyright 2012 BNL. All rights reserved.
*/

#ifndef _EventMCFilterABC_H_
#define _EventMCFilterABC_H_

#include <Rtypes.h>

namespace erhic {
   
   class EventMC;
   
   /** Abstract base class for a filter on Monte Carlo events. */
   
   class EventMCFilterABC {
      
   public:
      
      virtual ~EventMCFilterABC() { }
      
      /**
       Implement this method in a derived class such that it returns
       true when the event passes the filter's criteria.
       */
      virtual bool Accept(const EventMC&) const = 0;
      
      ClassDef(erhic::EventMCFilterABC, 1)
   };
   
} // namespace erhic

#endif
