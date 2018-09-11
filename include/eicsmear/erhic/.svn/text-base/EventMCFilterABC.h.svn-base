/**
 \file
 Declaration of class erhic::EventMCFilterABC.
 
 \author    Thomas Burton
 \date      2012-01-18
 \copyright 2012 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_ERHIC_EVENTMCFILTERABC_H_
#define INCLUDE_EICSMEAR_ERHIC_EVENTMCFILTERABC_H_

#include <Rtypes.h>

namespace erhic {

class VirtualEvent;

/**
 Abstract base class for a filter on Monte Carlo events.
 */
class EventMCFilterABC {
 public:
  /**
   Constructor.
   */
  virtual ~EventMCFilterABC() { }

  /**
   Implement this method in a derived class such that it returns
   true when the event passes the filter's criteria.
   */
  virtual bool Accept(const VirtualEvent&) const = 0;

  ClassDef(erhic::EventMCFilterABC, 1)
};

}  // namespace erhic

#endif  // INCLUDE_EICSMEAR_ERHIC_EVENTMCFILTERABC_H_
