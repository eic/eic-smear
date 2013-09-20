/**
 \file
 Declaration of class Smear::EventFactory.
 
 \author    Michael Savastio
 \date      2011-08-19
 \copyright 2011 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_SMEAR_EVENTFACTORY_H_
#define INCLUDE_EICSMEAR_SMEAR_EVENTFACTORY_H_

#include <string>

#include <TClass.h>
#include <TTree.h>
#include <TBranch.h>

#include "eicsmear/erhic/EventFactory.h"

namespace Smear {

/**
 Event factory for events of a particular type.
 Handles the Branch methods for that event type.
 The template argument can be any event type inheriting
 from erhic::VirtualEvent.
 */
template<typename T>
class EventFactory : public erhic::VirtualEventFactory {
 public:
  typedef T EventType;
  /**
   Destructor.
   */
  virtual ~EventFactory() { }

  virtual T* Create() = 0;

  /**
   Returns the name of the branch object type (Smear::Event).
   */
  virtual std::string EventName() const {
    return T::Class()->GetName();
  }

  /**
   Create and configure the output branch for smeared events in the tree.
   */
  virtual TBranch* Branch(TTree& tree, const std::string& name) {
    T* event(NULL);  // Temporary, just to pass to TTree::Branch()
    TBranch* branch = tree.Branch(name.c_str(), EventName().c_str(),
                                  &event, 32000, 99);
    branch->ResetAddress();
    if (event) {
      delete event;
    }  // if
    return branch;
  }

  /**
   Create() the new event and fill the requested tree branch with it.
   */
  virtual void Fill(TBranch& branch) {
    T* event = Create();
    branch.ResetAddress();
    branch.SetAddress(&event);
    branch.GetTree()->Fill();
    if (event) {
      delete event;
    }  // if
  }
};

}  // namespace Smear

#endif  // INCLUDE_EICSMEAR_SMEAR_EVENTFACTORY_H_
