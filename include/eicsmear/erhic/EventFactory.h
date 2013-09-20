/**
 \file
 Declaration of class erhic::EventFactory.
 
 \author    Thomas Burton
 \date      2011-10-31
 \copyright 2011 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_ERHIC_EVENTFACTORY_H_
#define INCLUDE_EICSMEAR_ERHIC_EVENTFACTORY_H_

#include <iostream>
#include <memory>
#include <string>

#include <TBranch.h>
#include <TTree.h>

#include "eicsmear/functions.h"
#include "eicsmear/erhic/VirtualEvent.h"
#include "eicsmear/erhic/EventPythia.h"
#include "eicsmear/smear/EventSmear.h"

namespace erhic {

class ParticleMC;

/**
 Abstract base class for event builders.
 Due to how ROOT handles allocating branch objects,
 we take care of creating and filling branches in this
 class also.
 */
class VirtualEventFactory : public TObject {
 public:
  /**
   Destructor.
   */
  virtual ~VirtualEventFactory() { }

  /**
   Returns a new event instance.
   */
  virtual VirtualEvent* Create() = 0;

  /**
   Returns a pointer to the event buffer.
   */
  virtual VirtualEvent* GetEvBufferPtr() { return 0; }

  /**
   Returns a string with the full (including namespace) class name
   of the event type produced.
   This is important for use with ROOT TTree to ensure the correct
   event type in branches.
   */
  virtual std::string EventName() const = 0;

  /**
   Add a branch named "name" for the event type generated
   by this factory to a ROOT TTree.
   Returns a pointer to the branch, or NULL in the case of an error.
   */
  virtual TBranch* Branch(TTree&, const std::string&) {
    return NULL;
  }
  /**
   Calls Create() to generate an event and fills the provided
   branch with that event.
   Resets the branch address in doing so.
   */
  virtual void Fill(TBranch&) { }

  ClassDef(VirtualEventFactory, 1)
};

/**
 Creates events from an input plain text file containing
 appropriately formatted data.
 Templated for all the types inheriting from EventMC
 (any event class implementing a Parse() method to
 populate the event's variables from a string will work.)
 */
template<typename T>
class EventFromAsciiFactory : public VirtualEventFactory {
 public:
  /**
   Constructor.
   */
  EventFromAsciiFactory() { }

  /**
   Destructor.
   */
  virtual ~EventFromAsciiFactory() { }

  /**
   Initialise the factory from an input stream. 
   */
  explicit EventFromAsciiFactory(std::istream& is)
  : mInput(&is)
  , mEvent(NULL) {
  }

  /**
   Returns a new event instance.
   */
  virtual T* Create();

  /**
   Returns the name of the event class created by this factory.
   */
  virtual std::string EventName() const;

  std::istream* mInput;  //!
  std::string mLine;  //!
  std::auto_ptr<T> mEvent;  //!

 protected:
  /**
   Returns true when an end-of-event marker is encountered in the input stream.
   */
  bool AtEndOfEvent() const;

  /**
   Perform end-of-event operations.
   */
  Int_t FinishEvent();

  /**
   Create a new particle from the last data read from the input stream.
   */
  bool AddParticle();

  // Warning: explicitly putting the erhic:: namespace before the class
  // name doesn't seen to work for template classes.
  ClassDef(EventFromAsciiFactory, 1)
};

}  // namespace erhic

#endif  // INCLUDE_EICSMEAR_ERHIC_EVENTFACTORY_H_
