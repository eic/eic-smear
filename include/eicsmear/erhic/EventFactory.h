#ifndef _EventFactory_H_
#define _EventFactory_H_

//
// EventFactory.h
//
// Created by TB on 10/31/11.
// Copyright 2011 BNL. All rights reserved.
//

#include <iostream>
#include <memory>
#include <string>

#include "eicsmear/functions.h"
#include "eicsmear/erhic/VirtualEvent.h"
#include "eicsmear/erhic/EventPythia.h"
#include "eicsmear/smear/EventSmear.h"

//class TBranch;
//class TTree;
#include <TBranch.h>
#include <TTree.h>
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
      
      virtual ~VirtualEventFactory() { }
      
      virtual VirtualEvent* Create() = 0;
      
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
      
      EventFromAsciiFactory() { }
      
      EventFromAsciiFactory(std::istream& is)
      : mInput(&is)
      , mEvent(NULL) {
      }
      
      virtual T* Create();
      
      virtual ~EventFromAsciiFactory() { }
      
      virtual std::string EventName() const;
      
      std::istream* mInput; //!
      std::string mLine; //!
      std::auto_ptr<T> mEvent; //!
      
   protected:
      
      bool AtEndOfEvent() const;
      
      Int_t FinishEvent();
      
      bool AddParticle();
      
      ClassDef(EventFromAsciiFactory, 1)
   };
   
   
} // namespace erhic

#endif
