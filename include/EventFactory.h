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

#include "functions.h"
#include "VirtualEvent.h"
#include "EventPythia.h"
namespace erhic {
   
//   class VirtualEvent;
   class ParticleMC;

   // template<typename ParticleType>
   class VirtualEventFactory : public TObject {
      
   public:
      
      virtual ~VirtualEventFactory() { }
      
      virtual VirtualEvent<ParticleMC>* Create() = 0;
      
      
      ClassDef(VirtualEventFactory, 1)
   };
   

   /*
    class EventFromAsciiFactory : public VirtualEventFactory {
    
    public:
    
    virtual EventMC* Create() const = 0;
    };
    
   class PythiaEventFactory : public EventFromAsciiFactory {
      
   public:
      
      virtual 
   };
   
   
   class PepsiEventFactory : public EventMCFactor {
      
   public:
      
   };
   */
   
   
//   class EventMC;
   
   /**
    Creates events from an input plain text file containing
    appropriately formatted data.
    Templated for all the types inheriting from EventMC
    (any event class implementing a Parse method to
    populate the event's variables from a string)
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
