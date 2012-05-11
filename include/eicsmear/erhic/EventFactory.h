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

namespace Smear {
   inline ParticleMCS* mcToSmear(const erhic::VirtualParticle& mc) {
      return new ParticleMCS(mc.Get4Vector(), mc.Id(), mc.GetStatus());
   }
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
      /** Destructor */
      virtual ~EventFactory() {
      }
      virtual T* Create() = 0;
      /** Returns the name of the branch object type (Smear::Event) */
      virtual std::string EventName() const { 
         return T::Class()->GetName();
      }
      /**
       Create and configure the output branch for smeared events in the tree.
      */
      virtual TBranch* Branch(TTree& tree, const std::string& name) {
         T* event(NULL); // Temporary, just to pass to TTree::Branch()
         TBranch* branch = tree.Branch(name.c_str(), EventName().c_str(),
                                       &event, 32000, 99);
         branch->ResetAddress();
         if(event) {
            delete event;
         } // if
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
         if(event) {
            delete event;
         } // if
      }
   };

   /**
    Factory class for smeared DIS events.
   */
   class EventDisBuilder : public EventFactory<Smear::Event> {
   public:
      /** Destructor */
      virtual ~EventDisBuilder() {
      }
      /**
       Constructor.
       Initialise with the Detector performing particle smearing and
       the tree branch providing DIS Monte Carlo events.
      */
      EventDisBuilder(const Detector& d, TBranch& mcBranch)
      : mDetector(d)
      , mMcEvent(NULL) {
         mcBranch.SetAddress(&mMcEvent);
      }
      /**
       Create a smeared event corresponding to the current DIS Monte Carlo
       event in the input branch passed to the constructor.
       The user should call TTree::GetEntry() on the tree with the
       input branch between calls to Create().
      */
      virtual Event* Create() {
         Event* event = new Event;
         ParticleIdentifier pid(mMcEvent->BeamLepton()->Id());
         for(unsigned j(0); j < mMcEvent->GetNTracks(); j++) {
            const erhic::VirtualParticle* ptr = mMcEvent->GetTrack(j);
            if(not ptr) {
               continue;
            } // if
            // It's convenient to keep the initial beams, unsmeared, in the
            // smeared event record, so copy their properties exactly
            if(mMcEvent->BeamLepton() == ptr or mMcEvent->BeamHadron() == ptr) {
               event->AddLast(mcToSmear(*ptr));
            } // if
            else {
               event->AddLast(mDetector.Smear(*ptr));
               // If this is the scattered lepton, record the index.
               // Check that the scattered lepton hasn't been set yet so we
               // don't replace it with a subsequent match.
               if(pid.isScatteredLepton(*ptr) and not event->ScatteredLepton()) {
                  event->mScatteredIndex = j;
               } // if
            } // else
         } // for
         // Fill the event-wise kinematic variables.
         mDetector.FillEventKinematics(*mMcEvent, event);
         return event;
      }
   protected:
      Detector mDetector;
      erhic::EventDis* mMcEvent;
   };
}

#endif
