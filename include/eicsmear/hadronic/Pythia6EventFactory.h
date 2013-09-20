/**
 \file
 Declaration of class erhic::hadronic::Pythia6EventFactory.
 
 \author    Thomas Burton 
 \date      2012-05-03
 \copyright 2012 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_HADRONIC_PYTHIA6EVENTFACTORY_H_
#define INCLUDE_EICSMEAR_HADRONIC_PYTHIA6EVENTFACTORY_H_

#include <memory>
#include <string>

#include <Rtypes.h>  // For ClassDef macro

#include <eicsmear/hadronic/EventPythia.h>
#include <eicsmear/erhic/EventFactory.h>
#include <eicsmear/erhic/EventMCFilterABC.h>

class TBranch;

namespace erhic {
namespace hadronic {

/**
 A realisation of erhic::VirtualEventFactory generating
 objects of type hadronic::EventPythiaPP.
 */
class Pythia6EventFactory : public erhic::VirtualEventFactory {
 public:
  /**
   Destructor.
   */
  virtual ~Pythia6EventFactory();

  /**
   Constructor.
   */
  explicit Pythia6EventFactory(EventMCFilterABC* filter);

  /**
   Returns a new event instance.
   */
  virtual EventPythiaPP* Create();

  /**
   Returns the name of the event class created by this factory.
   */
  virtual std::string EventName() const;

  /**
   Create a new branch on a tree.
   */
  virtual TBranch* Branch(TTree& tree, const std::string& branchName);

  /**
   Fill a tree branch with the current event.
   */
  virtual void Fill(TBranch&);

 protected:
  /**
   Construct an event from the current state of TPythia6.
   */
  virtual EventPythiaPP* BuildEvent();

  std::auto_ptr<erhic::EventMCFilterABC> mFilter;
  EventPythiaPP* mEvent;

  ClassDef(erhic::hadronic::Pythia6EventFactory, 1)
};

}  // namespace hadronic
}  // namespace erhic

#endif  // INCLUDE_EICSMEAR_HADRONIC_PYTHIA6EVENTFACTORY_H_
