/**
 \file
 Implementation of class erhic::hadronic::Pythia6EventFactory.
 
 \author    Thomas Burton 
 \date      2012-05-03
 \copyright 2012 Brookhaven National Lab
 */

#include <eicsmear/hadronic/Pythia6EventFactory.h>

#include <iostream>
#include <memory>
#include <string>

#include <TBranch.h>
#include <TClass.h>
#include <TCollection.h>  // For TIter
#include <TMCParticle.h>
#include <TObjArray.h>
#include <TProcessID.h>
#include <TPythia6.h>
#include <TTree.h>

namespace erhic {
namespace hadronic {

Pythia6EventFactory::~Pythia6EventFactory() {
}

Pythia6EventFactory::Pythia6EventFactory(erhic::EventMCFilterABC* filter)
: mFilter(filter)
, mEvent(NULL) {
}

EventPythiaPP* Pythia6EventFactory::Create() {
  std::auto_ptr<EventPythiaPP> event(BuildEvent());
  if (mFilter.get()) {
    while (!mFilter->Accept(*event)) {
      event.reset(BuildEvent());
    }  // while
  }  // if
  return event.release();
}

EventPythiaPP* Pythia6EventFactory::BuildEvent() {
  // Save current object count
  int objectNumber = TProcessID::GetObjectCount();
  // Generate a new PYTHIA event
  // Read the event kinematics from PYTHIA and create an event
  TPythia6* pythia = TPythia6::Instance();
  pythia->GenerateEvent();
  double Q2 = pythia->GetPARI(22);
  double x1 = pythia->GetPARI(33);
  double x2 = pythia->GetPARI(34);
  std::auto_ptr<EventPythiaPP> event(new EventPythiaPP(Q2, x1, x2));
  // Get the particles from the current PYTHIA event.
  // Build a ParticleMC from each and add to the event's list
  TObjArray* particles = pythia->ImportParticles("All");
  TIter iter(particles);
  TMCParticle* mc(NULL);
  // Populate particle list
  while ((mc = static_cast<TMCParticle*>(iter.Next()))) {
    if (mc) {
      std::auto_ptr<ParticleMC> p(new ParticleMC(*mc));
      p->SetParentIndex(mc->GetParent());
      event->Add(p.get());
    }  // if
  }  // while
  // Test against filter and exit loop if the event passes the filter.
  // Restore Object count
  // See example in $ROOTSYS/test/Event.cxx
  // To save space in the table keeping track of all referenced objects
  // we assume that our events do not address each other. We reset the
  // object count to what it was at the beginning of the event.
  TProcessID::SetObjectCount(objectNumber);
  return event.release();
}

std::string Pythia6EventFactory::EventName() const {
  return EventPythiaPP::Class()->GetName();
}

TBranch* Pythia6EventFactory::Branch(TTree& tree, const std::string& name) {
  EventPythiaPP* event(NULL);
  TBranch* branch =
  tree.Branch(name.c_str(), EventName().c_str(), &event, 32000, 99);
  tree.ResetBranchAddress(branch);
  return branch;
}

void Pythia6EventFactory::Fill(TBranch& branch) {
  if (mEvent) {
    branch.ResetAddress();
    delete mEvent;
    mEvent = NULL;
  }  // if
  mEvent = Create();
  branch.SetAddress(&mEvent);
}

}  // namespace hadronic
}  // namespace erhic
