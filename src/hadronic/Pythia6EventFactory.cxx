#include <eicsmear/hadronic/Pythia6EventFactory.h>

#include <memory>

#include <TBranch.h>
#include <TClass.h>
#include <TCollection.h> // For TIter
#include <TMCParticle.h>
#include <TObjArray.h>
#include <TPythia6.h>
#include <TTree.h>

namespace erhic {
namespace hadronic {
   Pythia6EventFactory::~Pythia6EventFactory() {
   }
   Pythia6EventFactory::Pythia6EventFactory() {
   }
   EventPythiaPP* Pythia6EventFactory::Create() {
      // Read the event kinematics from PYTHIA and create an event
      TPythia6* pythia = TPythia6::Instance();
      double Q2 = pythia->GetPARI(22);
      double x1 = pythia->GetPARI(33);
      double x2 = pythia->GetPARI(34);
      std::auto_ptr<EventPythiaPP> event(new EventPythiaPP(Q2, x1, x2));
      // Get the particles from the current PYTHIA event.
      // Build a ParticleMC from each and add to the event's list
      TObjArray* particles = pythia->ImportParticles("All");
      TIter iter(particles);
      TMCParticle* mc(NULL);
      while((mc = dynamic_cast<TMCParticle*>(iter.Next()))) {
         if(mc) {
            ParticleMC* p = new ParticleMC(*mc);
            p->SetParentIndex(mc->GetParent());
//            p->SetXFeynman(2. * p->GetPz() / event->GetCentreOfMassEnergy());
            event->Add(p);
         } // if
      } // while
      return event.release();
   }
   std::string Pythia6EventFactory::EventName() const {
      return EventPythiaPP::Class()->GetName();
   }
   TBranch* Pythia6EventFactory::Branch(TTree& tree, const std::string& name) {
      EventPythiaPP* event(NULL);
      std::cout << EventName() << std::endl;
      TBranch* branch =
         tree.Branch(name.c_str(), EventName().c_str(), &event, 32000, 99);
      tree.ResetBranchAddress(branch);
      return branch;
   }
   void Pythia6EventFactory::Fill(TBranch& branch) {
      EventPythiaPP* event(NULL);
      branch.SetAddress(&event);
      event = Create();
      branch.GetTree()->Fill();
      branch.ResetAddress();
      delete event;
   }
} // namespace hadronic
} // namespace erhic
