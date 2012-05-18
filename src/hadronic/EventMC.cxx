#include <eicsmear/hadronic/EventMC.h>

#include <TLorentzVector.h>

namespace erhic {
namespace hadronic {
   EventMC::~EventMC() {
      Clear("");
   }
   EventMC::EventMC() {
   }
   EventMC::EventMC(const EventMC& that)
   : VirtualEvent(that) {
      CopyTracks(that);
   }
   EventMC& EventMC::operator=(const EventMC& that) {
      if(this not_eq &that) {
         CopyTracks(that);
      } // if
      return *this;
   }
   void EventMC::Clear(Option_t*) {
      for(unsigned i(0); i < GetNTracks(); ++i) {
         delete mTracks.at(i);
         mTracks.at(i) = NULL;
      } // for
      mTracks.clear();
   }
   void EventMC::CopyTracks(const EventMC& that) {
      Clear("");
      for(unsigned i(0); i < that.GetNTracks(); ++i) {
         mTracks.push_back(
            dynamic_cast<ParticleMC*>(that.GetTrack(i)->Clone()));
      } // if
   }
   const ParticleMC* EventMC::GetTrack(UInt_t i) const {
      return mTracks.at(i);
   }
   ParticleMC* EventMC::GetTrack(UInt_t i) {
      return mTracks.at(i);
   }
   UInt_t EventMC::GetNTracks() const {
      return mTracks.size();
   }
   UInt_t EventMC::Add(ParticleMC* p) {
      if(p) {
         if(p->IsOnHeap()) {
            mTracks.push_back(p);
         } // if
      } // if
      return GetNTracks();
   }
   Double_t EventMC::GetCentreOfMassEnergy() const {
      double energy(NAN);
      if(GetTrack(0) and GetTrack(1)) {
         energy = (GetTrack(0)->Get4Vector() + GetTrack(1)->Get4Vector()).M();
      } // if
      return energy;
   }
} // namespace hadronic
} // namespace erhic
