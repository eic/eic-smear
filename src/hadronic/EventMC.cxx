#include <eicsmear/hadronic/EventMC.h>

#include <TLorentzVector.h>

namespace erhic {

namespace hadronic {

EventMC::~EventMC() {
   Clear("");
}

EventMC::EventMC()
: mTracks("erhic::hadronic::ParticleMC", 100) {
}

void EventMC::Clear(Option_t*) {
   
   mTracks.Clear();
}

const ParticleMC* EventMC::GetTrack(UInt_t i) const {
   return dynamic_cast<ParticleMC*>(mTracks.At(i));
}

ParticleMC* EventMC::GetTrack(UInt_t i) {
   return dynamic_cast<ParticleMC*>(mTracks.At(i));
}

UInt_t EventMC::GetNTracks() const {
   return mTracks.GetEntries();
}

UInt_t EventMC::Add(ParticleMC* p) {
   new(mTracks[GetNTracks()]) ParticleMC(*p);
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
