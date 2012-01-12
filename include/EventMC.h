//
// EventMC.h
// BuildTree
//
// Created by TB on 10/10/11.
// Copyright 2011 BNL. All rights reserved.
//

#ifndef _ERHIC_EVENTMC_H_
#define _ERHIC_EVENTMC_H_

#include <cmath>
#include <list>

#include <TLorentzVector.h>
#include <TObject.h>

#include "VirtualEvent.h"
#include "VirtualParticle.h"

class TTree;

//#ifdef USE_NAMESPACE_ERHIC
namespace erhic {
   
   class ParticleMC;
   
   /**
    Abstract base class for Monte Carlo events.
    Implements common event properties and methods.
    */
   class EventMC : public VirtualEvent<ParticleMC> {
      
   public:
      
      /**
       Constructor.
       */
      EventMC();
      
      /**
       Destructor.
       */
      virtual ~EventMC();
      
      /**
       Returns the index of this event.
       */
      virtual ULong64_t N() const;
      
      /**
       Returns Bjorken-x of the event.
       x<sub>B</sub> = Q<sup>2</sup>/(2p.q)
       */
      virtual Double_t X() const { return x; }
      
      /**
       Returns the four-momentum transfer (exchange boson mass) Q<sup>2</sup>.
       Q<sup>2</sup> = 2EE`(1+cos(theta)) = (e-e`)<sup>2</sup>
       */
      virtual Double_t Q2() const { return QSquared; }
      
      /**
       Returns the event inelasticity.
       y = (p.q)/(p.e)
       */
      virtual Double_t Y() const { return y; }
      
      /**
       Returns the invariant mass of the hadronic final state.
       W<sup>2</sup> = M<sup>2</sup> + Q<sup>2</sup>(1-x)/x
       */
      virtual Double_t W2() const { return WSquared; }
      
      /**
       Returns the exchange boson energy in the beam hadron rest frame.
       nu = q.p/M
       */
      virtual Double_t Nu() const { return nu; }
      
      /**
       Returns a code describing the production process of this event.
       */
      virtual Int_t Process() const;
      
      /**
       Returns the number of tracks in the event.
       */
      virtual UInt_t NTracks() const;
      
      /**
       Returns the nth track.
       Returns NULL if the track number is out of range.
       */
      virtual const ParticleMC* GetTrack(UInt_t) const;
      
      /**
       Returns the nth track.
       Returns NULL if the track number is out of range.
       */
      virtual ParticleMC* GetTrack(UInt_t);
      
      /**
       Add a new track to the end of the track list.
       The track must be allocated via new and is subsequently owned
       by the Event.
       */
      virtual void AddLast(TrackType*);
      
      /**
       Reset event properties to defaults.
       Clears particle list.
       */
      virtual void Reset();
      
      /**
       Compute event properties from the current particle list.
       */
      virtual Bool_t Compute();
      
      // Break off (Hadronic)FinalState(Momentum/Charge) methods into
      // some kind of "interrogator" class?
      
      typedef std::vector<const ::erhic::ParticleMC*> ParticlePtrList;
      
      void FinalState(ParticlePtrList&) const;
      
      /**
       Yields all particles that belong to the hadronic final state.
       This is the same as the result of FinalState(), minus the scattered
       beam lepton.
       */
      void HadronicFinalState(ParticlePtrList&) const;
      
      /**
       Returns the total momentum of the final state.
       */
      TLorentzVector FinalStateMomentum() const;
      
      /**
       Returns the total momentum of the hadronic final state.
       */
      TLorentzVector HadronicFinalStateMomentum() const;
      
      /**
       Returns the total charge of the final state in units of e
       */
      Double_t FinalStateCharge() const;
      
      /**
       Returns a pointer to the incident lepton beam particle.
       Returns a NULL pointer if the particle cannot be located in the event.
       IMPORTANT - DO NOT DELETE THE POINTER OR BAD THINGS WILL HAPPEN!
       
       In the standard eRHIC Monte Carlo format, the incident lepton beam
       is assumed to be the first particle in the particle list.
       Derived classes can implement other selection mechanisms depending on
       their data format.
       */
      virtual const TrackType* BeamLepton() const;
      
      /**
       Returns a pointer to the incident hadron beam particle.
       See also notes in BeamLepton().
       
       In the standard eRHIC Monte Carlo format, the incident hadron beam
       is assumed to be the second particle in the particle list.
       */
      virtual const TrackType* BeamHadron() const;
      
      /**
       Returns a pointer to the exchanged boson.
       See also notes in BeamLepton().
       
       In the standard eRHIC Monte Carlo format, the exchanged boson
       is assumed to be the third particle in the particle list.
       */
      virtual const TrackType* ExchangeBoson() const;
      
      /**
       Returns a pointer to the lepton beam particle after scattering.
       See also notes in BeamLepton().
       
       In the standard eRHIC Monte Carlo format, the scattered lepton beam
       is assumed to be the first final-state particle in the particle list
       with the same PDG code as the incident lepton beam.
       */
      virtual const TrackType* ScatteredLepton() const;
      
      /**
       Populates the event-wise variables from a string.
       Does not populate the particle list or compute
       derived quantities.
       */
      virtual bool Parse(const std::string&);
      
   protected:
      
      Int_t number; ///< Event number
      Int_t process; ///< PYTHIA code for the physics process producing the event
      Int_t nTracks; ///< Number of Particles in the event (intermediate + final)
      
      Double32_t ELeptonInNucl; ///< Incident lepton energy in the nuclear rest frame
      Double32_t ELeptonOutNucl; ///< Scattered lepton energy in the nuclear rest frame
      
      Double32_t x; ///< Bjorken scaling variable
      Double32_t QSquared; ///< Q<sup>2</sup> calculated from scattered electron
      Double32_t y; ///< Inelasticity
      Double32_t WSquared; ///< Invariant mass of the hadronic system
      Double32_t nu; ///< Energy transfer from the electron
      
      Double32_t yJB; ///< y calculated via the Jacquet-Blondel method
      Double32_t QSquaredJB; ///< Q2 calculated via the Jacquet-Blondel method
      Double32_t xJB; ///< x calculated via the Jacquet-Blondel method
      Double32_t WSquaredJB; ///< W2 calculated via the Jacquet-Blondel method
      
      Double32_t yDA; ///< y calculated via the double-angle method
      Double32_t QSquaredDA; ///< Q2 calculated via the double-angle method
      Double32_t xDA; ///< x calculated via the double-angle method
      Double32_t WSquaredDA; ///< W2 calculated via the double-angle method
      
      std::vector< ::erhic::ParticleMC*> particles; ///< Particle list
      
      /**
       Computes kinematics using Jacquet-Blondel method using current
       particle list.
       The argument should contain 4 particles - the incident lepton, hadron,
       the scattered lepton and the virtual photon.
       */
      void ComputeJaquetBlondel(const std::vector<const ::erhic::ParticleMC*>&);
      
      /**
       Computes kinematics using double-angle method using current
       particle list.
       */
      void ComputeDoubleAngle(const std::vector<const ::erhic::ParticleMC*>&);
      
   private:
      
      static int smCount;
      
      ClassDef(EventMC, 1)
   };
   
   inline bool EventMC::Parse(const std::string& ) { return false; }
   
   inline ULong64_t EventMC::N() const { return number; }
   
   inline Int_t EventMC::Process() const { return process; }
   
   inline UInt_t EventMC::NTracks() const { return particles.size(); }
   
   inline const EventMC::TrackType* EventMC::GetTrack(UInt_t u) const {
      return (u < particles.size() ? particles.at(u) : NULL);
   }
   
   inline EventMC::TrackType* EventMC::GetTrack(UInt_t u) {
      return (u < particles.size() ? particles.at(u) : NULL);
   }
   
#if 0
   /**
    Utiltiy methods for analysing the final-state particles of an event.
    */
   template<class T>
   class EventMCFinalState {
      
   public:
      
      typedef T Type;
      typedef std::vector<const T*> ParticlePtrList;
      
      /**
       Returns all the final state particles of this event
       i.e. those with status == 1.
       The returned pointers are those from the event (i.e. not copies)
       so the event must remain in existence while the particles are
       required.
       */
      ParticlePtrList FinalState(const ::erhic::VirtualEvent<Type>&) const;
      
      /**
       Yields all particles that belong to the hadronic final state.
       This is the same as the result of FinalState(), minus the scattered
       beam lepton.
       The returned pointers are those from the event (i.e. not copies)
       so the event must remain in existence while the particles are
       required.
       */
      ParticlePtrList HadronicFinalState(const ::erhic::VirtualEvent<Type>&) const;
      
      /**
       Returns the total momentum of the final state.
       */
      TLorentzVector FinalStateMomentum(const ::erhic::VirtualEvent<Type>&) const;
      
      /**
       Returns the total momentum of the hadronic final state.
       */
      TLorentzVector HadronicFinalStateMomentum(const ::erhic::VirtualEvent<Type>&) const;
      
      /**
       Returns the total charge of the final state in units of e
       */
      Double_t FinalStateCharge(const ::erhic::VirtualEvent<Type>&) const;
      
   };
#endif
   
   
   /**
    Wrapper for getting tree from file and event from tree.
    */
   class Reader {
      
   public:
      
      Reader(const std::string& treeName = "EICTree");
      
      EventMC* Read(Long64_t);
      EventMC* operator()(Long64_t i) { return Read(i); }
      
      TTree* GetTree() { return mTree; }
      
      EventMC* mEvent;
      TTree* mTree;
      
      ClassDef(Reader, 1)
   };
   
//#ifdef USE_NAMESPACE_ERHIC
} // namespace erhic
//#else

#endif
