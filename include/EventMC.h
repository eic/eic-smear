/**
 EventMC.h
 
 \file
 Declaration of class EventMC.
 
 \author TB
 \date 10/10/11
 \copyright 2011 BNL. All rights reserved.
 */

#ifndef _ERHIC_EVENTMC_H_
#define _ERHIC_EVENTMC_H_

#include <TLorentzVector.h>

#include "DisEvent.h"
#include "ParticleMC.h"

class TTree;

namespace erhic {
   
   /**
    Abstract base class for DIS Monte Carlo events.
    Implements common event properties and methods.
    */
   class EventMC : public DisEvent {
      
   public:
      
      /** Constructor. */
      EventMC();
      
      /** Destructor. */
      virtual ~EventMC();
      
      /** Returns a unique identifier for this event. */
      virtual ULong64_t GetN() const;
      
      /**
       Returns a code describing the production process of this event.
       */
      virtual Int_t GetProcess() const;
      
      /**
       Returns the number of tracks in the event.
       */
      virtual UInt_t GetNTracks() const;
      
      /**
       Returns the nth track.
       Returns NULL if the track number is out of the range [0, GetNTracks()).
       @param [in] The track index, in the range [0, GetNTracks()).
       */
      virtual const ParticleMC* GetTrack(UInt_t) const;
      
      /**
       \overload const ParticleMC* GetTrack(UInt_t) const
       */
      virtual ParticleMC* GetTrack(UInt_t);
      
      /**
       Returns a pointer to the incident lepton beam particle.
       Returns a NULL pointer if the particle cannot be located in the event.
       IMPORTANT - DO NOT DELETE THE POINTER OR BAD THINGS WILL HAPPEN!
       
       In the standard eRHIC Monte Carlo format, the incident lepton beam
       is assumed to be the first particle in the particle list.
       This is the behaviour implemented here.
       Derived classes can implement other selection mechanisms depending on
       their data format.
       */
      virtual const ParticleMC* BeamLepton() const;
      
      /**
       Returns a pointer to the incident hadron beam particle.
       See also notes in BeamLepton().
       
       In the standard eRHIC Monte Carlo format, the incident hadron beam
       is assumed to be the second particle in the particle list.
       */
      virtual const ParticleMC* BeamHadron() const;
      
      /**
       Returns a pointer to the exchanged boson.
       See also notes in BeamLepton().
       
       In the standard eRHIC Monte Carlo format, the exchanged boson
       is assumed to be the third particle in the particle list.
       */
      virtual const ParticleMC* ExchangeBoson() const;
      
      /**
       Returns a pointer to the lepton beam particle after scattering.
       See also notes in BeamLepton().
       
       In the standard eRHIC Monte Carlo format, the scattered lepton beam
       is assumed to be the first final-state particle in the particle list
       with the same PDG code as the incident lepton beam.
       */
      virtual const ParticleMC* ScatteredLepton() const;
      
      /**
       Populates the event-wise variables from a string.
       Does not populate the particle list or compute derived quantities.
       See also Compute().
       */
      virtual bool Parse(const std::string&) = 0;
      
      /**
       Add a new track to the end of the track list.
       The track must be allocated via new.
       It is subsequently owned and deleted by the Event.
       @param [in] Pointer to the track to add.
       */
      virtual void AddLast(ParticleMC* track);
      
      /**
       Resets event properties to defaults.
       Clears particle list (freeing any allocated memory).
       */
      virtual void Reset();
      
      /**
       Sets the code describing the production process of this event.
       @param [in] code The identifying code, an integer
       */
      virtual void SetProcess(int code);
      
      /**
       Sets the unique identifier for this event.
       @param [in] n The identifying number, an integer
       */
      virtual void SetN(int n);
      
      /**
       Sets the track count for this event.
       @param [in] n The track count, an integer
       */
      virtual void SetNTracks(int n);
      
      /** */
      virtual void SetELeptonInNuclearFrame(double e) { ELeptonInNucl = e; }
      
      /** */
      virtual void SetEScatteredInNuclearFrame(double e) { ELeptonOutNucl = e; }
      
      /**
       Stores pointers to all final state particles in the list.
       These pointers should not be deleted by the user.
       Any existing entries in the list are not changed.
       @param [out] particles The list in which to store particles.
       */
      void FinalState(ParticlePtrList& particles) const;
      
      /**
       Yields all particles that belong to the hadronic final state.
       This is the same as the result of FinalState(), minus the scattered
       beam lepton.
       */
      void HadronicFinalState(ParticlePtrList&) const;
      
      /** Returns the total momentum of the final state in GeV/c. */
      TLorentzVector FinalStateMomentum() const;
      
      /** Returns the total momentum of the hadronic final state in GeV/c. */
      TLorentzVector HadronicFinalStateMomentum() const;
      
      /** Returns the total charge of the final state in units of e */
      Double_t FinalStateCharge() const;
      
      //   protected:
      Int_t number; ///< Event number
      Int_t process; ///< PYTHIA code for the physics process producing the event
      Int_t nTracks; ///< Number of Particles in the event (intermediate + final)
      
      Double32_t ELeptonInNucl; ///< Incident lepton energy in the nuclear rest frame
      Double32_t ELeptonOutNucl; ///< Scattered lepton energy in the nuclear rest frame
      
      std::vector<ParticleMC*> particles; ///< Particle list
      
   private:
      
      ClassDef(EventMC, 1)
   };
   
   inline ULong64_t EventMC::GetN() const {
      return number;
   }
   
   inline Int_t EventMC::GetProcess() const {
      return process;
   }
   
   inline UInt_t EventMC::GetNTracks() const {
      return particles.size();
   }
   
   inline const ParticleMC* EventMC::GetTrack(UInt_t u) const {
      return (u < particles.size() ? particles.at(u) : NULL);
   }
   
   inline ParticleMC* EventMC::GetTrack(UInt_t u) {
      return (u < particles.size() ? particles.at(u) : NULL);
   }
   
   inline void EventMC::SetProcess(int code) {
      process = code;
   }
   
   inline void EventMC::SetN(int n) {
      number = n;
   }
   
   inline void EventMC::SetNTracks(int n) {
      nTracks = n;
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
      
      virtual ~Reader() { }
      
      EventMC* Read(Long64_t);
      EventMC* operator()(Long64_t i) { return Read(i); }
      
      TTree* GetTree() { return mTree; }
      
      EventMC* mEvent;
      TTree* mTree;
      
      ClassDef(Reader, 1)
   };
   
} // namespace erhic

#endif
