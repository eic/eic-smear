//
// Event.h
// BuildTree
//
// Created by TB on 10/10/11.
// Copyright 2011 BNL. All rights reserved.
//

#ifndef _ERHIC_EVENT_SMEARED_H_
#define _ERHIC_EVENT_SMEARED_H_

#include <cmath>
#include <list>

#include <TObject.h>

#include "VirtualParticle.h"
#include "ParticleMCSmeared.h"
#include "VirtualEvent.h"
/*
 *  erhic_event.h
 *  BuildTree
 *
 *  Created by TB on 8/19/11.
 *  Copyright 2011 BNL. All rights reserved.
 *
 */

//class Particle;

//#ifdef USE_NAMESPACE_ERHIC
//namespace erhic {
namespace Smear {
   //#endif

   //   class Event : public Event<ParticleMCSmeared> {
   //      
   //   };
   /*
    The smeared events are (currently) generator-independent.
    */
   //   typedef Event<ParticleMCSmeared> Event;
   
   //   class Event : public Event<ParticleMCSmeared> {
   class Event : public ::erhic::VirtualEvent< ::Smear::ParticleMCS> {
      
   public:
      
      void Grow(unsigned);
      
      Event() { }
      
      typedef ::Smear::ParticleMCS SPart;

      std::vector<SPart*> smearedParticles;
      
      virtual bool Parse(const std::string& ) { return false; }
      
      virtual void Reset();/* {
                            for( unsigned i(0); i < smearedParticles.size(); ++i ) {
                            if( GetTrack(i) ) {
                            delete GetTrack(i);
                            GetTrack(i) = NULL;
                            } // if
                            } // for
                            smearedParticles.clear();
                            }*/
      
      //      virtual void Compute();
      
//      typedef TrackType PTyp;
      
//      std::vector<SPart> particles;
      Int_t nTracks; ///< Number of Particles in the event (intermediate + final)
      virtual UInt_t GetNTracks() const { return smearedParticles.size(); }

      void ComputeJaquetBlondel(const std::vector<const ::Smear::ParticleMCS*>&) { }
      void ComputeDoubleAngle(const std::vector<const ::Smear::ParticleMCS*>&) { }
      
      virtual const ::Smear::ParticleMCS* GetTrack(UInt_t u) const {
         return (u < smearedParticles.size() ? smearedParticles.at(u) : NULL);
      }
      virtual ::Smear::ParticleMCS* GetTrack(UInt_t u) {
         return (u < smearedParticles.size() ? smearedParticles.at(u) : NULL);
      }
      
      virtual void AddTrack(TrackType* track) {
         smearedParticles.push_back(track);
      }
      
      void FinalState(std::vector<const ::Smear::ParticleMCS*>&) const {}
      
      
      //   typedef std::vector<const Particle*> ParticleCPtrList;
      
      /**
       Yields all particles that belong to the hadronic final state.
       This is the same as the result of FinalState(), minus the scattered
       beam lepton.
       */
      void HadronicFinalState(std::vector<const ::Smear::ParticleMCS*>&) const { }
      
      
      // Debugging:
      
      /**
       Returns the total momentum of the final state.
       */
      TLorentzVector FinalStateMomentum() const {return TLorentzVector();}
      
      /**
       Returns the total momentum of the hadronic final state.
       */
      TLorentzVector HadronicFinalStateMomentum() const {return TLorentzVector();}
      
      /**
       Returns the total charge of the final state in units of e
       */
      Double_t FinalStateCharge() const {return 0.;}
      /**
       Returns Bjorken-x of the event.
       x<sub>B</sub> = Q<sup>2</sup>/(2p.q)
       */
      virtual Double_t GetX() const { return x; }
      
      /**
       Returns the four-momentum transfer (exchange boson mass) Q<sup>2</sup>.
       Q<sup>2</sup> = 2EE`(1+cos(theta)) = (e-e`)<sup>2</sup>
       */
      virtual Double_t GetQ2() const { return QSquared; }
      
      /**
       Returns the event inelasticity.
       y = (p.q)/(p.e)
       */
      virtual Double_t GetY() const { return y; }
      
      /**
       Returns the invariant mass of the hadronic final state.
       W<sup>2</sup> = M<sup>2</sup> + Q<sup>2</sup>(1-x)/x
       */
      virtual Double_t GetW2() const { return WSquared; }
      
      /**
       Returns the exchange boson energy in the beam hadron rest frame.
       nu = q.p/M
       */
      virtual Double_t GetNu() const { return nu; }
      
      /**
       Add a new track to the end of the track list.
       The track must be allocated via new and is subsequently owned
       by the Event.
       */
      virtual void AddLast(TrackType*);
      
      virtual double GetXDoubleAngle() const { return xDA; }
      virtual double GetQ2DoubleAngle() const { return QSquaredDA; }
      virtual double GetYDoubleAngle() const { return yDA; }
      virtual double GetW2DoubleAngle() const { return WSquaredDA; }
      
      virtual double GetXJacquetBlondel() const { return xJB; }
      virtual double GetQ2JacquetBlondel() const { return QSquaredJB; }
      virtual double GetYJacquetBlondel() const { return yJB; }
      virtual double GetW2JacquetBlondel() const { return WSquaredJB; }
      
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
      
      ClassDef(Event, 1)
   };
   
   //#ifdef USE_NAMESPACE_ERHIC
} // namespace erhic
typedef ::Smear::Event EventS;
//#else
//typedef EventMC EventBase;
//typedef Event EventS;
//#endif


#endif