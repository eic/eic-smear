/**
 VirtualEvent.h
 
 \file
 Declaration of class EventSmear.
 
 \author Michael Savastio
 \date 10/10/11
 \copyright 2011 BNL. All rights reserved.
 
 \todo Tidy up
 */

#ifndef _ERHIC_EVENT_SMEARED_H_
#define _ERHIC_EVENT_SMEARED_H_

#include <cmath>
#include <list>

#include <TObject.h>

#include "ParticleMCS.h"
#include "VirtualEvent.h"
#include "VirtualParticle.h"

namespace Smear {
   
   /*
    The smeared events are (currently) generator-independent.
    \todo Add comments
    */
   class Event : public ::erhic::VirtualEvent<Smear::ParticleMCS> {
      
   public:
      
      Event() { }
      
      virtual void Reset();
      
      virtual UInt_t GetNTracks() const;
      
      virtual const Smear::ParticleMCS* GetTrack(UInt_t u) const;
      
      virtual Smear::ParticleMCS* GetTrack(UInt_t u);
      
      virtual void AddTrack(TrackType* track);
      
      /**
       Returns Bjorken-x of the event.
       x<sub>B</sub> = Q<sup>2</sup>/(2p.q)
       */
      virtual Double_t GetX() const;
      
      /**
       Returns the four-momentum transfer (exchange boson mass) Q<sup>2</sup>.
       Q<sup>2</sup> = 2EE`(1+cos(theta)) = (e-e`)<sup>2</sup>
       */
      virtual Double_t GetQ2() const;
      
      /**
       Returns the event inelasticity.
       y = (p.q)/(p.e)
       */
      virtual Double_t GetY() const;
      
      /**
       Returns the invariant mass of the hadronic final state.
       W<sup>2</sup> = M<sup>2</sup> + Q<sup>2</sup>(1-x)/x
       */
      virtual Double_t GetW2() const;
      
      /**
       Returns the exchange boson energy in the beam hadron rest frame.
       nu = q.p/M
       */
      virtual Double_t GetNu() const;
      
      /**
       Add a new track to the end of the track list.
       The track must be allocated via new and is subsequently owned
       by the Event.
       */
      virtual void AddLast(TrackType*);
      
      virtual double GetXDoubleAngle() const;
      virtual double GetQ2DoubleAngle() const;
      virtual double GetYDoubleAngle() const;
      virtual double GetW2DoubleAngle() const;
      
      virtual double GetXJacquetBlondel() const;
      virtual double GetQ2JacquetBlondel() const;
      virtual double GetYJacquetBlondel() const;
      virtual double GetW2JacquetBlondel() const;
      
//   protected:
      
      Int_t nTracks; ///< Number of particles (intermediate + final)
      
      Double32_t x;           ///< Bjorken scaling variable
      Double32_t QSquared;    ///< Q<sup>2</sup> calculated from
                              ///< scattered electron
      Double32_t y;           ///< Inelasticity
      Double32_t WSquared;    ///< Invariant mass of the hadronic system
      Double32_t nu;          ///< Energy transfer from the electron
      
      Double32_t yJB;         ///< y calculated via the Jacquet-Blondel method
      Double32_t QSquaredJB;  ///< Q2 calculated via the Jacquet-Blondel method
      Double32_t xJB;         ///< x calculated via the Jacquet-Blondel method
      Double32_t WSquaredJB;  ///< W2 calculated via the Jacquet-Blondel method
      
      Double32_t yDA;         ///< y calculated via the double-angle method
      Double32_t QSquaredDA;  ///< Q2 calculated via the double-angle method
      Double32_t xDA;         ///< x calculated via the double-angle method
      Double32_t WSquaredDA;  ///< W2 calculated via the double-angle method
      
      std::vector<TrackType*> particles;
      
      ClassDef(Event, 1)
   };
   
   inline UInt_t Event::GetNTracks() const {
      return particles.size();
   }
   
   inline const Smear::ParticleMCS* Event::GetTrack(UInt_t u) const {
      return (u < particles.size() ? particles.at(u) : NULL);
   }
   
   inline Smear::ParticleMCS* Event::GetTrack(UInt_t u) {
      return (u < particles.size() ? particles.at(u) : NULL);
   }
   
   inline void Event::AddTrack(TrackType* track) {
      particles.push_back(track);
   }
   
   inline Double_t Event::GetX() const {
      return x;
   }
   
   inline Double_t Event::GetQ2() const {
      return QSquared;
   }

   inline Double_t Event::GetY() const {
      return y;
   }
   
   inline Double_t Event::GetW2() const {
      return WSquared;
   }
   
   inline Double_t Event::GetNu() const {
      return nu;
   }
   
   inline double Event::GetXDoubleAngle() const {
      return xDA;
   }
   
   inline double Event::GetQ2DoubleAngle() const {
      return QSquaredDA;
   }
   
   inline double Event::GetYDoubleAngle() const {
      return yDA;
   }
   
   inline double Event::GetW2DoubleAngle() const {
      return WSquaredDA;
   }
   
   inline double Event::GetXJacquetBlondel() const {
      return xJB;
   }
   
   inline double Event::GetQ2JacquetBlondel() const {
      return QSquaredJB;
   }
   
   inline double Event::GetYJacquetBlondel() const {
      return yJB;
   }
   
   inline double Event::GetW2JacquetBlondel() const {
      return WSquaredJB;
   }
   
} // namespace Smear

typedef Smear::Event EventS;

#endif
