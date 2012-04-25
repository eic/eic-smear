/**
 VirtualEvent.h
 
 \file
 Declaration of class EventSmear.
 
 \author Michael Savastio
 \date 10/10/11
 \copyright 2011 BNL. All rights reserved.
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
    A generator-independent event with smeared kinematics and particles.
    */
   class Event : public ::erhic::VirtualEvent<Smear::ParticleMCS> {
      
   public:
      
      /** Default constructor */
      Event();
      
      /** Destructor */
      virtual ~Event();
      
      /** Clear the particle list, sets event properties to default values. */
      virtual void Reset();
      
      /** Clears particle array, leaves event variables unchanged */
      virtual void ClearParticles();
      
      /** Returns the number of tracks in the event. */
      virtual UInt_t GetNTracks() const;
      
      /**
       Returns the nth track.
       Returns NULL if the track number is out of the range [0, GetNTracks()).
       @param [in] The track index, in the range [0, GetNTracks()).
       */
      virtual const ParticleMCS* GetTrack(UInt_t) const;
      
      /**
       Returns the nth track.
       Returns NULL if the track number is out of the range [0, GetNTracks()).
       @param [in] The track index, in the range [0, GetNTracks()).
       */
      virtual ParticleMCS* GetTrack(UInt_t);
      
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
      
      /** Returns Bjorken x computed via the double-angle method */
      virtual double GetXDoubleAngle() const;

      /** Returns Q<sup>2</sup> computed via the double-angle method */
      virtual double GetQ2DoubleAngle() const;

      /** Returns inelasticity computed via the double-angle method */
      virtual double GetYDoubleAngle() const;

      /** Returns W<sup>2</sup> computed via the double-angle method */
      virtual double GetW2DoubleAngle() const;
      
      /** Returns Bjorken x computed via the Jacquet-Blondel method */
      virtual double GetXJacquetBlondel() const;

      /** Returns Q<sup>2</sup> computed via the Jacquet-Blondel method */
      virtual double GetQ2JacquetBlondel() const;

      /** Returns inelasticity computed via the Jacquet-Blondel method */
      virtual double GetYJacquetBlondel() const;

      /** Returns W<sup>2</sup> computed via the Jacquet-Blondel method */
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
      
      std::vector<TrackType*> particles; ///< The smeared particle list
      
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
