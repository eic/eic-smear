//
// ParticleMCSmeared.h
// BuildTree
//
// Created by TB on 10/10/11.
// Copyright 2011 BNL. All rights reserved.
//

#ifndef _EICSMEAR_SMEAR_PARTICLEMCS_H_
#define _EICSMEAR_SMEAR_PARTICLEMCS_H_

#include <cmath>

#include <TLorentzVector.h>

#include "eicsmear/erhic/Pid.h"
#include "eicsmear/erhic/VirtualParticle.h"

namespace Smear {

   class Event;
   
   /**
    A smeared Monte Carlo particle.
    */
   class ParticleMCS : public erhic::VirtualParticle {
      
   public:
      
      virtual ~ParticleMCS();
      
      /**
       Default constructor.
       Initialises the Particle from the argument string with the format
       I KS id orig daughter ldaughter px py pz m E xv yv zv
       */
      ParticleMCS();
      
      /** Construct from an E-p 4-vector, pdg code and status code */
      ParticleMCS(const TLorentzVector&, int pdg, int status);

      /**
       Returns the x component of 3-momentum.
       */
      virtual Double_t GetPx() const;
      
      /**
       Returns the y component of 3-momentum.
       */
      virtual Double_t GetPy() const;
      
      /**
       Returns the z component of 3-momentum.
       */
      virtual Double_t GetPz() const;
      
      /**
       Returns the energy of the particle in the lab frame.
       */
      virtual Double_t GetE() const;
      
      /**
       Returns the (E,p) 4-vector in the lab frame.
       */
      virtual TLorentzVector Get4Vector() const;
      
      /**
       Returns the (E,p) 4-vector in the lab frame.
       */
      virtual TLorentzVector PxPyPzE() const { return Get4Vector(); }
      
      /**
       Returns invariant mass (GeV/c<sup>2</sup>).
       */
      virtual Double_t GetM() const;
      
      /**
       Returns momentum transverse to the beam direction.
       \todo check this is set properly
       */
      virtual Double_t GetPt() const;
      
      /**
       Returns the origin point of the particle (cm).
       (0,0,0) indicates a particle originating in the collision.
       */
      virtual TVector3 GetVertex() const;
      
      /**
       Returns the total momentum (GeV).
       */
      virtual Double_t GetP() const;
      
      /**
       Returns the polar angle in the range [0,pi] radians.
       */
      virtual Double_t GetTheta() const;
      
      /**
       Returns the polar angle in the range [0,2pi] radians.
       */
      virtual Double_t GetPhi() const;
      
      /**
       Returns the rapidity.
       */
      virtual Double_t GetRapidity() const;
      
      /**
       Returns the pseudorapidity.
       */
      virtual Double_t GetEta() const;

      /** Returns a status code following the PYTHIA defintion, where
       21 indicates an initial-state particle and 1 indicates a final-
       state particle
      */
      virtual UShort_t GetStatus() const;

      /**
       Returns the ID of the particle.
       */
      virtual ::erhic::Pid Id() const;
      
      virtual void SetE(Double_t);
      
      virtual void SetP(Double_t);
      
      virtual void SetPt(Double_t);
      
      virtual void SetPz(Double_t);
      
      virtual void SetPhi(Double_t);
      
      virtual void SetTheta(Double_t);
      
      virtual UShort_t GetParentIndex() const { return 0; }
//   protected:
      
      UShort_t   status;      ///< Status code
      Int_t      id;				///< PDG particle code
      Double32_t px;          ///< x component of particle momentum
      Double32_t py;          ///< y component of particle momentum
      Double32_t pz;          ///< z component of particle momentum
      Double32_t E;           ///< Energy of particle
      Double32_t pt;          ///< Transverse momentum of particle
      Double32_t p;           ///< Total momentum of particle
      Double32_t theta;       ///< Polar angle
      Double32_t phi;         ///< Azimuthal angle
      
      ClassDef(ParticleMCS, 1)
   };
   
   inline Double_t ParticleMCS::GetPx() const {
      return p * sin(theta) * cos(phi);
   }
   
   inline Double_t ParticleMCS::GetPy() const {
      return p * sin(theta) * sin(phi);
   }
   
   inline Double_t ParticleMCS::GetPz() const { return pz; }
   
   inline Double_t ParticleMCS::GetE() const { return E; }
   
   inline Double_t ParticleMCS::GetM() const {
      return sqrt(pow(E, 2.) - pow(p, 2.));
   }
   
   inline Double_t ParticleMCS::GetPt() const { return pt; }
   
   inline TVector3 ParticleMCS::GetVertex() const { return TVector3(); }
   
   inline Double_t ParticleMCS::GetP() const { return p; }
   
   inline Double_t ParticleMCS::GetTheta() const { return theta; }
   
   inline Double_t ParticleMCS::GetPhi() const { return phi; }
   
   inline Double_t ParticleMCS::GetRapidity() const { return 0.; }
   
   inline Double_t ParticleMCS::GetEta() const { return 0.; }
 
   inline UShort_t ParticleMCS::GetStatus() const { return status; }

   inline void ParticleMCS::SetE(Double_t e) { E = e; }
   
   inline void ParticleMCS::SetP(Double_t momentum) { p = momentum; }
   
   inline void ParticleMCS::SetPt(Double_t momentum) { pt = momentum; }
   
   inline void ParticleMCS::SetPz(Double_t momentum) { pz = momentum; }
   
   inline void ParticleMCS::SetPhi(Double_t value) { phi = value; }
   
   inline void ParticleMCS::SetTheta(Double_t value) { theta = value; }
   
   inline erhic::Pid ParticleMCS::Id() const { return ::erhic::Pid(id); }
   
} // namespace erhic

#endif // _EICSMEAR_SMEAR_PARTICLEMCS_H_
