#ifndef _EICSMEAR_KINEMATICS_
#define _EICSMEAR_KINEMATICS_

//
// Kinematics.h
//
// Created by TB on 7/7/11.
// Copyright 2011 BNL. All rights reserved.
//

#include <list>

#include <Rtypes.h>
#include <TLorentzVector.h>

#include "eicsmear/erhic/VirtualParticle.h"
#include "eicsmear/erhic/EventBase.h"
//#include "eicsmear/erhic/EventMC.h"
#include "eicsmear/erhic/ParticleMC.h"
//class EventBase;
#include "eicsmear/erhic/ParticleIdentifier.h"
#include "eicsmear/erhic/BeamParticles.h"

//class EventBase;

namespace erhic {

   class EventDis;

   // Needed? Unless there is some functionality or interface associated
   // with Kinematics, just inherit subclasses directly from TObject
   struct Kinematics : public TObject {
      ClassDef(erhic::Kinematics, 1)
   };

   struct DisKinematics : public Kinematics {
      DisKinematics();
      DisKinematics(double x, double y, double nu, double Q2, double W2);
      Double32_t mX;
      Double32_t mQ2;
      Double32_t mW2;
      Double32_t mNu;
      Double32_t mY;
      ClassDef(erhic::DisKinematics, 1)
   };

   /**
    Abstract base class for computations of event kinematics.
   */
   class KinematicsComputer {
   public:
      virtual ~KinematicsComputer() { }
      virtual Kinematics* Calculate() = 0;
      ClassDef(KinematicsComputer, 1)
   };

   /**
    Computes DIS event kinematics from the scattered lepton.
   */
   class LeptonKinematicsComputer : public KinematicsComputer {
   public:
      virtual ~LeptonKinematicsComputer() { }
      /** Initialise with the beam info used to compute the kinematics */
      LeptonKinematicsComputer(const BeamParticles&);
      /** Determine the beam info from the input event */
      LeptonKinematicsComputer(const EventDis&);
      virtual DisKinematics* Calculate();
   protected:
      BeamParticles mBeams;
      ClassDef(LeptonKinematicsComputer, 1)
   };

   /**
    Computes DIS event kinematics from final-state hadrons using
    the Jacquet-Blondel method.
   */
   class JacquetBlondelComputer : public KinematicsComputer {
   public:
      virtual ~JacquetBlondelComputer() { }
      /**
       Initialise with the event to compute.
       If the second argument is non-NULL, use the beam information from it
       in the computation.
       If it is NULL, determine the beam information automatically from
       the event.
       This allows the same class to be used with smeared calculations, where
       the beam information isn't associated with the smeared event itself.
      */
      JacquetBlondelComputer(const EventDis& e, const BeamParticles*);
      virtual DisKinematics* Calculate();
   protected:
      const EventDis& mEvent;
      BeamParticles mBeams;
      ClassDef(JacquetBlondelComputer, 1)
   };

   /**
    Computes DIS event kinematics from a mixture of hadronic and lepton
    variables using the double-angle method.
   */
   class DoubleAngleComputer : public KinematicsComputer {
   public:
      virtual ~DoubleAngleComputer() { }
      /**
       Initialise with the event to compute.
       If the second argument is non-NULL, use the beam information from it
       in the computation.
       If it is NULL, determine the beam information automatically from
       the event.
       This allows the same class to be used with smeared calculations, where
       the beam information isn't associated with the smeared event itself.
      */
      DoubleAngleComputer(const EventDis& e, const BeamParticles*);
      virtual DisKinematics* Calculate();
   protected:
      const EventDis& mEvent;
      BeamParticles mBeams;
      ClassDef(DoubleAngleComputer, 1)
   };
} // namespace erhic

//	=================================================================================================
// Base class for calculation of event kinematics from the hadronic final state.
// Derived classes JacquetBlondel and DoubleAngle are declared below.
//	=================================================================================================
class KinematicsFromHadrons {
   
public:
   
   // This is just for testing/debug while I change the code.
   //   static void BeExact(bool b = true) { smBeExact = b; }
   
   KinematicsFromHadrons();
   virtual ~KinematicsFromHadrons() { }
   
   /**
    Deprecated, use setBeamLepton()
    */
//   virtual void setLeptonEnergy(const Double_t& );
   
   /**
    Deprecated, use setBeamLepton() and setBeamHadron()
    */
//   virtual void setMandelstamS(const Double_t& );
   virtual void addParticle(const TLorentzVector& );
   virtual void clearParticles();
   
   virtual Double_t getLeptonEnergy() const; // Electron beam energy
   virtual Double_t getMandelstamS() const; // Mandelstam s, squared centre-of-mass energy
   
   virtual Double_t computeY() const = 0;
   virtual Double_t computeQSquared() const = 0;
   virtual Double_t computeX() const = 0;
   
   virtual void setBeamLepton(const TLorentzVector& );
   virtual void setBeamHadron(const TLorentzVector& );
   
protected:
   
//   Double_t mLeptonEnergy;
//   Double_t mMandelstamS;
   TLorentzVector mBeamLepton;
   TLorentzVector mBeamHadron;
   std::list<TLorentzVector> mParticles;
};

inline void KinematicsFromHadrons::setBeamLepton(const TLorentzVector& v) {
   mBeamLepton = v;
}

inline void KinematicsFromHadrons::setBeamHadron(const TLorentzVector& v) {
   mBeamHadron = v;
}

//	=================================================================================================
// Class for calculating kinematic quantities using the Jacquet-Blondel method.
// Set the lepton beam energy (required for y) and the squared centre of mass energy (for Q2).
// Use addParticle() to store particles from the hadronic final state.
// y(), QSquared() and x() calculate the named variables using the current particle list.
// TODO tidy up
//	=================================================================================================
class JacquetBlondel : public KinematicsFromHadrons {
   
public:
   
   JacquetBlondel();
   virtual ~JacquetBlondel() { }
   
   virtual Double_t computeY() const;
   virtual Double_t computeQSquared() const;
   virtual Double_t computeX() const;
   
protected:
   
   virtual Double_t computeYExact() const;
   virtual Double_t computeQSquaredExact() const;
};

//	=================================================================================================
// Class for calculating kinematic quantities using the double-angle method.
// Set the lepton beam energy and Mandelstam s as for Jacquet-Blondel, plus the lepton
// scattering angle and the proton beam energy.
//	=================================================================================================
class DoubleAngle : public KinematicsFromHadrons {
   
public:
   
   DoubleAngle();
   virtual ~DoubleAngle() { }
   
   virtual void setLeptonAngle(const Double_t& a ) { mElectronAngle = a; }
   
   /**
    Deprecated, use setBeamHadron()
    */
//   virtual void setProtonEnergy(const Double_t& );//e ) { mProtonEnergy = e; }
   
   virtual void addParticle(const TLorentzVector& p );
   
   virtual Double_t getLeptonAngle() const { return mElectronAngle; }
   
   virtual Double_t getProtonEnergy() const;// { return mProtonEnergy; }
   
   virtual Double_t computeQuarkAngle() const; // Scattering angle of struck quark
   
   virtual Double_t computeY() const;
   virtual Double_t computeQSquared() const;
   virtual Double_t computeX() const;
   
protected:
   
   mutable Bool_t mHasChanged;
   mutable Double_t mAngle;
   Double_t mElectronAngle;
//   Double_t mProtonEnergy;
};

inline Double_t DoubleAngle::getProtonEnergy() const {
   return mBeamHadron.E();
}

class Kinematics {
   
public:
   
   Double32_t x;
   Double32_t QSquared;
   Double32_t y;
   Double32_t nu;
   Double32_t WSquared;
};

#if 0
namespace erhic {
   
   class Kinematics {
      
   public:
      
      double y;
      double Q2;
      double x;
      double W2;
      double nu;
   };
   
   /**
    Pure virtual base class defining methods for calculating event kinematics.
    */
   class KinematicsCalculation {
      
   public:
      
      virtual ~KinematicsCalculation() { }
#if 0
      virtual Double_t ComputeX(const EventBase&) const = 0;
      virtual Double_t ComputeQSquared(const EventBase&) const = 0;
      virtual Double_t ComputeY(const EventBase&) const = 0;
      virtual Double_t ComputeWSquared(const EventBase&) const = 0;
      virtual Double_t ComputeNu(const EventBase&) const = 0;
#endif
      typedef Double_t (::Particle::*Getter)() const;
      /**
       Compute Bjorken x from an event.
       The second argument allows selection of which variable to use as
       "energy" in calculations, as momentum is often the more sensible option
       in cases where detector resolution is superior for momentum than energy
       */
      virtual Double_t ComputeX(const EventBase&, Getter = &::Particle::GetE) const = 0;
      virtual Double_t ComputeQSquared(const EventBase&, Getter = &::Particle::GetE) const = 0;
      virtual Double_t ComputeY(const EventBase&, Getter = &::Particle::GetE) const = 0;
      virtual Double_t ComputeWSquared(const EventBase&, Getter = &::Particle::GetE) const = 0;
      virtual Double_t ComputeNu(const EventBase&, Getter = &::Particle::GetE) const = 0;
      
      /*
       Switch between using energy or momentum of particles for their
       "energy" in calculations.
       It makes sense to use momentum in cases where detector resolution for
       momentum is better than for energy.
       To really mimic the detector effects, when using the momentum, the
       energy should be calculated from the measured (smeared) momentum and
       the mass, based on the "smeared" id (i.e. potentially misidentified).
       Options for when there is no PID - use energy even if p requested;
       assume it is a pion.
       */
#if 0
      virtual void UseMomentum();
      virtual void UseEnergy();
      /*OR (more general input for E!*/
      virtual void UseForEnergy(double (Particle::*)() const);
      /*OR limited compared to above but perhaps simpler for novice users*/
      enum { kE, kP };
      virtual void UseForEnergy(int = kE);
      /*OR (keeps base class pure virtual*/
      typedef Double_t (Particle::*Getter)() const;
      /**
       Compute Bjorken x from an event.
       The second argument allows selection of which variable to use as
       "energy" in calculations, as momentum is often the more sensible option
       in cases where detector resolution is superior for momentum than energy
       */
      virtual Double_t ComputeX(const EventBase&, Getter = &Particle::GetE) const = 0;
      virtual Double_t ComputeQSquared(const EventBase&, Getter = &Particle::GetE) const = 0;
      virtual Double_t ComputeY(const EventBase&, Getter = &Particle::GetE) const = 0;
      virtual Double_t ComputeWSquared(const EventBase&, Getter = &Particle::GetE) const = 0;
      virtual Double_t ComputeNu(const EventBase&, Getter = &Particle::GetE) const = 0;
      virtual void Compute(const EventBase&, Getter = &Particle::GetE) const = 0;
      virtual Kinematics Compute(const EventBase&, Getter = &Particle::GetE) const = 0;
      /*OR*/
      virtual Double_t ComputeX(const EventBase&, int = kE) const = 0;
#endif
      
   protected:
      
   private:
      
   };
   
   /**
    Calculations of event kinematic quantities using the scattered lepton.
    */
   class LeptonKinematics : public KinematicsCalculation {
      
   public:
      
      LeptonKinematics();
      
      virtual ~LeptonKinematics();
#if 0
      virtual Double_t ComputeX(const EventBase&) const;
      virtual Double_t ComputeQSquared(const EventBase&) const;
      virtual Double_t ComputeY(const EventBase&) const;
      virtual Double_t ComputeWSquared(const EventBase&) const;
      virtual Double_t ComputeNu(const EventBase&) const;
#endif 
      /**
       Compute Bjorken x from an event.
       The second argument allows selection of which variable to use as
       "energy" in calculations, as momentum is often the more sensible option
       in cases where detector resolution is superior for momentum than energy
       */
      virtual Double_t ComputeX(const EventBase&, Getter = &::Particle::GetE) const = 0;
      virtual Double_t ComputeQSquared(const EventBase&, Getter = &::Particle::GetE) const = 0;
      virtual Double_t ComputeY(const EventBase&, Getter = &::Particle::GetE) const = 0;
      virtual Double_t ComputeWSquared(const EventBase&, Getter = &::Particle::GetE) const = 0;
      virtual Double_t ComputeNu(const EventBase&, Getter = &::Particle::GetE) const = 0;
   protected:
      
   private:
      
   };
   
   
   /**
    Calculations of event kinematics using the Jacquet-Blondel method.
    Requires the first particle in the event to be the beam lepton and the
    second particle to be the beam hadron.
    */
   class JacquetBlondel : public KinematicsCalculation {
      
   public:
      
      JacquetBlondel();
      
      virtual ~JacquetBlondel();
#if 0
      virtual Double_t ComputeX(const EventBase&) const;
      virtual Double_t ComputeQSquared(const EventBase&) const;
      virtual Double_t ComputeY(const EventBase&) const;
      virtual Double_t ComputeWSquared(const EventBase&) const;
      virtual Double_t ComputeNu(const EventBase&) const;
#endif
      virtual Double_t ComputeX(const EventBase&, Getter = &::Particle::GetE) const;
      virtual Double_t ComputeQSquared(const EventBase&, Getter = &::Particle::GetE) const;
      virtual Double_t ComputeY(const EventBase&, Getter = &::Particle::GetE) const;
      virtual Double_t ComputeWSquared(const EventBase&, Getter = &::Particle::GetE) const;
      virtual Double_t ComputeNu(const EventBase&, Getter = &::Particle::GetE) const;
      
   protected:
      
   private:
      
   };
   
   
   /**
    Calculations of event kinematics using the Jacquet-Blondel method.
    */
   class DoubleAngle : public KinematicsCalculation {
      
   public:
      
      DoubleAngle();
      
      virtual ~DoubleAngle();
#if 0
      virtual Double_t ComputeX(const EventBase&) const;
      virtual Double_t ComputeQSquared(const EventBase&) const;
      virtual Double_t ComputeY(const EventBase&) const;
      virtual Double_t ComputeWSquared(const EventBase&) const;
      virtual Double_t ComputeNu(const EventBase&) const;
#endif
      virtual Double_t ComputeX(const EventBase&, Getter = &::Particle::GetE) const = 0;
      virtual Double_t ComputeQSquared(const EventBase&, Getter = &::Particle::GetE) const = 0;
      virtual Double_t ComputeY(const EventBase&, Getter = &::Particle::GetE) const = 0;
      virtual Double_t ComputeWSquared(const EventBase&, Getter = &::Particle::GetE) const = 0;
      virtual Double_t ComputeNu(const EventBase&, Getter = &::Particle::GetE) const = 0;
      
   protected:
      
   private:
      
   };
   
} // namespace erhic

#endif

#endif
