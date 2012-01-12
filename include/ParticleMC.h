//
// ParticleMC.h
// BuildTree
//
// Created by TB on 10/10/11.
// Copyright 2011 BNL. All rights reserved.
//

#ifndef _ERHIC_ParticleMC_H_
#define _ERHIC_ParticleMC_H_

//#include <functional>
//#include <iostream>
//#include <memory>
#include <string>

#include <TLorentzVector.h>
#include <TRef.h>
//#include <TRefArray.h>
#include <TVector3.h>

#include "Pid.h"
#include "VirtualParticle.h"

//class TRefArray;

namespace erhic {
   
   
   class EventMC;
   
   
   class ParticleMC : public VirtualParticle {
   public:
      
      virtual ~ParticleMC();
      
      static Long64_t smNInstances;
      
      /**
       Default constructor.
       Initialises the Particle from the argument string with the format
         I KS id orig daughter ldaughter px py pz m E xv yv zv
       */
      ParticleMC(const std::string& = "");
      
      /**
       Print the contents of Particle to standard output.
       The format is that of the input Monte Carlo i.e.
       I KS id orig daughter ldaughter px py pz E m xv yv zv.
       Inherited from TObject. The argument is unused.
       */
      virtual void Print(Option_t* = "") const;
      
      /**
       Returns the particle index in an event.
       */
      virtual UInt_t Index() const;
      
      /**
       Returns the status of the particle.
       See the description of variable K(I,1) in the PYTHIA manual for the
       meanings of codes.
       */
      virtual UShort_t Status() const;
      
      /**
       Returns the index of this particles parent in an event.
       */
      virtual UShort_t ParentIndex() const;
      
      /**
       Returns a pointer to the parent of this particle.
       This is the particle with index ParentIndex() returned from the
       event containing this particle (obtainable via GetEvent()).
       Returns NULL if this particle has no parent or it cannot
       be accessed via GetEvent().
       */
      virtual const ParticleMC* GetParent() const;
      
      /**
       Returns the index of this particle's first child particle.
       Returns 0 if this particle has no children.
       */
      virtual UShort_t Child1() const;
      
      /**
       Returns the index of this particle's last child particle.
       Returns 0 if this particle has zero or one children.
       */
      virtual UShort_t ChildN() const;
      
      /**
       Returns the number of children of this particle.
       Returns 0 if the particle did not decay.
       */
      virtual UInt_t NChildren() const;
      
      /**
       Returns a pointer to the nth child particle of this particle,
       where n is in the range [0, NChildren()).
       GetChild(0) returns the child whose index in this particle's
       event is Child1().
       GetChild(NChildren()-1) returns the child whose index
       is ChildN().
       Returns NULL if there is no such child particle or it
       cannot be accessed via the event (see GetEvent()).
       */
      virtual const ParticleMC* GetChild(UShort_t) const;
      
      /**
       Returns true if n is <= the number of children of this particle.
       Returns false otherwise.
       Equivalent to NULL != GetChild(n).
       */
      virtual Bool_t HasChild(Int_t) const;
      
      /**
       Returns the x component of 3-momentum.
       */
      virtual Double_t Px() const;
      
      /**
       Returns the y component of 3-momentum.
       */
      virtual Double_t Py() const;
      
      /**
       Returns the z component of 3-momentum.
       */
      virtual Double_t Pz() const;
      
      /**
       Returns invariant mass (GeV/c<sup>2</sup>).
       */
      virtual Double_t M() const;
      
      /**
       Returns momentum transverse to the beam direction.
       */
      virtual Double_t Pt() const;
      
      /**
       Returns the origin point of the particle (cm).
       (0,0,0) indicates a particle originating in the collision.
       */
      virtual TVector3 Vertex() const;
      
      /**
       Returns the PDG code of this particle's parent.
       
       */
      virtual erhic::Pid ParentId() const;
      
      /**
       Returns the total momentum (GeV).
       */
      virtual Double_t P() const;
      
      /**
       Returns the polar angle in the range [0,pi] radians.
       */
      virtual Double_t Theta() const;
      
      /**
       Returns the polar angle in the range [0,2pi] radians.
       */
      virtual Double_t Phi() const;
      
      /**
       Returns the rapidity.
       */
      virtual Double_t Y() const;
      
      /**
       Returns the pseudorapidity.
       */
      virtual Double_t Eta() const;
      
      /**
       Returns the variable z.
       z = (P.p_h)/(P.q).
       */
      virtual Double_t Z() const;
      
      /**
       Returns Feynman-x.
       x<sub>F</sub> = 2*p<sub>z</sub>/sqrt(s).
       */
      virtual Double_t Xf() const;
      
      /**
       Returns the angle with respect to the exchange boson.
       Defined in the beam hadron's rest frame.
       Given in the range [0,pi] radians.
       */
      virtual Double_t AngleVsBoson() const;
      
      /**
       Returns the p<sub>T</sub> with respect to the exchange boson.
       Defined in the beam hadron's rest frame.
       */
      virtual Double_t PtVsBoson() const;
      
      /**
       Returns a pointer to the event containing this particle.
       Returns NULL if this particle has not been associated
       with an event.
       */
      const EventMC* GetEvent() const;
      
      void SetEvent(EventMC*);
      
      /**
       Returns the (E,p) 4-vector in the lab frame.
       */
      virtual TLorentzVector PxPyPzE() const;
      virtual TLorentzVector Get4Vector() const { return PxPyPzE(); }
      /**
       Returns the (E,p) 4-vector in the hadron-boson frame.
       This frame is defined such that
       <UL>
       <LI> the beam hadron is at rest </LI>
       <LI> the z direction is the exchange boson momentum vector </LI>
       <LI> the y direction is defined as q x e`, where q is the boson and
       e` is the scattered lepton momentum </LI>
       <LI> x is defined to complete the right-handed coordinate system </LI>
       </UL>
       */
      virtual TLorentzVector Get4VectorInHadronBosonFrame() const;
      
      /**
       Returns the invariant mass of the particle.
       */
      virtual Double_t GetM() const;
      
      /**
       Returns the energy of the particle in the lab frame.
       */
      virtual Double_t GetE() const;
      
//      virtual void SetE(Double_t);
      
//      virtual void SetP(Double_t);
      
//      virtual void SetPt(Double_t);
      
//      virtual void SetPz(Double_t);
      
//      virtual void SetPhi(Double_t);
      
//      virtual void SetTheta(Double_t);
      
//      virtual void SetId(Int_t);
      
//      virtual void SetStatus(UShort_t);
      
      /**
       Returns the ID of the particle.
       */
      virtual Pid Id() const;
      
      /**
       Sets quantities derived from (E, px, py, pz).
       Namely:
       total momentum, transverse momentum, rapidity, pseudorapidity, theta,
       phi. This should be called if (E, px, py, pz) are manually altered.
       */
      virtual void ComputeDerivedQuantities();
      
      /**
       Sets quantities that depend on the contents of the whole event.
       Namely:
       z, Feynman x, angle and pt with respect to the exchanged boson, 
       azimuthal angle around the exchanged boson, parent pdg code.
       Important: this particle is assumed to be in the same frame of
       reference as those in the event argument.
       */
      virtual void ComputeEventDependentQuantities(EventMC&);
      
//   protected:
      
      UShort_t    I;             ///< Particle index in event
      UShort_t    KS;				///< Particle status code: see PYTHIA manual
      Int_t       id;				///< PDG particle code
      UShort_t    orig;          ///< I of parent particle
      UShort_t    daughter;      ///< I of first child particle
      UShort_t    ldaughter;     ///< I of last child particle
      
      Double32_t px;           ///< x component of particle momentum
      Double32_t py;           ///< y component of particle momentum
      Double32_t pz;           ///< z component of particle momentum
      
      // Reducing file size.
      // The following variables could be removed, at the cost of increased
      // running time to calculate them on-the-fly:
      // p - compute from px/y/z
      // m - accessible via PID with TDatabasePDG
      // E - compute from p and m
      // parentId - accessible via event record and parent index IF a
      //          persistent pointer to the containing event can be stored
      // daughter - IF daughter particles always start immediately after this one
      // pt - from px, py
      // theta - from pt, pz
      // phi - from px, py
      // rapidity - from E, pz
      // eta - from theta
      // xFeynman - if event is accessible
      // thetaGamma - duplicated in pHadronBoson
      // ptVsGamma - duplicated in pHadronBoson
      // phiPrf - duplicated in pHadronBoson
      
      Double32_t E;            ///< Energy of particle
      Double32_t m;            ///< Invariant mass of particle
      Double32_t pt;           ///< Transverse momentum of particle
      Double32_t xv;           ///< x coordinate of particle production vertex
      Double32_t yv;           ///< y coordinate of particle production vertex
      Double32_t zv;           ///< z coordinate of particle production vertex
      
      // Can be deprecated if parent particle itself is available
      Int_t parentId;          ///< PDG code of this particle's parent
      
      Double32_t p;            ///< Total momentum of particle
      Double32_t theta;        ///< Polar angle
      Double32_t phi;          ///< Azimuthal angle
      Double32_t rapidity;     ///< Rapidity of particle
      Double32_t eta;          ///< Pseudorapidity of particle
      Double32_t z;            ///< Fraction of virtual photon energy
                               ///< carried by particle
      Double32_t xFeynman;     ///< Feynman x = p<sub>z</sub>/(2sqrt(s))
      Double32_t thetaGamma;   ///< Angle between particle and the exchange
                               ///< boson in the hadron beam rest frame
      Double32_t ptVsGamma;    ///< pt w.r.t. the virtual photon in the
                               ///< hadron beam rest frame
      Double32_t phiPrf;       ///< Azimuthal angle around virtual
                               ///< photon in hadron beam rest frame
      
      TRef event; ///< Persistent reference to the event containing
                  ///< this particle.
      
      ClassDef(ParticleMC, 1)
   };
   
   
   inline UInt_t ParticleMC::Index() const { return I; }
   
   inline UShort_t ParticleMC::Status() const { return KS; }
   
   inline UShort_t ParticleMC::ParentIndex() const { return orig; }
   
   inline UShort_t ParticleMC::Child1() const { return daughter; }
   
   inline UShort_t ParticleMC::ChildN() const { return ldaughter; }
   
   inline Double_t ParticleMC::Px() const { return px; }
   
   inline Double_t ParticleMC::Py() const { return py; }
   
   inline Double_t ParticleMC::Pz() const { return pz; }
   
   inline Double_t ParticleMC::M() const { return m; }
   
   inline Double_t ParticleMC::Pt() const { return pt; }
   
   inline TVector3 ParticleMC::Vertex() const {
      return TVector3(xv, yv, zv);
   }
   
   inline erhic::Pid ParticleMC::ParentId() const {
      return erhic::Pid(parentId);
   }
   
   inline Double_t ParticleMC::P() const { return p; }
   
   inline Double_t ParticleMC::Theta() const { return theta; }
   
   inline Double_t ParticleMC::Phi() const { return phi; }
   
   inline Double_t ParticleMC::Y() const { return rapidity; }
   
   inline Double_t ParticleMC::Eta() const { return eta; }
   
   inline Double_t ParticleMC::Z() const { return z; }
   
   inline Double_t ParticleMC::Xf() const { return xFeynman; }
   
   inline Double_t ParticleMC::AngleVsBoson() const { return thetaGamma; }
   
   inline Double_t ParticleMC::PtVsBoson() const { return ptVsGamma; }
   
   inline Pid ParticleMC::Id() const { return Pid(id); }
   
   inline UInt_t ParticleMC::NChildren() const {
      if(0 == daughter) return 0;
      if(0 == ldaughter) return 1;
      return ldaughter - daughter + 1;
   }
   
   inline Double_t ParticleMC::GetM() const { return M(); }
   
   inline Double_t ParticleMC::GetE() const { return E; }
   
//   inline void ParticleMC::SetE(Double_t e) { E = e; }
   
//   inline void ParticleMC::SetP(Double_t momentum) { p = momentum; }
   
//   inline void ParticleMC::SetPt(Double_t momentum) { pt = momentum; }
   
//   inline void ParticleMC::SetPz(Double_t momentum) { pz = momentum; }
   
//   inline void ParticleMC::SetPhi(Double_t value) { phi = value; }
   
//   inline void ParticleMC::SetTheta(Double_t value) { theta = value; }
   
//   inline void ParticleMC::SetId(Int_t theId) { id = theId; }
   
//   inline void ParticleMC::SetStatus(UShort_t status) { KS = status; }
   
} // namespace erhic

#endif
