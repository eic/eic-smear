/**
 \file
 Declaration of class erhic::ParticleMC.

 \author    Thomas Burton
 \date      2011-10-10
 \copyright 2011 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_ERHIC_PARTICLEMC_H_
#define INCLUDE_EICSMEAR_ERHIC_PARTICLEMC_H_

#include <string>

#include <TLorentzVector.h>
#include <TRef.h>
#include <TVector3.h>

#include "eicsmear/erhic/Pid.h"
#include "eicsmear/erhic/VirtualParticle.h"

namespace erhic {

class EventMC;
/**
 A particle produced by a Monte Carlo generator.
 */
class ParticleMC : public VirtualParticle {
 public:
  /**
   Default constructor.
   Optionally pass a string with particle information in a HEPEVT
   format, namely:
     "index status id parent firstChild lastChild px py pz E m xv yv zv"
   */
  explicit ParticleMC(const std::string& = "");

  /**
   Destructor
   */
  virtual ~ParticleMC();

  /**
   Print the contents of Particle to standard output.
   The format is that of the input Monte Carlo i.e.
   I KS id orig daughter ldaughter px py pz E m xv yv zv.
   Inherited from TObject. The argument is unused.
   */
  virtual void Print(Option_t* = "") const;

  /**
   Returns the particle index in an event, in the range [1, N].
   */
  virtual UInt_t GetIndex() const;

  /**
   Returns the status of the particle.
   The meaning of the status code depends on the generator.
   For PYTHIA, see the description of variable K(I,1) in the manual.
   */
  virtual UShort_t GetStatus() const;

  /**
   Returns the index of this particle's parent in an event.
   */
  virtual UShort_t GetParentIndex() const;

  /**
   Returns a pointer to the parent of this particle.
   This is the particle with index GetParentIndex() in the
   event containing this particle (obtainable via GetEvent()).
   Returns NULL if this particle has no parent or if it cannot
   be accessed via GetEvent().
   */
  virtual const ParticleMC* GetParent() const;

  /**
   Returns the index of this particle's first child particle.
   Returns 0 if this particle has no children.
   */
  virtual UShort_t GetChild1Index() const;

  /**
   Returns the index of this particle's last child particle.
   Returns 0 if this particle has zero or one children.
   */
  virtual UShort_t GetChildNIndex() const;

  /**
   Returns the number of children of this particle.
   Returns 0 if the particle did not decay.
   */
  virtual UInt_t GetNChildren() const;

  /**
   Returns a pointer to the nth child particle of this particle,
   where n is in the range [0, GetNChildren()).
   GetChild(0) returns the child whose index in this particle's
   event is GetChild1Index().
   GetChild(GetNChildren()-1) returns the child whose index
   is GetChildNIndex().
   Returns NULL if there is no such child particle or it
   cannot be accessed via the event for some reason (see GetEvent()).
   */
  virtual const ParticleMC* GetChild(UShort_t) const;

  /**
   Returns true if n in the range [0, N), where N is the number
   of children of this particle.
   Returns false otherwise.
   Equivalent to GetChild(n) != NULL.
   */
  virtual Bool_t HasChild(Int_t) const;

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
   Returns invariant mass (GeV/c<sup>2</sup>).
   */
  virtual Double_t GetM() const;

  /**
   Returns momentum transverse to the beam direction.
   */
  virtual Double_t GetPt() const;

  /**
   Returns the origin point of the particle (cm).
   (0,0,0) indicates a particle originating in the collision.
   */
  virtual TVector3 GetVertex() const;

  /**
   Returns the identity information of this particle's parent.
   
   */
  virtual Pid GetParentId() const;

  /**
   Returns the total momentum (GeV).
   */
  virtual Double_t GetP() const;

  /**
   Returns the polar angle in the range [0, pi] radians.
   */
  virtual Double_t GetTheta() const;

  /**
   Returns the polar angle in the range [0, 2pi] radians.
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

  /**
   Returns the variable z.
   z = (P.p_h)/(P.q).
   */
  virtual Double_t GetZ() const;

  /**
   Returns Feynman-x.
   x<sub>F</sub> = 2*p<sub>z</sub>/sqrt(s).
   */
  virtual Double_t GetXFeynman() const;

  /**
   Returns the angle with respect to the exchange boson.
   Defined in the beam hadron's rest frame.
   Given in the range [0,pi] radians.
   */
  virtual Double_t GetThetaVsGamma() const;

  /**
   Returns the p<sub>T</sub> with respect to the exchange boson.
   Defined in the beam hadron's rest frame.
   */
  virtual Double_t GetPtVsGamma() const;

  /**
   Returns a pointer to the event containing this particle.
   Returns NULL if this particle has not been associated
   with an event.
   */
  const EventMC* GetEvent() const;

  /**
   Set the event with which to associate this particle.
   */
  void SetEvent(EventMC* event);

  /**
   Returns the (E,p) 4-vector in the lab frame.
   */
  virtual TLorentzVector Get4Vector() const;

  /**
   Returns the (E,p) 4-vector in the lab frame.
   */
  virtual TLorentzVector PxPyPzE() const { return Get4Vector(); }

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
   
   \note Due to details of the implementation and how ROOT handles reading
         events in a TTree, this function will not work for the exchange
         boson when using a TTree::Draw (or similar) statement. It does
         work for the exchange boson if reading events manually via
         TTree::GetEntry.
   */
  virtual TLorentzVector Get4VectorInHadronBosonFrame() const;

  /**
   Returns the energy of the particle in the lab frame.
   */
  virtual Double_t GetE() const;

  virtual void SetE(Double_t);

  virtual void SetM(Double_t);

  virtual void SetP(Double_t);

  virtual void SetPt(Double_t);

  virtual void SetPz(Double_t);

  virtual void SetPhi(Double_t);

  virtual void SetTheta(Double_t);

  virtual void SetStatus(UShort_t);

  /**
   Returns the ID of the particle.
   */
  virtual Pid Id() const;

  /**
   Returns the ID of the particle.
   */
  virtual Pid GetPdgCode() const { return Id(); }

  /**
   Sets quantities derived from the four-momentum (E, px, py, pz), namely
   <UL>
   <LI> total momentum </LI>
   <LI> transverse momentum </LI>
   <LI> rapidity </LI>
   <LI> pseudorapidity </LI>
   <LI> theta </LI>
   <LI> phi </LI>
   </UL>
   This should be called if (E, px, py, pz) are manually altered
   in order to propagate the changes to these other quantities.
   */
  virtual void ComputeDerivedQuantities();

  /**
   Sets quantities that depend on the properties of the event or
   associations of one particle with another, namely
   <UL>
   <LI> z </LI>
   <LI> Feynman x </LI>
   <LI> angle and pt with respect to the exchanged boson </LI> 
   <LI> azimuthal angle around the exchanged boson </LI>
   <LI> parent pdg code </LI>
   </UL>
   Important: this particle is assumed to be in the same frame of
   reference as those contained in the event that is passed as an
   argument.
   */
  virtual void ComputeEventDependentQuantities(EventMC&);

  /**
   Sets the index of the particle i.e. its position in the track list
   (in principle this can be any
   integer you require to associated with the particle).
   */
  virtual void SetIndex(int i) { I = i; }

  /** Sets the status code of the particle (generally final state
   particles are given status == 1 */
  virtual void SetStatus(int i) { KS = i; }

  /** Sets the ID of the particle. In order to make use of class Pid this
   should be the PDG code of the particle, but in principle can be any
   value you wish to use to identify it. */
  virtual void SetId(int i) { id = i; }

  /** Sets the index of this particle's parent if it has one.
   By default this is zero, indicating no parent. */
  virtual void SetParentIndex(int i) { orig = i; }

  /** Sets the index of this particle's first child. By default this
   is zero, indicating no children. */
  virtual void SetChild1Index(int i) { daughter = i; }

  /** Sets the index of this particle's last child.
   By default this is zero, indication zero or one children. */
  virtual void SetChildNIndex(int i) { ldaughter = i; }

  /**
   Sets the four-momentum of the particle.
   Changes are propagated to derived quantities.
   */
  virtual void Set4Vector(const TLorentzVector&);

  /**
   Sets the origin coordinates
   */
  virtual void SetVertex(const TVector3&);

  /**
   Sets the ID of this particle's parent. See comments in SetId()
   */
  virtual void SetParentId(int i) { parentId = i; }

// protected:

  UShort_t    I;             ///< Particle index in event
  UShort_t    KS;            ///< Particle status code: see PYTHIA manual
  Int_t       id;            ///< PDG particle code
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

  TRef event;  ///< Persistent reference to the event containing
               ///< this particle.

  ClassDef(ParticleMC, 1)
};

inline UInt_t ParticleMC::GetIndex() const {
  return I;
}

inline UShort_t ParticleMC::GetStatus() const {
  return KS;
}

inline UShort_t ParticleMC::GetParentIndex() const {
  return orig;
}

inline UShort_t ParticleMC::GetChild1Index() const {
  return daughter;
}

inline UShort_t ParticleMC::GetChildNIndex() const {
  return ldaughter;
}

inline Double_t ParticleMC::GetPx() const {
  return px;
}

inline Double_t ParticleMC::GetPy() const {
  return py;
}

inline Double_t ParticleMC::GetPz() const {
  return pz;
}

inline Double_t ParticleMC::GetM() const {
  return m;
}

inline Double_t ParticleMC::GetPt() const {
  return pt;
}

inline TVector3 ParticleMC::GetVertex() const {
  return TVector3(xv, yv, zv);
}

inline erhic::Pid ParticleMC::GetParentId() const {
  return erhic::Pid(parentId);
}

inline Double_t ParticleMC::GetP() const {
  return p;
}

inline Double_t ParticleMC::GetTheta() const {
  return theta;
}

inline Double_t ParticleMC::GetPhi() const {
  return phi;
}

inline Double_t ParticleMC::GetRapidity() const {
  return rapidity;
}

inline Double_t ParticleMC::GetEta() const {
  return eta;
}

inline Double_t ParticleMC::GetZ() const {
  return z;
}

inline Double_t ParticleMC::GetXFeynman() const {
  return xFeynman;
}

inline Double_t ParticleMC::GetThetaVsGamma() const {
  return thetaGamma;
}

inline Double_t ParticleMC::GetPtVsGamma() const {
  return ptVsGamma;
}

inline Pid ParticleMC::Id() const {
  return Pid(id);
}

inline UInt_t ParticleMC::GetNChildren() const {
  if (0 == daughter) return 0;
  if (0 == ldaughter) return 1;
  return ldaughter - daughter + 1;
}

inline Double_t ParticleMC::GetE() const {
  return E;
}

inline void ParticleMC::SetE(Double_t e) {
  E = e;
}

inline void ParticleMC::SetM(Double_t mass) {
  m = mass;
}

inline void ParticleMC::SetP(Double_t momentum) {
  p = momentum;
}

inline void ParticleMC::SetPt(Double_t momentum) {
  pt = momentum;
}

inline void ParticleMC::SetPz(Double_t momentum) {
  pz = momentum;
}

inline void ParticleMC::SetPhi(Double_t value) {
  phi = value;
}

inline void ParticleMC::SetTheta(Double_t value) {
  theta = value;
}

inline void ParticleMC::SetStatus(UShort_t status) {
  KS = status;
}

}  // namespace erhic

#endif  // INCLUDE_EICSMEAR_ERHIC_PARTICLEMC_H_
