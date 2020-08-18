/**
 \file
 Declaration of class Smear::ParticleMCS.

 \author    Michael Savastio
 \date      2011-08-19
 \copyright 2011 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_SMEAR_PARTICLEMCS_H_
#define INCLUDE_EICSMEAR_SMEAR_PARTICLEMCS_H_

#include <cmath>
#include <TLorentzVector.h>

#include "eicsmear/erhic/Pid.h"
#include "eicsmear/erhic/VirtualParticle.h"
#include "eicsmear/smear/SmearConstants.h"
#include <iostream>

namespace Smear {

class Event;

/**
 A smeared Monte Carlo particle.
 */
class ParticleMCS : public erhic::VirtualParticle {

 public:
  /**
   Destructor.
   */
  virtual ~ParticleMCS();

  /**
   Default constructor.
   Initialises the Particle from the argument string with the format
   I KS id orig daughter ldaughter px py pz m E xv yv zv
   */
  ParticleMCS();

  /**
   Construct from an E-p 4-vector, pdg code and status code.
   */
  ParticleMCS(const TLorentzVector&, int pdg, int status);

  // Let's kill all else
  ParticleMCS(const ParticleMCS&) =delete;
  ParticleMCS& operator=(const ParticleMCS&) =delete;
  

  // ---------------
  // --- Getters ---
  // ---------------
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
   Returns the apparent mass of the smeared particle.
   \todo Consider the implementation here.
   If the particle is identified (either correctly or incorrectly),
   it could return the PDG mass of that particle.
   If the particle is not identified, should it return sqrt(E^2 - p^2)?
   What about if p or E are not known? Return E or p? Or zero (e.g.
   that would work for photons, E known, p not).
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
   Returns the pdg ID of the particle.
   */
  virtual ::erhic::Pid Id() const;

  /**
     Returns the numSigma deviation of this particle
  */
  double GetNumSigma() const;

  /**
     Returns the numSigma deviation type (pi_k, k_p) of this particle
  */
  int GetNumSigmaType() const;

  /** should always be true for a ParticleMCS
      This replaces the brittle mechanism of checking values against 0
      If false, it should indicate an old tree was used.
  */
  virtual bool IsSmeared() const;
  virtual bool IsESmeared() const; //< E smeared?
  virtual bool IsPSmeared() const; //< Total momentum smeared?
  virtual bool IsPtSmeared() const; //< P_t smeared?
  virtual bool IsPxSmeared() const; //< P_x smeared?
  virtual bool IsPySmeared() const; //< P_y smeared?
  virtual bool IsPzSmeared() const; //< P_z smeared?
  virtual bool IsThetaSmeared() const; //< &theta; smeared?
  virtual bool IsPhiSmeared() const; //< &phi; smeared?
  virtual bool IsIdSmeared() const; //< pdg Id smeared?
  virtual bool IsNumSigmaSmeared() const; //< PID numSigma smeared?

  // ---------------
  // --- Setters ---
  // ---------------
  /** Set energy. By default, marks it as smeared and checks that it wasn't before.
      @param CheckSetSmearFlag=false disables this check (e.g. for adjustments)
  */
  virtual void SetE(const Double_t value, const bool CheckSetSmearFlag=true);

  /** Set total momentum P. By default, marks it as smeared and checks that it wasn't before.
      @param CheckSetSmearFlag=false disables this check (e.g. for adjustments)
  */
  virtual void SetP(const Double_t value, const bool CheckSetSmearFlag=true);

  /** Set transverse momentum Pt. By default, marks it as smeared and checks that it wasn't before.
      @param CheckSetSmearFlag=false disables this check (e.g. for adjustments)
  */
  virtual void SetPt(const Double_t value, const bool CheckSetSmearFlag=true);

  /** Set P_x. By default, marks it as smeared and checks that it wasn't before.
      @param CheckSetSmearFlag=false disables this check (e.g. for adjustments)
  */
  virtual void SetPx(const Double_t value, const bool CheckSetSmearFlag=true);

  /** Set P_y. By default, marks it as smeared and checks that it wasn't before.
      @param CheckSetSmearFlag=false disables this check (e.g. for adjustments)
  */
  virtual void SetPy(const Double_t value, const bool CheckSetSmearFlag=true);

  /** Set P_z. By default, marks it as smeared and checks that it wasn't before.
      @param CheckSetSmearFlag=false disables this check (e.g. for adjustments)
  */
  virtual void SetPz(const Double_t value, const bool CheckSetSmearFlag=true);

  /** Set azimuth &phi;. By default, marks it as smeared and checks that it wasn't before.
      @param CheckSetSmearFlag=false disables this check (e.g. for adjustments)
  */
  virtual void SetPhi(const Double_t value, const bool CheckSetSmearFlag=true);

  /** Set polar angle &theta;. By default, marks it as smeared and checks that it wasn't before.
      @param CheckSetSmearFlag=false disables this check (e.g. for adjustments)
  */
  virtual void SetTheta(const Double_t value, const bool CheckSetSmearFlag=true);

  /** Set particle id (pdg code). By default, marks it as smeared and checks that it wasn't before.
      @param CheckSetSmearFlag=false disables this check (e.g. for adjustments)
  */
  virtual void SetId(Int_t value, const bool CheckSetSmearFlag=true);

  /** Set the numSigma deviation of this particle
  */
  virtual void  SetNumSigma( const double d, const bool CheckSetSmearFlag=true);

  /** Set the numSigma deviation type (pi_k, k_p) of this particle
  */
  virtual void  SetNumSigmaType( const int i);

  virtual void SetSmeared( bool flag=true);  //< Particle smeared
  virtual void SetESmeared( bool flag=true); //< E smeared
  virtual void SetPSmeared( bool flag=true); //< Total momentum smeared
  virtual void SetPtSmeared( bool flag=true); //< P_t smeared
  virtual void SetPxSmeared( bool flag=true); //< P_x smeared
  virtual void SetPySmeared( bool flag=true); //< P_y smeared
  virtual void SetPzSmeared( bool flag=true); //< P_z smeared
  virtual void SetThetaSmeared( bool flag=true); //< &theta; smeared
  virtual void SetPhiSmeared( bool flag=true); //< &phi; smeared
  virtual void SetIdSmeared( bool flag=true); //< pdg Id smeared
  virtual void SetNumSigmaSmeared( bool flag=true); //< PID numSigma smeared
  

  /**
   Dummy one; just need to compile;
   */
  void Set4Vector(const TLorentzVector&) { }

  virtual void SetStatus(Int_t);

  virtual UShort_t GetParentIndex() const { return 0; }

  /**
   Prints the attributes of this particle to standard output.
   The output format is "status id px py pz E".
   */
  virtual void Print(Option_t* = "") const;

  /**
   Sets the origin coordinates.
   No-op, just required for compilation.
   */
  virtual void SetVertex(const TVector3&) { }

  /**
     Stores z in the ParticleS.K where K is the kinematic variable associated with kin.
  */
  virtual void SetVariable(const double z, const KinType kin);

  /**
     This dictates how the class deals with positive definite variables
     which have been smeared to negative values.
  */
  void HandleBogusValues( const KinType kin );


 protected:

  /** should always be true for a ParticleMCS
      This replaces the brittle mechanism of checking values against 0
      If false, it should indicate an old tree was used.
   */
  bool kParticleSmeared=false;
  bool kESmeared=false;
  bool kPSmeared=false;
  bool kPtSmeared=false;
  bool kPxSmeared=false;
  bool kPySmeared=false;
  bool kPzSmeared=false;
  bool kThetaSmeared=false;
  bool kPhiSmeared=false;
  bool kIdSmeared=false;
  bool kNumSigmaSmeared=false;


  UShort_t   status;      ///< Status code
  Int_t      id;          ///< PDG particle code
  Double32_t px;          ///< x component of particle momentum
  Double32_t py;          ///< y component of particle momentum
  Double32_t pz;          ///< z component of particle momentum
  Double32_t E;           ///< Energy of particle
  Double32_t pt;          ///< Transverse momentum of particle
  Double32_t p;           ///< Total momentum of particle
  Double32_t theta;       ///< Polar angle
  Double32_t phi;         ///< Azimuthal angle

  double  numSigma;       ///< PID: nSigma deviation
  int numSigmaType;       ///< PID: nSigma deviation type. pi_k == 1, k_p == 2. \TODO: This should be agreed upon and fixed better.

  ClassDef(Smear::ParticleMCS, 3)
};


}  // namespace Smear

#endif  // INCLUDE_EICSMEAR_SMEAR_PARTICLEMCS_H_
