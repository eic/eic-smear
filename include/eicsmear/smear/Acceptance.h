/**
 \file
 Declaration of class Smear::Acceptance.
 
 \author    Michael Savastio
 \date      2011-08-19
 \copyright 2011 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_SMEAR_ACCEPTANCE_H_
#define INCLUDE_EICSMEAR_SMEAR_ACCEPTANCE_H_

#include <set>
#include <vector>

#include <TFormula.h>
#include <TMath.h>

#include "eicsmear/smear/Smear.h"  // Definition of KinType

class TString;

namespace erhic {

class VirtualParticle;

}  // namespace erhic

namespace Smear {

/**
 Defines a range of acceptance in one or more of: theta, phi, E, p, pt, pz.
 Comprises one or more acceptance zones, which may or may not overlap.
 \todo Implement data hiding
 */
class Acceptance {
 public:
  /**
   A (min, max) range in some variable evaluated as an arbitrary function
   of theta, phi, E and p (up to two supported).
   For example, if you want to set the acceptance in pT to [0.,100.]
   CustomCut("P*sin(theta)", 0., 100.);
   */
  class CustomCut {
   public:
    virtual ~CustomCut();
    CustomCut();
    CustomCut(const TString&, double min, double max);
    virtual bool Contains(const erhic::VirtualParticle&) const;
   protected:
    TFormula mFormula;
    int dim;
    KinType Kin1;
    KinType Kin2;
    double Min;
    double Max;
  };

  /**
   A single contiguous region of acceptance.
   */
  class Zone {
   public:
    /** Destructor */
    virtual ~Zone();

    /**
     Constructor.
     Define accepted ranges in theta, phi, E, p, pT and pz.
     Ranges in each variable are combined via boolean AND.
     By default accepts all particles.
     */
    Zone(double theta = 0., double = TMath::Pi(),
         double phi = 0., double = TMath::TwoPi(),
         double E = 0., double = TMath::Infinity(),
         double p = 0., double = TMath::Infinity(),
         double pt = 0., double = TMath::Infinity(),
         double pz = -TMath::Infinity(), double = TMath::Infinity());

    /**
     Add a CustomCut to the list of acceptance tests.
     */
    virtual void Add(const CustomCut&);

    /**
     Returns true if the particle lies in this zone, false if not.
     */
    virtual Bool_t Contains(const erhic::VirtualParticle&) const;

   protected:
    double thetaMin;
    double thetaMax;
    double phiMin;
    double phiMax;
    double EMin;
    double EMax;
    double PMin;
    double PMax;
    double pTMin;
    double pTMax;
    double pZMin;
    double pZMax;
    std::vector<Smear::Acceptance::CustomCut> CustomCuts;
    // We want to be able to write Acceptance objects to a ROOT file,
    // so nested classes need a dictionary as well.
    ClassDef(Smear::Acceptance::Zone, 1)
  };

  /** Destructor */
  virtual ~Acceptance();

  /**
   Default constructor.
   The genre sets which particle types are accepted (see Smear::EGenre).
   By default, the device has 4pi coverage,
   and accepts particles with all energy and momenta.
   */
  explicit Acceptance(int genre = kAll);

  /**
   Add a new zone with user-specified coverage.
   Particles will be accepted if they fall within any acceptance zone.
   */
  void AddZone(const Zone&);

  /**
   Returns the number of acceptance zones.
   */
  UInt_t GetNZones() const;

  /**
   Returns the "genre" of the particle (em, hadronic or all).
   */
  Int_t GetGenre() const;

  /**
   Select the class(es) of particles to accept.
   */
  void SetGenre(int genre);

  /**
   Select the charges of particles to accept.
   */
  void SetCharge(ECharge charge);

  /**
   Returns the charge of particles to accept.
   */
  ECharge GetCharge() const;

  /**
   Add a particle type to the list of particles to be smeared.
   If you never add anything, the device will
   smear all stable particles within its acceptance.
   Must use PDG particle codes.
   */
  void AddParticle(int particle);

  /**
   This function determines if the particle provided lies within
   the acceptance of the
   detector.  Default acceptance is a full 4*pi solid angle,
   with E and p in (0.,1e12) GeV.
   This function automatically fixes polar and azimuthal angles
   which are not within their proper range.
   */
  bool Is(const erhic::VirtualParticle& prt) const;

 protected:
  int mGenre;
  ECharge mCharge;  // Particle charges accepted (neutral, charged or all)
  std::vector<Zone> mZones;
  std::set<int> mParticles;

  ClassDef(Smear::Acceptance, 1)
};

inline UInt_t Acceptance::GetNZones() const {
  return mZones.size();
}

inline Int_t Acceptance::GetGenre() const {
  return mGenre;
}

inline ECharge Acceptance::GetCharge() const {
  return mCharge;
}

}  // namespace Smear

#endif  // INCLUDE_EICSMEAR_SMEAR_ACCEPTANCE_H_
