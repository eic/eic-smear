/**
 \file
 Declaration of class Smear::PerfectId.
 
 \author    Michael Savastio
 \date      2011-08-19
 \copyright 2011 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_SMEAR_PERFECTID_H_
#define INCLUDE_EICSMEAR_SMEAR_PERFECTID_H_

#include <set>
#include <vector>

#include "eicsmear/smear/Smearer.h"

namespace Smear {

/**
 Smearer that copies the PDG ID of a particle to a smeared particle
 with no modification.
 Use this to represent perfect particle identification performance,
 or if you don't care about modelling ID performance but want PID copied to
 the output smeared tree.
 */
class PerfectID : public Smearer {
 public:
  /**
   Constructor
   @param pdg [in] An optional list of PDG codes with which the PerfectID
   works. The ID will not be copied for particles with other PDG codes.
   If the list is empty the PerfectID operates on all PDG codes.
   */
  PerfectID(const std::vector<Int_t>& pdg = std::vector<int>());

  /**
   Destructor.
   */
  virtual ~PerfectID();

  /**
   Returns a new copy of this object.
   The argument has no effect.
   */
  virtual PerfectID* Clone(const char* = "") const;

  /**
   Copies the PDG code from the ParticleMC to the ParticleMCS if either
   the PDG code is in the initialisation list or if the list was empty.
   Otherwise, does nothing.
   */
  virtual void Smear(const erhic::VirtualParticle&, ParticleMCS&);

  /**
   Prints information about this object.
   */
  virtual void Print(Option_t* = "") const;

  /**
   Add a PDG code to the list.
   */
  virtual void Insert(Int_t);

 protected:
  /// PDG codes to copy. Does not operate on particles with other codes.
  std::set<Int_t> mPdg;

  ClassDef(Smear::PerfectID, 1)
};

}  // namespace Smear

#endif  // INCLUDE_EICSMEAR_SMEAR_PERFECTID_H_
