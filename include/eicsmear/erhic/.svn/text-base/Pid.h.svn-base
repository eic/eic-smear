/**
 \file
 Declaration of class erhic::Pid.
 
 \author    Thomas Burton
 \date      2011-08-12
 \copyright 2011 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_ERHIC_PID_H_
#define INCLUDE_EICSMEAR_ERHIC_PID_H_

#include <Rtypes.h>

class TParticlePDG;

namespace erhic {

/**
 Particle identity.
 */
class Pid {
 public:
  /**
   Returns a code indicating 'not valid'.
   Use this value to consistently indicate particles without an ID number.
   */
  static Int_t InvalidCode();

  /**
   Constructor.
   @param [in] pdg The particle's PDG code, an integer
   */
  Pid(Int_t pdg = Pid::InvalidCode());

  /** Destructor */
  virtual ~Pid();

  /** Returns the PDG code, an integer. */
  Int_t Code() const;

  /**
   Sets the integer PDG code.
   @param [in] pdg The particle's PDG code, an integer
   */
  void Set(Int_t pdg);

  /**
   Returns the particle information object corresponding to this PDG code.
   From this, properties such as particle name, charge and mass
   can be accessed.
   */
  TParticlePDG* Info() const;

  /**
   Type conversion operator to integer PDG code.
   int i = Pid(211) is equivalent to int i = Pid(211).Code();
   */
  operator Int_t() const;

  /**
   Comparison operator with integer.
   Compares the argument with this object's integer PDG code.
   */
  bool operator==(Int_t i) const;

  /**
   Comparison operator with integer.
   Compares the argument with this object's integer PDG code.
   */
  bool operator<(Int_t i) const;

  /**
   Comparison operator with integer.
   Compares the argument with this object's integer PDG code.
   */
  bool operator>(Int_t i) const;

  /**
   Comparison operator with integer.
   Compares the argument with this object's integer PDG code.
   */
  bool operator!=(Int_t i) const;

  /**
   Comparison operator with integer.
   Compares the argument with this object's integer PDG code.
   */
  bool operator<=(Int_t i) const;

  /**
   Comparison operator with integer.
   Compares the argument with this object's integer PDG code.
   */
  bool operator>=(Int_t i) const;

 protected:
  Int_t mCode;  ///< The integer PDG code

  ClassDef(erhic::Pid, 1)
};

inline Pid::operator Int_t() const { return mCode; }

inline Int_t Pid::Code() const { return mCode; }

inline void Pid::Set(Int_t code) { mCode = code; }

inline bool Pid::operator==(Int_t i) const { return mCode == i; }

inline bool Pid::operator<(Int_t i) const { return mCode < i; }

inline bool Pid::operator>(Int_t i) const { return mCode > i; }

inline bool Pid::operator!=(Int_t i) const { return !operator==(i); }

inline bool Pid::operator<=(Int_t i) const { return mCode <= i; }

inline bool Pid::operator>=(Int_t i) const { return mCode >= i; }

}  // namespace erhic

#endif  // INCLUDE_EICSMEAR_ERHIC_PID_H_
