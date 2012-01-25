#ifndef _ERHIC_BUILDTREE_Pid_H_
#define _ERHIC_BUILDTREE_Pid_H_

// Pid.h
//
// Created by TB on 8/12/11.
// Copyright 2011 BNL. All rights reserved.

#include <Rtypes.h>
#include <TObject.h>

class TParticlePDG;

namespace erhic {
   
   /**
    Particle identity.
    */
   class Pid {
      
   public:
      
      /**
       Returns a code indicating 'not valid'.
       Use this value to consistently indicate particles without a PDG code.
       */
      static Int_t InvalidCode();
      
      /**
       Constructor.
       */
      Pid(Int_t = Pid::InvalidCode());
      
      virtual ~Pid();
      
      /**
       Returns the integer PDG code.
       */
      Int_t Code() const;
      
      /**
       Sets the integer PDG code.
       */
      void Set(Int_t);
      
      /**
       Returns the particle information object corresponding to this PDG code.
       From this, properties such as particle class, name, charge and mass
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
       Compares the argument with this objects integer PDG code.
       */
      bool operator==(Int_t i) const;
      
      /**
       Comparison operator with integer.
       Compares the argument with this objects integer PDG code.
       */
      bool operator<(Int_t i) const;
      
      /**
       Comparison operator with integer.
       Compares the argument with this objects integer PDG code.
       */
      bool operator>(Int_t i) const;
      
      /**
       Comparison operator with integer.
       Compares the argument with this objects integer PDG code.
       */
      bool operator!=(Int_t i) const;
      
      /**
       Comparison operator with integer.
       Compares the argument with this objects integer PDG code.
       */
      bool operator<=(Int_t i) const;
      
      /**
       Comparison operator with integer.
       Compares the argument with this objects integer PDG code.
       */
      bool operator>=(Int_t i) const;
      
   protected:
      
      Int_t mCode; ///< The integer PDG code
      
      ClassDef(Pid, 1)
   };
   
   inline Pid::operator Int_t() const { return mCode; }
   
   inline Int_t Pid::Code() const { return mCode; }
   
   inline void Pid::Set(Int_t code) { mCode = code; }
   
   inline bool Pid::operator==(Int_t i) const { return mCode == i; }
   
   inline bool Pid::operator<(Int_t i) const { return mCode < i; }
   
   inline bool Pid::operator>(Int_t i) const { return mCode > i; }
   
   inline bool Pid::operator!=(Int_t i) const { return ! operator==(i); }
   
   inline bool Pid::operator<=(Int_t i) const { return mCode <= i; }
   
   inline bool Pid::operator>=(Int_t i) const { return mCode >= i; }
   
} // namespace

#endif
