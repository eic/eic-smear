/**
 \file
 Declaration of class Smear::FormulaString.
 
 \author    Michael Savastio
 \date      2011-08-19
 \copyright 2011 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_SMEAR_FORMULASTRING_H_
#define INCLUDE_EICSMEAR_SMEAR_FORMULASTRING_H_

#include <string>
#include <vector>

#include <TObject.h>

#include "eicsmear/smear/Smear.h"  // For KinType enum

class TFormula;

namespace Smear {

/**
 A formula described by a string with up to four variables.
 */
class FormulaString : public TObject {
 public:
  /**
   Destructor.
   */
  virtual ~FormulaString();

  /**
   Default construtor.
   Only present to meet ROOT requirements for objects to have a default
   constructor.
   Always evaluates to zero.
   */
  FormulaString();

  /**
   Initialise with a string describing the formula.
   Formula syntax is as for ROOT::TFormula.
   Valid variable names are:
   <ul>
   <li>"E"</li>
   <li>"P"</li>
   <li>"theta"</li>
   <li>"phi"</li>
   <li>"pZ"</li>
   <li>"pT"</li>
   </ul>
   Variable names ARE case sensitive.
   Energy (momentum) are in GeV(/c), angles are in radians.
   e.g. for the function p/sin(theta)
   Smear::FormulaString("P/sin(theta)");
   */
  explicit FormulaString(const std::string&);

  /**
   Evaluate the formula with the provided arguments.
   Arguments should be listed in the order they were
   named in the constructor string.
   e.g. if the input function was "P/sin(theta)" pass a vector with
   element 0 == 0.5 and element 1 == 3 to evaluate at P == 0.5 and
   theta == 3.
   */
  virtual double Eval(const std::vector<double>&) const;

  /**
   Returns a vector of Smear::KinType corresponding to the variables
   named in the constructor string.
   e.g. If initialised with "P/sin(theta)", would result in a vector with
   two elements: kP, kTheta.
   */
  virtual std::vector<Smear::KinType> Variables() const;

  /**
   Returns the processed formula string with variables substituted.
   */
  virtual std::string GetString() const;

  /**
   Returns the unprocessed input formula string.
   */
  virtual std::string GetInputString() const;

  /**
   Returns the name corresponding the a Smear::KinType.
   */
  static std::string GetKinName(KinType);

  /**
   Returns the KinType corresponding to the variable name,
   or kInvalidKinType if the name is invalid.
   */
  static KinType GetKinType(const std::string&);

 protected:
  /**
   Process the input string, containing "P", "theta" etc into a version
   with "x", "y", "z", "t" substituted, compatible with TFormula.
   */
  std::string Parse(const std::string&);

  TFormula* mFormula;
  std::string mInput;  ///< Original formula (before parsing)
  std::vector<Smear::KinType> mVariables;

  ClassDef(Smear::FormulaString, 1)
};

}  // namespace Smear

#endif  // INCLUDE_EICSMEAR_SMEAR_FORMULASTRING_H_
