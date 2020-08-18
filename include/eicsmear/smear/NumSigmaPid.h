/**
 \file
 Declaration of base class Smear::numSigmaPid.
 
 \author    Roberto Preghenella, K. Kauder
 \date      2020-08-17
 */

#ifndef INCLUDE_EICSMEAR_SMEAR_NUMSIGMAPID_H
#define INCLUDE_EICSMEAR_SMEAR_NUMSIGMAPID_H

#include "eicsmear/smear/Smearer.h"
#include "PID.h"

// Based on https://gitlab.com/preghenella/pid

// A smearer needs to implement
//   virtual void Smear(const erhic::VirtualParticle&, ParticleMCS&) = 0;
// And it has an Acceptance object. 



//  Hello PID Fans:
//
//  This is the base class for a simple set of objects that 
//  are used to evaluate various detector technology choices for the EIC.
//
//  The idea is simple.  Regardless of PID detector technology choice 
//  one must be able to answer some simple questions.  This virtual base class
//  organizes the method of asking and answering these simple questions to
//  allow for efficient and apples-to-apples comparisons.
//  
//  The base class requires that any derived class procide responses to several queries:
//   -- valid   (double eta, double p                      );
//   -- numSigma(double eta, double p,        numSigmaPid::type PID);
//   -- maxP    (double eta, double numSigma, numSigmaPid::type PID);
//   -- minP    (double eta, double numSigma, numSigmaPid::type PID);
//   -- name    ();
//
//  Here numSigmaPid::type is an enumerated constant set allowing one to choose pi-vs-k or k-vs-p etc...
//
//  The detector types that inherit from numSigmaPid will clearly have parameters that define their 
//  performance and these are expected to be arguments of the constructor of those derived classes.
//  For example:
//     numSigmaPid* tof = new tofBarrel(radius, etaLow, etaHigh, sigmaT);
//
//

#include <string>
#include <memory>
	
namespace Smear {

  /**  Base class for a simple set of objects that 
   *  are used to evaluate various PID detector technology choices for the EIC.
   */
  class NumSigmaPid : public Smearer {
  public:
    /**
       Default constructor.
    */
    NumSigmaPid();

    /**
       Destructor.
    */
    virtual ~NumSigmaPid() {}
    
    /**
       Smears the input ParticleMC and stores the result(s) in the ParticleMCS.
    */
    void Smear(const erhic::VirtualParticle&, ParticleMCS&);

    /** Set the numSigma type
	PID::type is an enumerated constant set 
	allowing one to choose pi-vs-k or k-vs-p etc.
    */
    void SetNumSigmaType( const int i );

    /** Get the numSigma type
	PID::type is an enumerated constant set 
	allowing one to choose pi-vs-k or k-vs-p etc.
     */
    int GetNumSigmaType( ) const;

    /** Needed for TObject
     */
    NumSigmaPid* Clone(const char* = "") const;

    /** Passing through methods of PID
     */
    bool valid  (double eta, double p ) const;

    /** Passing through methods of PID
     */
    double maxP (double eta, double numSigma, PID::type PID) const;

    /** Passing through methods of PID
     */
    double minP (double eta, double numSigma, PID::type PID) const;

    /** Passing through methods of PID
     */
    std::string name() const;

    /** Passing through methods of PID
     */
    void description() const;

  protected:
    std::shared_ptr<PID> ThePidObject;
    int NumSigmaType=-1;
    PID::type EnumType;
  };

}
#endif // INCLUDE_EICSMEAR_SMEAR_NUMSIGMAPID_H
  
