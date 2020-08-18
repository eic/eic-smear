#ifndef __PID_H__
#define __PID_H__
	
//
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
//   -- numSigma(double eta, double p,        PID::type PID);
//   -- maxP    (double eta, double numSigma, PID::type PID);
//   -- minP    (double eta, double numSigma, PID::type PID);
//   -- name    ();
//
//  Here PID::type is an enumerated constant set allowing one to choose pi-vs-k or k-vs-p etc...
//
//  The detector types that inherit from PID will clearly have parameters that define their 
//  performance and these are expected to be arguments of the constructor of those derived classes.
//  For example:
//     PID* tof = new tofBarrel(radius, etaLow, etaHigh, sigmaT);
//
//

#include <string>
	
class PID
{
public:
  PID() {}
  virtual ~PID() {}
	
  enum type {pi_k , k_p};

  virtual bool   valid    (double eta, double p                      )=0;
  virtual double numSigma (double eta, double p,        PID::type PID)=0;
  virtual double maxP     (double eta, double numSigma, PID::type PID)=0;
  virtual double minP     (double eta, double numSigma, PID::type PID)=0;
  virtual std::string name()=0;
  virtual void description()=0;
		
protected:
	
};
	
#endif /* __PID_H__ */
