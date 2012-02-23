#ifndef _ERHIC_BUILDTREE_DEVICE_
#define _ERHIC_BUILDTREE_DEVICE_

/*
 Change log
 
 2011/08/15
 Added .cxx file with implementation.
 Moved inline function definitions out of class declaration.
 Made Clone() const.
 Made Get...() methods const.
// Replaced TString arguments with const TString&.
 Replaced Particle arguments with const Particle&.
 Added arguments with default values to constructor.
 
 TODO Remove separate TF1, TF2 - have Param pointer allocated appropriately
 TODO depending on the number of arguments. Should be simple enough to extend
 TODO to TF3 also.
 */
#include <cmath>

#include <TF1.h>
#include <TF2.h>
#include <Math/ParamFunctor.h> // For ROOT::TMath::ParamFunctor
#include <TRandom3.h>
#include <TString.h>

#include "Acceptance.h"
#include "Smear.h"

class TRandom3;

namespace Smear {
   
	/**
	 The device class.
    This is an object that will smear a single particle-wise kinematic
	 variable according to a user provided parametrization.
    Devices can be added a detector
	 which can be used to smear data with a series of devices.
    Each device has an 
	 acceptance, which is an open subset of (E,p,theta,phi) momentum space.
    Only particles whose momentum lie in the acceptance can be smeared.
    \todo Add a kElectromagnetic/kHadronic enum of genre
    \todo Add genre as an optional third argument (default = all)
    \todo Inherit from TObject
    \todo Get rid of ParDim. Instead select TF1/2/3 via template?
	 */
	struct Device : public TObject {
		
		/**
		 Default constructor.  Sets a 1-dimensional parametrization f(x)=0 in
       range [0.,1.e6] by default.
       Default input kinematics E, theta, default variable to be smeared is E.
		 */
		Device(KinType = kE, TString parameterisation = "0.", int genre = 0);
		
		/**
		 This determines what sorts of particles the device is able to smear.
		 
		 0: (Default)The device doesn't care what the particle is, as long as
       it is stable, 
		 (i.e. in the Pythia sense), in the final state, and not a gluon.
		 1: The device detects final state photons and leptons only.
		 2: The device detects final state hadrons only.
		 */
		void SetGenre(int);
      
//		void SetGenre(const TString&);
		
		/**
		 Set the number of arguments of the smearing parametrization.  
		 Must be either 1 or 2, or the function call does nothing.
		 */
		void SetParDim(int);
		
		/**
		 Set the kinematics to be used as arguments in the smearing
       parametrization
		 (In1,In2), as well as the variable to be smeared (Out).
       Must use the 
		 KinType enumerators, kE, kP, kTheta, kPhi, kPz, kPt.
       Automatically sets the number
		 of arguments of the parametrization to be 2.
		 */
		void SetInOutKinematics(KinType In1, KinType In2, KinType Out);
		
		/**
		 Set the kinematics to be used as arguments in the smearing
       parametrization (In1)
		 as well as the variable to be smeared (Out).
       Must use the KinType enumerators 
		 kE, kP, kTheta, kPhi, kPz, kPt.
       Automatically sets the number of arguments of the 
		 parametrization to be 1.
		 */
		void SetInOutKinematics(KinType In, KinType Out);
		
		/** 
		 Set the kinematics to be used as arguments in the smearing
       parametrization.
		 Must use the KinType enumerators kE, kP, kTheta, kPhi, kPz, kPt.
       Automatically sets the 
		 number of arguments of the parametrization to be 2.
		 */
		void SetInputKinematics(KinType In1, KinType In2);
		
		/**
		 Set the kinematics to be used as an argument in the smearing
       parametrization.
		 Must use the KinType enumerators kE, kP, kTheta, kPhi, kPz, kPt.
       Automatically sets the 
		 number of arguments of the parametrization to be 1.
		 */
		void SetInputKinematics(KinType In); 
		
		/** 
		 Set the kinematic variable to be smeared.
       Must use the KinType enumerators
		 kE,kP,kTheta,kPhi,kPz,kPt.  
		 */
		virtual void SetSmearedKinematics(KinType Out);
		
		/**
		 Set the probability distribution function on which smeared values
       are to be generated.  By default, the distribution
		 is Gaussian.  If you changed the distribution and want to revert
       to a Gaussian, your input string should contain "GAUSS".
		 Your string should take the form of a C++ function, and can contain
       "[0]" or "[1]".  The parameter [0] will be the original
		 value of the variable the device smears, (plus a bias, if you set one),
       and the parameter [1] will be given by the device's
		 parametrization.  Exempli gratia, for the default Gaussian, [0] is
       the mean and [1] is the standard deviation.
		 
		 A custom distribution is automatically normalized to 1.
       Custom PDF's may cause smearing to run extremely slowly.
		 */
		void SetDistribution(TString& s);
		
		/**
		 Set the probability distribution function on which smeared values
       are generated to a general C++ function.  
		 Your function should be of the form f(double *x, double *par).
       The first argument is the argument of your distribution,
		 the second can contain up to two parameters par[0] and par[1]
       where par[0] will be the original (unsmeared) value of the 
		 variable to be smeared and par[1] is given by the device
       parametrization.  See Root TF1 documentation for more information.
		 */
		void SetDistribution(Distributor::Function f, int npars);
		
		/** 
		 Set the range of validity of the custom probability distribution.
       It is recomended that you set a distribution range whenever
		 you use a custom probability distribution.
		 */ 
		void SetDistributionRange(double min, double max);
		
		/**
		 Set the range of validity of the custom probability distribution
       to depend on the original (unsmeared) value x.  The range
		 in which values are generated will then be [x-minus,x+plus].
       This is useful if you want to generate new values in a flat
		 distribution in a neighborhood of the original value, which can
       be an important troubleshooting tool.  Setting a range
		 will automatically activate this feature, but not the use of the
       custom probability distribution.
		 */
		void SetMoveableRange(double plus, double minus);
		
		/**
		 Turn the dependence of the range of a custom probability
       distribution on the original (unsmeared) value on or off.  
		 See SetMoveableRange documentation.
		 */
		void SetMoveable(bool b);
      
		/**
		 Set a bias to the mean of your probability distribution.
       The effective mean will then become the original value of the
		 smeared variable plus the bias you set.
		 */
		void SetDistributionBias(double bias);
		
		KinType GetTargetVariable();
		
		const TString& GetDeviceType();
		
		/**
		 Provides a string with the smearing parametrization to be used.
       This will 
		 parametrize the standard deviation of the smeared variable.
       Arguments must
		 be of the form "f(x)" using C++ operators for 1d parametrizations, and 
		 "f(x,y)" for 2d parametrizations.  If you do not have the proper
       number of arguments
		 it will throw an exception, or crash.  
		 
		 For example, if you have set the parametrization to take E as an
       argument, and
		 smear the variable E and do Dev.SetParametrization("0.01*x"), the
       standard deviation
		 of the Gaussian distribution on which new E values are generated
       will go as "0.01*E".
		 
		 As another example, if you have set the parametrization to take two
       arguments, p and theta
		 and smear the variable theta, for instacne using
       Dev.SetInputKinematics(kP,kTheta) then
		 Dev.SetSmearedKinematics(kTheta) and then you do
       Dev.SetParametrization("0.01*x+0.01*sin(y)"),
		 the standard deviation of the Gaussian distribution on which new
       theta values are generated will
		 be "0.01*p+0.01*sin(theta)".  The first variable you set as an
       input argument will be called
		 x, and the second will be called y.
		 
		 Alternatively, you can write your functions to be parsed and
       automatically set input variables.
		 For example, instead of the above you can do "0.01*P+0.01*sin(theta)".
       Variable names to be used
		 are "E", "P", "theta", "phi", "pZ" and "pT".  Of course you will
       still need to use SetSmearedKinematics.
		 
		 See Root TF1 documentation for more information.
		 */
		void SetParametrization(TString par, bool parse=true);
		
		void SetParametrization(void *f, int npars);
		
		/**
		 Set the device to have exactly the same acceptance as the device
       provided.
		 */
		void GetAcceptanceFromDevice(Device dev);
      
		/** 
		 Set the range in which the parametrization is valid.
       By default, it is valid in a range
		 [0, 1e6] so hopefully you won't ever have to change this.
		 */
		void SetParRange(double min, double max);
		
		/**
		 Set the range in which the parametrization is valid for the
       first argument and the second
		 argument seperately.  Defaults are [0,1e6] for both.
		 */
		void SetParRange(double min1, double min2, double max1, double max2);
		
		/**
		 Set up to 6 parameters for the TF1, TF2 function class.
       See Root TF1,TF2 documentation.
		 For the most part, you shouldn't ever need to use this,
       it is just to provide a convenience
		 in case you have a particularly complicated parametrization.
		 */
		void SetParams(double xi0, double xi1=0., double xi2=0.,
                     double xi3=0., double xi4=0., double xi5=0.);
		
		void SetParams(int n, double xi);
		
		/**
		 This returns the Root TRandom3 random number generator instance
       being used by the device.
		 */
		TRandom3 *GetRandomGenerator();
		
		/**
		 Set the random number generator seed for the Root TRandom3 random
       number generator instance
		 being used by the device.  See Root TRandom3 documentation for
       more information.
		 */
		void SetRanSeed(int n);
		
		//note: every inherited class should contain its own version of
      //this function to work properly with Detector
		virtual Device* Clone();
		
		virtual double EvaluateRes(double x);
		
		virtual double EvaluateRes(double x, double y);
      
		/**
		 Device level particle kinematic smearing.  This function will take
       a particle of class Particle
		 and a special, stripped-down particle of class ParticleS and smear
       kinematics in the ParticleS
		 in accordance with the manner defined by the user.
       Particles which do not lie within the acceptance
		 of the device, or are not of the genre (id est, photons and leptons
       or hadrons) associated with the 
		 device will not be smeared, and the ParticleS provided will be
       unchanged.  Also, only the variable 
		 of the type which the device smeares (E, p, theta, phi, pz, pt)
       will be affected, others will not be changed.
		 In particular, this function does NOT copy unaltered values from
       the Particle to the ParticleS.  So,
		 if your device is set to smear E, this function will not in any
       way change the value of p, theta or
		 phi in the ParticleS.
		 
		 Smearing works in the following way.  If we smear a variable
       Particle.X=x with parametrization f(xi), 
		 then z[x,f(xi)] will be written to ParticleS.X where z is a
       randomly generated number from a Gaussian
		 distribution of which x is the mean and f(xi) is the standard
       deviation, where xi is a set of up to two
		 variables from the original Particle.
		 
		 If the smearing process results in negative values of positive
       definite quantities E and p, they will
		 instead be set to -999.
		 */
		virtual void DevSmear(const Particle& prt, ParticleS& prtOut);
      
      // TODO implement data-hiding
//   protected:
      
		Distributor Distribution;
		
		//the input/output kinematics 
		KinType InKin1;
      KinType InKin2;
      KinType OutKin;
		
		//dimensionality of parametrization.  Must be either 1 or 2
		int ParDim;
		//minimum and maximum of parametrizations
		double ParMin1;
      double ParMax1;
      double ParMin2;
      double ParMax2;
      // TODO ROOT allows up to 3 variables (TF3) - look into support.
		//the TF1 class for parametrizing the uncertainty
		TF1 Param1;
		//likewise for TF2
		TF2 Param2;
      
		TString name;
		
		/**
		 This is the instance of the acceptance class associated with
       the device.  Use this object to set the acceptance
		 of the device.
		 */
      Acceptance Accept;
      
      // TODO look into using functor for generic user-defined function.
      ROOT::Math::ParamFunctor* mFunctor;
//      TF1*
      
      /**
       Pure virtual base class for a generic smearing function.
       The derived class must be copy-constructible to be compatible
       with Device.
       */
      struct Function : public TObject {
         
         virtual ~Function() { }
         
         /**
          Implement the function operator to conduct the desired smearing.
          The first argument is an input particle with Monte Carlo values.
          The second argument is the smeared particle whose properties are
          to be set.
          It should return true in the event of success or false in the
          event of an error.
          */
         virtual bool operator()(const Particle&, ParticleS&) = 0;
         
         ClassDef(Function, 0)
      };
      
      Function* mFunction;
      
      void SetFunction(const Function& f) {
         mFunction = dynamic_cast<Function*>(f.Clone());
      }
      
   private:
      
		ClassDef(Device, 1)
	}; // class Device
   
   /**
    Example custom function.
    Copies values straight from a Particle to a ParticleS.
    */
   struct Copier : public Device::Function {
      virtual bool operator()(const Particle& p, ParticleS& s) {
         s.p = p.GetP();
         s.E = p.GetE();
         s.theta = p.GetTheta();
         return true;
      }
      
      ClassDef(Copier, 0)
   };
   
   /**
    Example custom function.
    Applies Gaussian smearing to the Particle energy.
    */
   struct Gauss : public Device::Function {
      virtual bool operator()(const Particle& p, ParticleS& s) {
         if(not gRandom) return false;
         s.E = gRandom->Gaus(p.GetE(), 0.35 * ::sqrt(p.GetE()));
         return true;
      }
      
      ClassDef(Gauss, 0)
   };
   
} // namespace Smear

#endif