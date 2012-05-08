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
#include <map>
#include <vector>

#include <TF1.h>
#include <TF2.h>
#include <TFormula.h>
#include <Math/ParamFunctor.h> // For ROOT::TMath::ParamFunctor
#include <TRandom3.h>
#include <TString.h>

#include "eicsmear/smear/Acceptance.h"
#include "eicsmear/smear/Smear.h"
#include "eicsmear/smear/Smearer.h"

class TRandom3;

namespace erhic {
   class VirtualParticle;
}

namespace Smear {

   class ParticleMCS;
   
   /**
    Performs smearing of a single kinematic variable according to a simple
    expression defined via a string.
   */
   class ParamSimple : public Smearer {
   public:

      /**
       Constructor.
       The first argument is the type of kinematic variable to smear.
       The second argument is a formula giving the standard deviation of the
       resolution in the variable selected with the first argument,
          sigma(A) = f([B, C...])
       where A, B, C.. are selected from: E, P, theta, phi, pZ, pT.
       For example, for resolution in pT of 1% of pT times sin of polar angle:
       ParamSimple(kPt, "0.01 * pT * sin(theta)");
       See ROOT::TFormula for valid expressions.
       Formulae can be a function of up to four variables.
       The third argument allows selection of the types of particles that are
       smeared: electromagnetic, hadronic or all [default].
      */
      ParamSimple(KinType = kE, const TString& formula = "0", EGenre = kAll);

      /** Copy constructor */
      ParamSimple(const ParamSimple&);

      /** Destructor */
      virtual ~ParamSimple();

      /**
       Returns a dynamically allocated copy of this object.
       The argument is unused and is present for compatibility with
       ROOT::TObject::Clone().
       */
      virtual ParamSimple* Clone(const char* = "") const;

      /**
       Smear the kinematic quantity of the input ParticleMC and store the
       result in the ParticleMCS.
      */
      virtual void DevSmear(const erhic::VirtualParticle&, ParticleMCS&);

   protected:

      /**
       Process the formula, replacing patterns such as "E", "p" with variables
       recognised by ROOT::TFormula ("x", "y", "z", "t").
       Populate the vector with the KinType corresponding to pattern.
       e.g. "0.3/sqrt(E)" would yield "0.3/sqrt(x)" and store Smear::kE in the
       vector.
      */
      virtual int Parse(TString& formula, std::vector<Smear::KinType>&);

      KinType mSmeared; ///< Smeared variable
      TFormula* mFormula; ///< Expression for resolution standard deviation
      std::vector<Smear::KinType> mDimensions; ///< Variables on which smearing
                                               ///< is dependent (up to 4)
//		Distributor Distribution; ///< Random distribution

      /** Fills pattern-KinType map. */
      static void FillKinTypeTable();
      static const std::vector<TString> smPatterns;
      static std::map<TString, KinType> smKinTypes; ///< Keyed by pattern string

   private:

      ParamSimple& operator=(const ParamSimple&) { return *this; }

      ClassDef(Smear::ParamSimple, 1)
   };


	/**
	 The device class.
    This is an object that will smear a single particle-wise kinematic
	 variable according to a user provided parametrization.
    Devices can be added a detector
	 which can be used to smear data with a series of devices.
    Each device has an 
	 acceptance, which is an open subset of (E,p,theta,phi) momentum space.
    Only particles whose momentum lie in the acceptance can be smeared.
    \todo Add genre as an optional third argument (default = all)
    \todo Inherit from TObject
    \todo Get rid of ParDim. Instead select TF1/2/3 via template?
    \todo A lot of duplication of SetInOutKinematics methods. Would it make
    more sense to have different classes for 1 and 2D functions? Or something
    more sophisticated than now?
	 */
	struct Device : public Smearer {
		
		/**
		 Default constructor.  Sets a 1-dimensional parametrization f(x)=0 in
       range [0.,1.e6] by default.
       Default input kinematics E, theta, default variable to be smeared is E.
		 */
		Device(KinType = kE, TString parameterisation = "0.", int genre = kAll);

		/**
		 This determines what sorts of particles the device is able to smear.
		 
		 0: (Default)The device doesn't care what the particle is, as long as
       it is stable, 
		 (i.e. in the Pythia sense), in the final state, and not a gluon.
		 1: The device detects final state photons and leptons only.
		 2: The device detects final state hadrons only.
		 */
		void SetGenre(int);
		
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
       be an important troubleshooting tool. Setting a range
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
		
      /** \overload */
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
		
      /** \overload */
		void SetParams(int n, double xi);
		
		/**
		 This returns the Root TRandom3 random number generator instance
       being used by the device.
		 */
		TRandom3& GetRandomGenerator();
		
		/**
		 Set the random number generator seed for the Root TRandom3 random
       number generator instance
		 being used by the device.  See Root TRandom3 documentation for
       more information.
		 */
		void SetRanSeed(int n);
		
		//note: every inherited class should contain its own version of
      //this function to work properly with Detector
		virtual Device* Clone(const char* newname = "") const;
		
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
		virtual void DevSmear(const erhic::VirtualParticle& prt,
                              ParticleMCS& prtOut);
      
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

		ClassDef(Device, 1)
	}; // class Device
} // namespace Smear

#endif
