#ifndef _ERHIC_BUILDTREE_ACCEPTANCE_
#define _ERHIC_BUILDTREE_ACCEPTANCE_

#include <vector>

#include <TF1.h>
#include <TF2.h>
#include <TString.h>

#include "Particle.h"
#include "Smear.h"

/*
 2011/08/16
 Forced RemoveZone() argument to be with number of zones, which fixed
 erroneous deletion of Zones previously.
 */

namespace Smear {
   
	/**
	 Determines acceptance for a Device.
    Every Device contains an instance of this structure called Accept.
    Comprises one or more (potentially non-overlapping) acceptance zones.
	 */
   class Acceptance {
		
   public:
      
		/**
		 Default constructor.
       By default, the device has 4pi coverage,
       and accepts particles with energy and momenta up
		 to 1e12 GeV.  Also, the genre is neutral, meaning the acceptance
       function doesn't care what type of particle it 
		 sees.
		 */
		Acceptance(int genre = kE);
      
		/**
		 Add a new default acceptance zone with 4pi coverage.
       Particles will be accepted if they fall within any 
		 acceptance zone.
		 */
      // Replace with default arguments for next method
      //		void AddZone();
		
		/**
		 Add a new zone with use specified coverage.
       Particles will be accepted if they fall within any acceptance
		 zone.
		 */
		void AddZone(double thetamin = 0., double thetamax = pi,
                   double phimin = 0., double phimax = 2. * pi,
                   double Emin = 0., double Emax = 1.e12,
                   double pmin = 0., double pmax = 1.e12);
		
      /**
       Returns the number of acceptance zones.
       */
      UInt_t GetNZones() const;
      
		/**
		 Remove an acceptance zone.
       The nth zone you added is Zone n, the original zone is Zone 0.
		 */
		void RemoveZone(unsigned n);
		
		/**
		 Set the acceptance in theta. 
		 */
		void SetTheta(double min, double max, int n=0);		
		/**
		 Set the acceptance in phi.  
		 */
		void SetPhi(double min, double max, int n=0);      
		/** 
		 Set the acceptance in E.  
		 */
		void SetE(double min, double max, int n=0);		
		/**
		 Set the acceptance in p.
		 */
		void SetP(double min, double max, int n=0);		
		/**
		 Set acceptance in the kinematic variable assiciated with type.
       Only valid for
		 kE, kP, kTheta, kPhi.
		 */
		void Set(KinType type, double min, double max, int n=0);		
		/**
		 Set the acceptance in pT.
       Note that this uses the AddCustomAcceptance method.
       Therefore, if for some reason
		 you want to call this more than once (with the purpose of changing
       rather than adding a window) you should
		 use ClearCustomAcceptance.  
		 */
		void SetPt(double min, double max, int n=0);		
		/**
		 Set the acceptance int pz.
       Note that this uses the AddCustomAcceptance method.
       Therefore, if for some reason
		 you want to call this more than once (with the purpose of changing
       rather than adding a window) you should 
		 use ClearCustomAcceptance.
		 */
		void SetPz(double min, double max, int n=0);
      
		void SetGenre(int n);
		
//		void SetGenre(TString genre);
      
		/**
		 Allows you to set the acceptance in some function of E,p,theta,phi
       for zone n.  For example, if you want to set the
		 acceptance in pT to [0.,100.] you can do
		 
		 Accept.AddCustomAcceptance("P*sin(theta)",0.,100.);
		 
		 For this you can use up to 2 input variables (must be E,p,theta,phi)
       using the same syntax as for the 
		 Device::SetParametrization method.
       Calling AddCustomAcceptance more
       than once will require that the particle
		 fall within all of the custom acceptances you added.
		 */
		void AddCustomAcceptance(TString s, double min, double max, int n=0);		
      
		/**
		 Clear all custom acceptance ranges from zone n.
       Acceptance will again be determined soley by the E,p,theta,phi
		 acceptance windows.
		 */
		void ClearCustomAcceptance(int n=0);	
      
		/**
		 Add a particle type to the list of particles to be smeared.
       If you never add anything, the device will
		 smear all stable particles within its acceptance.
       Must use PDG particle codes.
		 */
      void AddParticle(int);
      
		/**
		 Clear the list of particles to be smeared by this device.
       The device will then smear all stable particles within
		 its acceptance.
		 */
		void ClearParticles();
      
		/**
		 Remove particle from the list of particles smeared by this device.
       Must use PDG particle code.
		 */
		void RemoveParticle(int);
      
      //======================================================================
      // TODO Does this need to be public?
      //======================================================================
		struct CustomCut {
			
			CustomCut();
			
			TF1 F1;
			TF2 F2;
			int dim;
			KinType Kin1;
         KinType Kin2;
			double Min;
         double Max;
		};
		
      //======================================================================
      // TODO Does this need to be public?
		//acceptance parameters
      //======================================================================
		struct Zone {
         
			Zone(double = 0., double = pi, double = 0., double = 2. * pi,
              double = 0., double = 1.e12, double = 0., double = 1.e12);
         
         Bool_t Contains(Particle&) const;
         
         Bool_t IsCustomAccepted(CustomCut, const Particle&) const;
			
			double thetaMin;
         double thetaMax;
			double phiMin;
         double phiMax;
			double EMin;
         double EMax;
			double PMin;
         double PMax;
			
			std::vector<CustomCut> CustomCuts;
		};
      
      bool IsCustomAccepted(CustomCut& C, const Particle& prt) const;
      
		/**
		 This function determines if the particle provided lies within
       the acceptance of the
		 detector.  Default acceptance is a full 4*pi solid angle,
       with E and p in (0.2,1e12) GeV.
		 This function automatically fixes polar and azimuthal angles
       which are not within their proper range.
		 */
		bool Is(const Particle& prt);
      
//   protected: TODO implement data hiding
		std::vector<Zone> mZones;
		std::vector<int> Particles;
		
		int Genre;
      
   private:
      
//#ifdef __CINT__
		ClassDef(Acceptance, 1)
//#endif
   };
   
   inline UInt_t Acceptance::GetNZones() const { return mZones.size(); }
   
} // namespace Smear

#endif
