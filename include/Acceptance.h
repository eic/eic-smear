#ifndef _ERHIC_BUILDTREE_ACCEPTANCE_
#define _ERHIC_BUILDTREE_ACCEPTANCE_

#include <set>
#include <vector>

#include <TF1.h>
#include <TF2.h>
#include <TString.h>

#include "Smear.h"

namespace erhic {
   class ParticleMC;
} // namespace erhic

namespace Smear {
   
	/**
	 Determines acceptance for a Device.
    Every Device contains an instance of this structure called Accept.
    Comprises one or more (potentially non-overlapping) acceptance zones.
    
    \todo Implement data hiding
	 */
   class Acceptance {
		
   public:
      
      //========================================================================
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
      
      //========================================================================
		struct Zone {
         virtual ~Zone() { }

			Zone(double theta = 0., double = pi,
              double phi = 0., double = 2. * pi,
              double E = 0., double = 1.e12,
              double p = 0., double = 1.e12,
              double pt = 0., double = 1.e12,
              double pz = -1.e12, double = 1.e12);
         
         Bool_t Contains(erhic::ParticleMC&) const;
         
         Bool_t IsCustomAccepted(CustomCut, const erhic::ParticleMC&) const;
			
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
			
			std::vector<CustomCut> CustomCuts; //!

         // We want to be able to write Acceptance objects to a ROOT file,
         // so nested classes need a dictionary as well.
         ClassDef(Smear::Acceptance::Zone, 1)
		};
      
		/**
		 Default constructor.
       By default, the device has 4pi coverage,
       and accepts particles with energy and momenta up
		 to 1e12 GeV.  Also, the genre is neutral, meaning the acceptance
       function doesn't care what type of particle it sees.
		 */
		Acceptance(int genre = kAll);

      /** Destructor */
      virtual ~Acceptance() { }

		/**
		 Add a new zone with user-specified coverage.
       Particles will be accepted if they fall within any acceptance zone.
		 */
      void AddZone(const Zone&);

      /** Returns the number of acceptance zones. */
      UInt_t GetNZones() const;

      /** Returns the "genre" of the particle (em, hadronic, any) */
      Int_t GetGenre() const;

      /** Select the class(es) of particles to accept */
		void SetGenre(int);

		/**
		 Allows you to set the acceptance in some function of E,p,theta,phi
       for zone n.  For example, if you want to set the
		 acceptance in pT to [0.,100.] you can do
		 
		 AddCustomAcceptance("P*sin(theta)", 0., 100.);
		 
		 For this you can use up to 2 input variables (must be E,p,theta,phi)
       using the same syntax as for the 
		 Device::SetParametrization method.
       Calling AddCustomAcceptance more
       than once will require that the particle
		 fall within all of the custom acceptances you added.
		 */
		void AddCustomAcceptance(TString s, double min, double max, int n=0);		

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

		/**
		 This function determines if the particle provided lies within
       the acceptance of the
		 detector.  Default acceptance is a full 4*pi solid angle,
       with E and p in (0.,1e12) GeV.
		 This function automatically fixes polar and azimuthal angles
       which are not within their proper range.
		 */
		bool Is(const erhic::ParticleMC& prt);

   protected:

      /**
       Returns true if the ParticleMC is accepted by the acceptance defined
       in the CustomCut.
      */
      bool IsCustomAccepted(CustomCut&, const erhic::ParticleMC&) const;

		int mGenre;
		std::vector<Zone> mZones;
		std::set<int> mParticles;

		ClassDef(Acceptance, 1)
   };
   
   inline UInt_t Acceptance::GetNZones() const {
      return mZones.size();
   }
   
   inline Int_t Acceptance::GetGenre() const {
      return mGenre;
   }
} // namespace Smear

#endif
