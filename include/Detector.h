#ifndef _ERHIC_BUILDTREE_DETECTOR_
#define _ERHIC_BUILDTREE_DETECTOR_

#include <vector>

#include <TObject.h>
#include <TString.h>

#include "SmearEvent.h" // For EventKinematicsComputer

namespace erhic {

   class EventMC;
   class ParticleMC;

} // namespace erhic

namespace Smear {
	
   class Event;
   class Smearer;

	/**
	 The detector structure.
    A detector posseses a collection of Smearer objects, each smearing
    some variable(s).
    It contains a detector-level smearing function which
    applies smearing of all its devices to the provided ParticleMCS.
    The detector can also generate event-wise smeared
    kinematics if provided with smeared particles from that event.
    
    \todo Add consts where possible.
    \todo Implement data hiding
	 */
	class Detector : public TObject {
		
   public:
      
		/** 
		 Default contructor.
       By default particle ID is off, and the detector thinks the beam
		 lepton is an electron.
       \todo Add optional arguments to construtor to set PID, lepton species,
         event kinematics calculations
         Detector(pid, lepton, kinematics)
		 */
		Detector();
      
      /** Copy constructor */
      Detector(const Detector&);
		
      /** Destructor */
      virtual ~Detector();
      
		/**
		 Delete all devices referenced by the detector.
		 */
		void DeleteAllDevices();

		/**
		 Adds a copy of the device to this detector.
       The detector will use all its devices when applying smearing. 
		 Detectors are labeled by integers in the order in which you added them
       minus 1.
		 */
		void AddDevice(Smearer& dev);

		/**
		 Set the beam lepton.
       This is needed for event kinematic smearing.
       Must use PDG particle codes.
		 The default is electrons (11).
       \remark Simply accessing the EventKinematicsComputer may be more
       straightforward
		 */
		void SetPDGLeptonCode(int);

      /**
       Equivalent to SetEventKinematicsCalculator()
       \remark I think this syntax is confusing (it's nothing to do with streams)
       */
		Detector& operator<<(TString EKCalc);

      /**
       \remark Simply accessing the EventKinematicsComputer may be more
       straightforward
      */
		void SetMissingEnergyTolerance(double);
		
      /**
       \remark Simply accessing the EventKinematicsComputer may be more
       straightforward
       */
		void TolerateBadEvents(double);
		
      /**
       \remark Simply accessing the EventKinematicsComputer may be more
       straightforward
       */
		void SetSupressEventWarnings(bool);
		
		/** 
		 Set the method for calculating event kinematics if FillEventKinematics
       is used. String must contain "NM" for null momentum approximation
       (using scattered lepton), "JB" for Jacquet Blondel method,
		 or "DA" for double angle method.
       Strings not containing one of these turns the method off.
       \remark Simply accessing the EventKinematicsComputer may be more
       straightforward
		 */
		void SetEventKinematicsCalculator(TString);
		
		/**
		 Remove device number n from the detector.
       Devices are labeled in the order in which they are added
		 minus 1.
       \remark A Device* argument may make more sense
		 */
		void RemoveDevice(int);

		/**
		 Return a pointer to device number n from the detector.
       Devices are labeled in the order in which they are added minus 1.
       Do not delete the returned pointer.
		 */
		Smearer* GetDevice(int);

		/**
		 Calculate event-wise smeared kinematics for an event which has already
       had its particles smeared and stored in eventS.
       Newly calculated values are stored in eventS.
		 
		 This uses the null momentum (m->0) approximation for the scattered
       lepton, so you should be careful about using low energy muons, and if
       for some strange reason you want to use taus, this probably won't work
       too well.
		 
		 Also, the smeared lepton momentum (as opposed to energy) is used in
       the assumption that its smearing is less severe.
		 */
		void FillEventKinematics(const erhic::EventMC*, Event*);
		
		/**
		 Detector level particle smearing.  This is intended to be the primary
       method for smearing particles. The function will output a pointer to a
       particleS, which is a version of the input Particle which has been
       smeared by all the detectors devices, and identified with the
       detector's particle ID (if on). Note that the ParticleS class contains
       only E,p,theta,phi,pz,pt and id.  This smearing will only be applied
		 particles which are stable (i.e. in the Pythia sense).
		 
		 If the input particle is unstable, an initial state or intermediate
       state particle, the detector smearing will
		 return a null ParticleS pointer.  
		 */
		ParticleMCS* DetSmear(const erhic::ParticleMC&);
		
		std::vector<Smearer*> Devices;
		
		bool useNM;
      bool useJB;
      bool useDA;
      
		EventKinematicsComputer EventKinComp;
		
      UInt_t GetNDevices() const;
      std::vector<Smear::Smearer*> CopyDevices() const;

   private:

      Detector& operator=(const Detector&);

		ClassDef(Detector, 1 )
	};

   inline Detector& Detector::operator << (TString EKCalc) {
      SetEventKinematicsCalculator(EKCalc);
      return *this;
   }

   inline UInt_t Detector::GetNDevices() const {
      return Devices.size();
   }

   inline Detector& Detector::operator=(const Detector&) {
      return *this;
   }
}

#endif
