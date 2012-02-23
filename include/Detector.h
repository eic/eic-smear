#ifndef _ERHIC_BUILDTREE_DETECTOR_
#define _ERHIC_BUILDTREE_DETECTOR_

#include <vector>

#include <TString.h>

#include "VirtualParticle.h"
#include "SmearEvent.h"

// TODO mixture of owned and non-ownded Devices makes memory management hard.
// TODO Aside from keeping separate "owned" & "non-owned" lists, should be
// TODO easier to keep a single list and always add a clone (i.e. own all).
// TODO I think the ParticleID class should inherit from Device and be treated
// TODO the same as any other device - they essentially just "smear" id.
// TODO Add consts where possible.
// TODO Make Detector inherit from TObject (other classes, like Device, will
// TODO also need this. This will allow the Detector itself to be written to
// TODO the smeared file, which will be useful for documentation purposes.

namespace Smear {
	
   class Device;
   class ParticleID;
   
	/** 
	 The detector structure.
    A detector can be equiped with devices, and contains its own
	 ParticleID instance.  It contains a detector level smearing function which
    applies smearing of all its devices to the provided ParticleS and
    generates particle ID. The detector can also generate event-wise smeared
    kinematics if provided with smeared particles from that event.
	 
	 Detector level particle and event smearing is used by the SmearTree
    function.
    
    \todo Do something about the mix of object ownership with added Devices -
    they may or may not need to be deleted, and there is no tracking of which
    is which. Either make all or none owned by the Detector.
	 */
	class Detector {
		
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
		
		/**
		 Delete all devices referenced by the detector.
		 */
		void DeleteAllDevices();
		
      /**
       Clears device list but does not free memory.
       */
		void RemoveAllDevices();
		
		/**
		 Delete all particle identifiers referenced by from the detector.
		 */
		void DeleteAllIdentifiers();
		
      /**
       Clears particle identifier list but does not free memory.
       */
		void RemoveAllIdentifiers();
		
      /**
       Replace this with reference-only version
       */
		void AddDevice(Device*);
      
		/**
		 Add a device to the detector.
       The detector will use all its devices when applying smearing. 
		 Detectors are labeled by integers in the order in which you added them
       minus 1.
		 
		 Note that when you add a device, the detector will create a clone of
       the device you added, and then use it instead of your original device.
       This is done so that you don't have to worry about including devices in
       a program that terminates before you are done with your detector.
       To remove the clone used by the Detector, use RemoveDevice.
		 
		 To add a device you created, and not a clone, use AddDevice(Device*).
		 */
		void AddDevice(Device &dev);
		
      /**
       Replace this with reference-only version
       */
		void AddDevice(ParticleID*);
		
		/**
		 Add a particle identifier to the detector.  The detector will use all
       of its particle identifiers when applying smearing.
       Particle identifiers are cloned, as with devices, see AddDevice(Device) 
		 documentation for more information. 
		 */
		void AddDevice(ParticleID&);
		
		/**
		 Set the beam lepton.
       This is needed for event kinematic smearing.
       Must use PDG particle codes.
		 The default is electrons (11).
		 */
		void SetPDGLeptonCode(int);
		
		/**
		 Set the acceptance in the variable associated with type of all equiped
       devices to [min,max].
		 Useful if you want to quickly set your detector to pick up everything.
		 */
		void SetMasterAccept(KinType, double min, double max);
		
		/**
		 Set the parametrization of all equiped devices. 
       Very useful for troubleshooting.
		 */
		void SetMasterParametrization(TString);
		
		/**
		 Set the TRandom3 random number generator seed of all equiped devices.
       Very useful for troubleshooting.
		 */
		void SetMasterRanSeed(int);
		
      /**
		 Set the Root TF1/TF2 parameters for the parametrization of all equiped
       devices simultaneously.
		 Useful for advanced troubleshooting.
		 */
		void SetMasterParams(int n, double xi);
		
		/**
		 Set the smearing probability distribution function for all member
       devices.
       Useful for troubleshooting.
		 */
		void SetMasterDistribution(TString);
		
		/**
		 Set the custom smearing probability distribution range for all
       member devices.
       Useful for troubleshooting.
		 */
		void SetMasterDistributionRange(double min, double max);
		
      /**
       Equivalent to AddDevice().
       */
		Detector& operator<< (Device&);
		
      /**
       Equivalent to AddDevice().
       */
		Detector& operator<< (Device*);
		
      /**
       Equivalent to AddDevice().
       */
		Detector& operator<< (ParticleID&);
		
      /**
       Equivalent to AddDevice().
       */
		Detector& operator<< (ParticleID*);
		
      /**
       Equivalent to SetEventKinematicsCalculator()
       */
		Detector& operator << (TString EKCalc);
		
      /**
       Equivalent to GetDevice()
       */
		Device* operator[] (int);
		
		void SetMissingEnergyTolerance(double);
		
		void TolerateBadEvents(double);
		
		void SetSupressEventWarnings(bool);
		
		/** 
		 Set the method for calculating event kinematics if FillEventKinematics
       is used. String must contain "NM" for null momentum approximation
       (using scattered lepton), "JB" for Jacquet Blondel method,
		 or "DA" for double angle method.
       Strings not containing one of these turns the method off.
		 */
		void SetEventKinematicsCalculator(TString);
		
		/**
		 Remove device number n from the detector.
       Devices are labeled in the order in which they are added
		 minus 1.
		 */
		void RemoveDevice(int);
		
		/**
		 Remove particle identifier number n from the detector.
       Identifiers are labeled in the order in 
		 which they are added minus 1.
		 */
		void RemoveIdentifier(int);
		
		/**
		 Return a pointer to device number n from the detector.
       Devices are labeled in the order in which
		 they are added minus 1.
		 */
		Device* GetDevice(int);
		
		/** 
		 Return a pointer to particle identifier number n from the detector.
       Identifiers are labeled in the order in which they are
		 added minus 1.
		 */
		ParticleID* GetIdentifier(int);
		
		/**
		 This returns the device used to smear variable type kin of the provided
       particle when using detector level smearing.
		 */
		Device* GetDeviceUsedFor(Particle, KinType);
		
		/**
		 Returns the number of the device used to smear variable of type kin of
       the provided particle when using detector
		 level smearing.
		 */
		int WhichDeviceUsedFor(Particle, KinType);
		
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
		void FillEventKinematics(const EventBase*, EventS*);
		
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
		ParticleS* DetSmear(const Particle&);
		
      // TODO implement data-hiding
//   protected:
      
		std::vector<Device*> Devices;
		std::vector<ParticleID*> Identifiers;
		
		bool useNM;
      bool useJB;
      bool useDA;
      
		EventKinematicsComputer EventKinComp;
		
      UInt_t GetNDevices() const { return Devices.size(); }
//		int NDevices; 
		int NIdentifiers;
		
//		bool bPID; 
		
   private:
      
		ClassDef(Detector, 1 )
	}; 
   
   inline Detector& Detector::operator<< (Device &dev) {
      AddDevice(dev);
      return *this;
   }
   
   inline Detector& Detector::operator<< (Device *dev) {
      AddDevice(dev);
      return *this;
   }
   
   inline Detector& Detector::operator<< (ParticleID &ident) {
      AddDevice(ident);
      return *this;
   }
   
   inline Detector& Detector::operator<< (ParticleID *ident) {
      AddDevice(ident);
      return *this;
   }
   
   inline Detector& Detector::operator << (TString EKCalc) {
      SetEventKinematicsCalculator(EKCalc);
      return *this;
   }
   
   inline Device* Detector::operator[] (int n) {
      return GetDevice(n);
   }
}

#endif
