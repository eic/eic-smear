/**
 \file
 Declaration of class Smear::Detector.
 
 \author    Michael Savastio
 \date      2011-08-19
 \copyright 2011 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_SMEAR_DETECTOR_H_
#define INCLUDE_EICSMEAR_SMEAR_DETECTOR_H_

#include <list>
#include <vector>

#include <TObject.h>
#include <TString.h>

namespace erhic {

class EventDis;
class VirtualParticle;

}  // namespace erhic

namespace Smear {

class Event;
class ParticleMCS;
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
   */
  Detector();

  /**
   Copy constructor.
   */
  Detector(const Detector&);

  /**
   Assignment operator.
   */
  Detector& operator=(const Detector&);

  /**
   Destructor.
   */
  virtual ~Detector();

  /**
   Adds a copy of the smearing device to this detector.
   The detector will use all its devices when applying smearing. 
   */
  void AddDevice(Smearer& device);

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
   Delete all devices in the detector.
   */
  void DeleteAllDevices();

  /**
   Return a pointer to device number n from the detector.
   Devices are labeled in the order in which they are added minus 1.
   Do not delete the returned pointer.
   */
  Smearer* GetDevice(int index);

  /** Returns the number of devices in the detector */
  UInt_t GetNDevices() const;

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
  void FillEventKinematics(Event* event);

  /**
   Detector level particle smearing.
   Returns a pointer to a new smeared particle, which is a version of
   the input Particle that has been smeared by all the detector's devices.
   Smearing is only be applied to final-state particles.
   If the input particle is unstable or an initial- or intermediate-
   state particle, returns a null ParticleS pointer.  
   */
  ParticleMCS* Smear(const erhic::VirtualParticle&) const;

  /**
   Print information about all smearers to standard output.
   */
  virtual void Print(Option_t* = "") const;

  /**
   Returns the list of devices in this detector that accept a particle.
   */
  std::list<Smear::Smearer*> Accept(const erhic::VirtualParticle&) const;

 protected:
  /**
   Returns pointers to new copies of all devices.
   */
  std::vector<Smear::Smearer*> CopyDevices() const;

  bool useNM;
  bool useJB;
  bool useDA;
  std::vector<Smearer*> Devices;

  ClassDef(Smear::Detector, 1)
};

inline UInt_t Detector::GetNDevices() const {
  return Devices.size();
}

}  // namespace Smear

#endif  // INCLUDE_EICSMEAR_SMEAR_DETECTOR_H_
