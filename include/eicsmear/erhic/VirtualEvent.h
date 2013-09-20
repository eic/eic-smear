/**
 \file
 Declaration of class erhic::VirtualEvent.
 
 \author    Thomas Burton
 \date      2011-08-19
 \copyright 2011 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_ERHIC_VIRTUALEVENT_H_
#define INCLUDE_EICSMEAR_ERHIC_VIRTUALEVENT_H_

#include <TObject.h>

#include <vector>

namespace erhic {

class VirtualParticle;

/**
 \brief   Abstract base class for a physics event.
 \details Defines an "event" simply as a collection of tracks.
 */
class VirtualEvent : public TObject {
 public:
  /**
   Destructor
   */
  virtual ~VirtualEvent() { }

  /**
   Returns the nth track from the event.
   
   Indices run from 0 to (n-1).
   */
  virtual const VirtualParticle* GetTrack(UInt_t /* number */) const = 0;

  /**
   \overload
   */
  virtual VirtualParticle* GetTrack(UInt_t /*number*/) = 0;

  /**
   Returns the number of tracks in the event.
   */
  virtual UInt_t GetNTracks() const = 0;

  /**
   typedef for a track pointer collection.
   */
  typedef std::vector<const erhic::VirtualParticle*> ParticlePtrList;

  /**
   Populate a track list with the hadronic final-state.
   */
  virtual void HadronicFinalState(ParticlePtrList&) const { }

  ClassDef(erhic::VirtualEvent, 1)
};

}  // namespace erhic

#endif  // INCLUDE_EICSMEAR_ERHIC_VIRTUALEVENT_H_
