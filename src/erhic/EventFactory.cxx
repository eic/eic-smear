/**
   \file
   Implementation of class erhic::EventFactory.
 
   \author    Thomas Burton
   \date      2011-10-31
   \copyright 2011 Brookhaven National Lab
*/

#include "eicsmear/erhic/EventFactory.h"

#include <memory>
#include <stdexcept>
#include <string>

#include <TClass.h>
#include <TProcessID.h>

#include "eicsmear/erhic/BeamParticles.h"
#include "eicsmear/erhic/EventPythia.h"
#include "eicsmear/erhic/EventHepMC.h"
#include "eicsmear/erhic/EventMilou.h"
#include "eicsmear/erhic/EventDjangoh.h"
#include "eicsmear/erhic/EventDpmjet.h"
#include "eicsmear/erhic/EventRapgap.h"
#include "eicsmear/erhic/EventPepsi.h"
#include "eicsmear/erhic/EventGmcTrans.h"
#include "eicsmear/erhic/EventSimple.h"
#include "eicsmear/erhic/EventSartre.h"
#include "eicsmear/functions.h"  // For getFirstNonBlank()
#include "eicsmear/erhic/Kinematics.h"
#include "eicsmear/erhic/ParticleIdentifier.h"
#include "eicsmear/erhic/ParticleMC.h"

#include <TVector3.h>
#include <TParticlePDG.h>
#include <TLorentzVector.h>
#include <TDatabasePDG.h>

#include<map>

using std::cout;
using std::cerr;
using std::endl;
using std::map;

namespace erhic {

  template<typename T>
  bool EventFromAsciiFactory<T>::AtEndOfEvent() const {
    return !mLine.empty()
      && mLine.find("finished") != std::string::npos;
  }

  // Use this struct to automatically reset TProcessID object count.
  struct TProcessIdObjectCount {
    // Initialse object with current TProcessID object count.
    TProcessIdObjectCount() {
      count = TProcessID::GetObjectCount();
    }
    // Restore object count to the value at initialisation.
    // See example in $ROOTSYS/test/Event.cxx
    // To save space in the table keeping track of all referenced objects
    // we assume that our events do not address each other.
    ~TProcessIdObjectCount() {
      TProcessID::SetObjectCount(count);
    }
    int count;
  };

  template<typename T>
  T* EventFromAsciiFactory<T>::Create() {
    // Save current object count. Will reset it when this function returns.
    TProcessIdObjectCount objectCount;
    mEvent.reset(new T);
    // We use this flag to check input doesn't end mid-event.
    // Initialised finished flag to "success" in case of no input.
    int finished(0);
    std::string error;
    // Read line-by-line until the stream is not good or we break out.
    while (std::getline(*mInput, mLine).good()) {
      // Reached end-of-event marker
      if (AtEndOfEvent()) {
	// If we built a good event (i.e. no errors reading strings)
	// perform end-of-event calculations.
	if (!error.empty()) {
	  throw std::runtime_error(error);
	} else {
	  finished = FinishEvent();  // 0 upon success
	  break;
	}  // if
      } else if ('0' == getFirstNonBlank(mLine)) {
	// '0' indicates the event header line
	// An event started, set finished flag to "unfinished".
	finished = -1;
	// Parse string and check for validity.
	if (!mEvent->Parse(mLine)) {
	  // Set error message based on bad event input.
	  // Don't break out of the loop yet, so that we continue reading
	  // lines until the end-of-event marker. That way we stay
	  // "aligned" with the input data ready for the next event.
	  error = "Bad event input: " + mLine;
	}  // if
      } else if ('=' != getFirstNonBlank(mLine)) {
	// Anything remaining other than a line of '=' is a particle line
	// Don't raise an exception for a failed track, as the event in
	// general may be OK. AddParticle will print a message though.
	if (!AddParticle()) {
	  error = "Bad particle input in event";
	}  // if
      }  // if
    }  // if
    // Check for events that started but did not finish.
    // We should have ended with an end-of-event marker, meaning the finished
    // flag is set to zero. If not, the finished flag will be non-zero.
    if (finished != 0) {
      mEvent.reset(NULL);
      throw std::runtime_error("Ended mid-event");
    }  // if
    // Return a NULL event *without* throwing an exception to indicate
    // end-of-file. We shouldn't have hit eof if we have read a good event,
    // as we won't have yet tried to read past the end (mLine will still be
    // at the end-of-file marker line).
    if (mInput->eof()) {
      mEvent.reset(NULL);
    }  // if
    return mEvent.release();
  }

  template<typename T>
  Int_t EventFromAsciiFactory<T>::FinishEvent() {
    std::unique_ptr<DisKinematics> nm(
				      LeptonKinematicsComputer(*mEvent).Calculate());
    std::unique_ptr<DisKinematics> jb(
				      JacquetBlondelComputer(*mEvent).Calculate());
    std::unique_ptr<DisKinematics> da(
				      DoubleAngleComputer(*mEvent).Calculate());
    if (nm.get()) {
      mEvent->SetLeptonKinematics(*nm);
    }  // if
    for (unsigned n(0); n < mEvent->GetNTracks(); ++n) {
      mEvent->GetTrack(n)->ComputeEventDependentQuantities(*mEvent);
    }  // for
    if (jb.get()) {
      mEvent->SetJacquetBlondelKinematics(*jb);
    }  // if
    if (da.get()) {
      mEvent->SetDoubleAngleKinematics(*da);
    }  // if
    // We also have to set the remaining variables not taken care of
    // by the general DIS event kinematic computations.
    // Find the beams, exchange boson, scattered lepton.
    BeamParticles beams;
    if (!ParticleIdentifier::IdentifyBeams(*mEvent, beams)) {
      std::cerr <<
	"EventFromAsciiFactory::FinishEvent(): failed to find beams"
		<< std::endl;
      return 0;
    }  // if
    const TLorentzVector h = beams.BeamHadron();
    TLorentzVector l = beams.BeamLepton();
    TLorentzVector s = beams.ScatteredLepton();
    TVector3 boost = -h.BoostVector();
    l.Boost(boost);
    s.Boost(boost);
    mEvent->SetELeptonInNuclearFrame(l.E());
    mEvent->SetEScatteredInNuclearFrame(s.E());
    return 0;
  }

  template<typename T>
  bool EventFromAsciiFactory<T>::AddParticle() {
    //return true;
    try {
      if (mEvent.get()) {
	ParticleMC particle(mLine, mEvent->RequiresEaParticleFields());  // Throws if the string is bad
	particle.SetEvent(mEvent.get());
	mEvent->AddLast(&particle);
	//ParticleMCeA *particle = new ParticleMCeA(mLine);  // Throws if the string is bad
	//particle->SetEvent(mEvent.get());
	//mEvent->AddLast(particle);
	//delete particle;
      }  // if
      return true;
    }  // try
    catch(std::exception& error) {
      std::cerr << "Exception building particle: " << error.what() << std::endl;
      return false;
    }
  }

  template<typename T>
  std::string EventFromAsciiFactory<T>::EventName() const {
    return T::Class()->GetName();
  }

  template<typename T>
  void EventFromAsciiFactory<T>::FindFirstEvent()  {
    for (int i=0; i<5; i++)
      {
	std::getline(*mInput,mLine);
      }
  }

    
}  // namespace erhic

namespace {

  // Need this to generate the CINT code for each version
  erhic::EventFromAsciiFactory<erhic::EventDjangoh> ed;
  erhic::EventFromAsciiFactory<erhic::EventDpmjet> ej;
  erhic::EventFromAsciiFactory<erhic::EventPepsi> ee;
  erhic::EventFromAsciiFactory<erhic::EventMilou> em;
  erhic::EventFromAsciiFactory<erhic::EventHepMC> eh;
  erhic::EventFromAsciiFactory<erhic::EventRapgap> er;
  erhic::EventFromAsciiFactory<erhic::EventPythia> ep;
  erhic::EventFromAsciiFactory<erhic::EventGmcTrans> eg;
  erhic::EventFromAsciiFactory<erhic::EventSimple> es;
  erhic::EventFromAsciiFactory<erhic::EventSartre> esa;
  erhic::EventFromAsciiFactory<erhic::EventBeagle> eb;

}  // namespace
