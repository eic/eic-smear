/**
   \file
   Implementation of class erhic::EventFactoryHepMC.

   \author    Chris Pinkenburg, Kolja Kauder
   \date      2020-07-07
   \copyright 2020 Brookhaven National Lab
*/

#include "eicsmear/erhic/EventFactoryHepMC.h"

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
#include "eicsmear/erhic/EventDEMP.h"
#include "eicsmear/erhic/EventSartre.h"
#include "eicsmear/functions.h"  // For getFirstNonBlank()
#include "eicsmear/erhic/Kinematics.h"
#include "eicsmear/erhic/ParticleIdentifier.h"
#include "eicsmear/erhic/ParticleMC.h"

#include <HepMC3/ReaderAsciiHepMC2.h>
#include <HepMC3/ReaderAscii.h>
#include "HepMC3/GenVertex.h"
// The newer ReaderFactory is header-only and can be used for older versions
// This file is copied verbatim from https://gitlab.cern.ch/hepmc/HepMC3
// Copyright (C) 2014-2020 The HepMC collaboration, licensed under GPL3
#include <HepMC3/Version.h>
#if HEPMC3_VERSION_CODE < 3002004
#include <eicsmear/HepMC_3_2_4_ReaderFactory.h>
#else
#include <HepMC3/ReaderFactory.h>
#endif

#include <TVector3.h>
#include <TParticlePDG.h>
#include <TLorentzVector.h>
#include <TDatabasePDG.h>
#include <TObjArray.h>
#include <TObjString.h>

#include <map>
#include <vector>
#include <algorithm> // for *min_element, *max_element

using std::cout;
using std::cerr;
using std::endl;
using std::map;
using std::vector;

namespace erhic {

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

  template<>
  bool EventFromAsciiFactory<erhic::EventHepMC>::AtEndOfEvent() const {return false;}

  template<>
  void EventFromAsciiFactory<erhic::EventHepMC>::FindFirstEvent()  {
    // nothing to do
  }

  template<>
  bool EventFromAsciiFactory<erhic::EventHepMC>::AddParticle() {
    try {
      // open only once
      // otherwise this gets called for every event
      // and we can't make it a member since we're based on the EventFromAsciiFactory template
      static std::shared_ptr<HepMC3::Reader> adapter = HepMC3::deduce_reader(*mInput); 

      if (mEvent.get()) {
        HepMC3::GenEvent evt(HepMC3::Units::GEV,HepMC3::Units::MM);
	adapter->read_event(evt);
	if ( adapter->failed() ) return false;

	// We could take all particles "in order"
	// This isn't well defined - depends on traversal of either particles()
	// or for (auto& v : evt.vertices() ){ for (auto& p : v->particles_out() ) {}}
	// Then we would also have to overwrite
	// EventMC::BeamLepton(), EventMC::BeamHadron(), EventMC::ExchangeBoson(), EventMC::ScatteredLepton().
	// and worry about the logic when the file is used.

	// Alternatively, we can enforce standard order
	// BeamLepton, BeamHadron, ScatteredLepton, ExchangeBoson
	// That way we only need to do this logic once.
	// Need to keep track of daughter-mother relationship for that.

	int particleindex;
	particleindex = 1;
	// Can't use GenParticle::children() because they don't have indices assigned yet
	// map each HepMC particle onto its corresponding particleindex
	std::map < HepMC3::GenParticlePtr, int > hepmcp_index;
	hepmcp_index.clear();
	
	// start with the beam
	// Boring determination of the order:
	auto beams = evt.beams();
	assert ( beams.size() == 2 );
	auto hadron = beams.at(0); // or nucleon
	auto lepton = beams.at(1);

	// Before going on to handle normal events, 
	// catch special beam-gas-only events
	if ( hadron->pid()==22 || lepton->pid()==22 ){
	  if ( evt.particles().size() > 2 ){
	    cout << "Warning: Treating the event as beam gas but found "
		 << evt.particles().size() << " particles" << endl;
	  }
	  auto photon=hadron; // for clarity use the right name
	  if ( lepton->pid()==22 ) {
	    std::swap (lepton, photon);
	  }

	  
	  // Still need beam particles (otherwise fun4all is confused)
	  auto beam_e_mom = lepton->momentum() + photon->momentum();
	  // cout << beam_e_mom.x() << " "  
	  //      << beam_e_mom.y() << " "
	  //      << beam_e_mom.z() << " "
	  //      << beam_e_mom.e() << endl;
	  auto beam_e = std::make_shared<HepMC3::GenParticle>( beam_e_mom, 11, 21 );
	  auto v = lepton->production_vertex();
	  v->add_particle_in  (beam_e);

	  // make up a proton
	  auto pz_p = -10.0 *beam_e_mom.z();
	  auto e_p = sqrt ( pow(pz_p,2) + pow(0.938,2) );
	  HepMC3::FourVector beam_p_mom ( 0,0,pz_p, e_p);
	  auto beam_p = std::make_shared<HepMC3::GenParticle>( beam_p_mom, 2212, 21 );

	  HandleHepmcParticle( beam_e, hepmcp_index, particleindex, mEvent );
	  HandleHepmcParticle( beam_p, hepmcp_index, particleindex, mEvent );
	  
	  HandleHepmcParticle( lepton, hepmcp_index, particleindex, mEvent );
	  HandleHepmcParticle( photon, hepmcp_index, particleindex, mEvent );

	  // x-setion etc.
	  UpdateRuninfo( mObjectsToWriteAtTheEnd, evt );
	  // We're done. Avoid FinishEvent()
	  return true;
	}
	     
	// Find the lepton
	auto pdgl = TDatabasePDG::Instance()->GetParticle( lepton->pid() );
	auto pdgh = TDatabasePDG::Instance()->GetParticle( hadron->pid() );
	// Sigh. TParticlePDG uses C strings for particle type. I refuse to use strcmp
	// Could test individual pid's instead.
	bool b0islepton = (string(pdgl->ParticleClass()) == "Lepton");
	bool b1islepton = (string(pdgh->ParticleClass()) == "Lepton");
	if ( b0islepton == b1islepton ){
	  // exactly one of them should be a lepton.
	  throw std::runtime_error ("Exactly one beam should be a lepton - please contact the authors for ff or hh beams");
	}
	if ( !b0islepton ) {
	  std::swap (lepton, hadron);
	}
	// careful, don't try to use pdg[lh], b[01]islepton from here on out;
	// they don't get swapped along

	// now find the scattered lepton and potentially the gamma
	// some processes (ff2ff) don't have a gamma!
	HepMC3::GenParticlePtr scatteredlepton=0;
	HepMC3::GenParticlePtr photon=0;

	// we can handle one or two only
	if ( lepton->children().size() >2 || lepton->children().size()==0 ){
	  cerr << "electron has " << lepton->children().size() << " daughters." << endl;
	  throw std::runtime_error ("Wrong number of lepton daughters; should be 1 (lepton) or 2 (lepton+boson).");
	}

	// if one exists, it's a daughter of the beam lepton, and its sibling is the lepton (which isn't necessarily final yet though)
	if ( lepton->children().size() == 2 ){
	  scatteredlepton = lepton->children().at(0);
	  photon = lepton->children().at(1);
	  if ( scatteredlepton->pid() != 22 && photon->pid() != 22 ){
	    cerr << "lepton child 1 pid = " << scatteredlepton->pid() << endl;
	    cerr << "lepton child 2 pid = " << photon->pid() << endl;
	    throw std::runtime_error ("Found two beam lepton daughters, none or both of them a photon.");
	  }
	  if ( photon->pid() != 22 ){
	    std::swap ( photon, scatteredlepton );
	  }
	}

	// no exchange boson
	if ( lepton->children().size() == 1 ){
	  scatteredlepton = lepton->children().at(0);
	  auto pdgs = TDatabasePDG::Instance()->GetParticle( scatteredlepton->pid() );
	  if ( !TString(pdgs->ParticleClass()).Contains("Lepton") ){
	    throw std::runtime_error ("Found one beam lepton daughter, and it is not a lepton.");
	  }

	  // We could play games and make one up:
	  // auto photonmom = lepton->momentum() - scatteredlepton->momentum();
	  // photon = std::make_shared<HepMC3::GenParticle>( photonmom, 22, 13 );
	  // photon->set_momentum(photonmom);
	  // photon->set_pid( 22 );
	  // // And add to the vertex (meaning the photon now has the lepton as a mother)
	  // scatteredlepton->production_vertex()->add_particle_out(photon);
	  // But the cleaner solution seems to be to save a 0 particle in the photon spot
	  photon=std::make_shared<HepMC3::GenParticle>();
	   // we'll catch this in the ParticleIdentifier
	  photon->set_pid( 0 );
	  photon->set_status( 0 );
	  // Need to give a production vertex to avoid segfaulting
	  scatteredlepton->production_vertex()->add_particle_out(photon);
	}

	// Follow the lepton
	// Assumptions (otherwise this code will break...)
	// - it can onlhy branch into l -> l + X, where X is not a lepton
	// -> we can traverse children and immediately follow any lepton,
	//    we won't miss a "second" lepton
	// -> we can go until we run out of children, the lepton won't decay
	//    or we can go until the status is final, the result should be the same
	// There's no clear-cut final() method though, so we go through the end and hope!
	auto spid = scatteredlepton->pid();
	// avoid endless loop
	bool foundbranch=true;
	while ( scatteredlepton->children().size() > 0 ){
	  if ( !foundbranch ){
	    throw std::runtime_error ("none of this lepton's children is a lepton.");
	  }
	  foundbranch=false;
	  for ( auto& c : scatteredlepton->children() ){
	    if ( c->pid() == spid ){
	      // found the correct branch,
	      // update and break out of for loop,
	      // resume while loop with new candidate
	      scatteredlepton = c;
	      foundbranch=true;
	      break;
	    }
	  }
	}

	if ( scatteredlepton->pid() != lepton->pid() ){
	  cerr << "lepton  pid = " << lepton->pid() << endl;
	  cerr << "scattered lepton  pid = " << scatteredlepton->pid() << endl;
	  throw std::runtime_error ("Scattered lepton pid mismatch.");
	}

	// Now add all four in the right order
	HandleHepmcParticle( lepton, hepmcp_index, particleindex, mEvent );
	HandleHepmcParticle( hadron, hepmcp_index, particleindex, mEvent );
	HandleHepmcParticle( scatteredlepton, hepmcp_index, particleindex, mEvent );
	HandleHepmcParticle( photon, hepmcp_index, particleindex, mEvent );

	// Now go over all vertices and handle the rest.
	// Note that by default this could double-count what we just did
	// Instead of trying to find out which vertex to skip, just use the lookup table
	// (inside HandleHepmcParticle) to avoid this.
	for (auto& v : evt.vertices() ){
	  for (auto& p : v->particles_out() ) {
	    HandleHepmcParticle( p, hepmcp_index, particleindex, mEvent );
	  }
	}

	// Now the map has built up full 1-1 correspondence between all hepmc particles
	// and the ParticleMC index.
	// So we can loop over the particles again, find their parents and offspring, and map accordingly.
	// We explicitly take advantage of particleid = Event entry # +1
	// If that changes, we'll need to maintain a second map
	// Note: the beam proton appears twice; that's consistent with the behavior of pythia6
	for (auto& v : evt.vertices() ){
	  for (auto& p : v->particles_out() ) {
	    // corresponding ParticleMC is at
	    int treeindex = hepmcp_index[p]-1;
	    assert ( treeindex >=0 ); // Not sure if that can happen. If it does, we could probably just continue;
	    auto track = mEvent->GetTrack(treeindex);

	    // collect all parents
	    vector<int> allparents;
	    for (auto& parent : p->parents() ) {
	      allparents.push_back( hepmcp_index[parent] );
	    }
	    // orig and orig1 aren't very intuitively named...
	    // For pythia6, orig is the default, and orig1 seems purely a placeholder
	    // trying to mimick that here.
	    if ( allparents.size() == 1){
	      track->SetParentIndex( allparents.at(0) );
	    }
	    if ( allparents.size() >= 2){
	      // smallest and highest are stored
	      track->SetParentIndex( *min_element(allparents.begin(), allparents.end() ) );
	      track->SetParentIndex1( *max_element(allparents.begin(), allparents.end() ) );
	    }

	    // same for children
	    vector<int> allchildren;
	    for (auto& child : p->children() ) {
	      allchildren.push_back( hepmcp_index[child] );
	    }
	    if ( allchildren.size() == 1){
	      track->SetChild1Index( allchildren.at(0) );
	    }
	    if ( allchildren.size() >= 2){
	      // smallest and highest are stored
	      track->SetChild1Index( *min_element(allchildren.begin(), allchildren.end() ) );
	      track->SetChildNIndex( *max_element(allchildren.begin(), allchildren.end() ) );
	    }
	  }
	}
	// x-setion etc.
	UpdateRuninfo( mObjectsToWriteAtTheEnd, evt );
      }  // if

      auto finished = FinishEvent(); // 0 is success
      return (finished==0);
    }  // try
    catch(std::exception& error) {
      std::cerr << "Exception building particle: " << error.what() << std::endl;
      return false;
    }
  }

  template<>
  erhic::EventHepMC* EventFromAsciiFactory<erhic::EventHepMC>::Create()
  {
    TProcessIdObjectCount objectCount;
    mEvent.reset(new erhic::EventHepMC());
    if (!AddParticle()) {
      mEvent.reset(nullptr);
    }  // if

    return mEvent.release();
  }

  void HandleHepmcParticle( const HepMC3::GenParticlePtr& p, std::map < HepMC3::GenParticlePtr, int >& hepmcp_index, int& particleindex, std::unique_ptr<erhic::EventHepMC>& mEvent ){
    // do nothing if we already used this particle
    auto it = hepmcp_index.find(p);
    if ( it != hepmcp_index.end() ) return;

    // Create the particle
    ParticleMC particle;
    auto v = p->production_vertex();
    TVector3 vertex;
    if ( v ) vertex = TVector3(v->position().x(),v->position().y(),v->position().z());

    particle.SetVertex(vertex);
    // takes care of  xv, yv, zv;

    TLorentzVector lovec(p->momentum().x(),p->momentum().y(),p->momentum().z(),p->momentum().e());
    particle.Set4Vector(lovec);
    // takes care of  E, px, py, pz, m, and
    // derived quantities: pt, p, rapidity, eta, theta, phi

    // fill up the missing parts
    particle.SetId( p->pid() );
    particle.SetStatus(p->status());

    // Index: Runs from 1 to N.
    particle.SetIndex ( particleindex );

    // remember this HepMC3::GenParticlePtr <-> index connection
    hepmcp_index[p] = particleindex;

    particleindex++;

    particle.SetEvent(mEvent.get());
    mEvent->AddLast(&particle);
  }

  // -----------------------------------------------------------------------
  void UpdateRuninfo( std::vector<VirtualEventFactory::NamedObjects>& mObjectsToWriteAtTheEnd, 
		      const HepMC3::GenEvent& evt ){
    // Wasteful, since it's only written out for the final event, but
    // this class doesn't know when that happens
    TString xsec = ""; 
    TString xsecerr = "";
    if ( auto xsecinfo=evt.cross_section() ){
      // not all generators have this info
      xsec += xsecinfo->xsec();
      xsecerr += xsecinfo->xsec_err();
    }
    TObjString* Xsec;
    TObjString* Xsecerr;
    
    // // not saved in the examples I used.
    // TString trials = ""; trials+= evt.cross_section()->get_attempted_events();
    // TString nevents = ""; nevents += evt.cross_section()->get_accepted_events();
    // TObjString* Trials;
    // TObjString* Nevents;
    
    if ( mObjectsToWriteAtTheEnd.size() == 0 ){
      Xsec = new TObjString;
      Xsecerr = new TObjString;
      mObjectsToWriteAtTheEnd.push_back ( VirtualEventFactory::NamedObjects( "crossSection", Xsec) );
      mObjectsToWriteAtTheEnd.push_back ( VirtualEventFactory::NamedObjects( "crossSectionError", Xsecerr) );
      // Trials = new TObjString;
      // Nevents = new TObjString;
      // mObjectsToWriteAtTheEnd.push_back ( VirtualEventFactory::NamedObjects( "nTrials", Trials) );
      // mObjectsToWriteAtTheEnd.push_back ( VirtualEventFactory::NamedObjects( "nEvents", Nevents) );
    }
    Xsec = (TObjString*) mObjectsToWriteAtTheEnd.at(0).second;
    Xsec->String() = xsec;
    Xsecerr = (TObjString*) mObjectsToWriteAtTheEnd.at(1).second;
    Xsecerr->String() = xsecerr;
    // Trials = (TObjString*) mObjectsToWriteAtTheEnd.at(2).second;
    // Trials->String() = trials;
    // Nevents = (TObjString*) mObjectsToWriteAtTheEnd.at(3).second;
    // Nevents->String() = nevents;
  }
  
}  // namespace erhic

namespace {

  // Need this to generate the CINT code for each version
  erhic::EventFromAsciiFactory<erhic::EventHepMC> eh;
}  // namespace
