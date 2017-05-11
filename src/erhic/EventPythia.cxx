/**
 \file
 Implementation of class erhic::EventRapgap.
 
 \author    Thomas Burton
 \date      2011-08-31
 \copyright 2011 Brookhaven National Lab
 */

#include "eicsmear/erhic/EventPythia.h"

#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace erhic {

EventPythia::EventPythia(const std::string& /* Unused */)
: nucleon(std::numeric_limits<Int_t>::max())
, tgtparton(std::numeric_limits<Int_t>::max())
, beamparton(std::numeric_limits<Int_t>::max())
, genevent(-1)
, xtgtparton(NAN)
, xbeamparton(NAN)
, thetabeamparton(NAN)
, leptonphi(NAN)
, F1(NAN)
, sigma_rad(NAN)
, t_hat(NAN)
, u_hat(NAN)
, Q2_hat(NAN)
, SigRadCor(NAN)
, EBrems(NAN)
, photonflux(NAN)
, trueY(NAN)
, trueQ2(NAN)
, trueX(NAN)
, trueW2(NAN)
, trueNu(NAN)
, F2(NAN)
, R(NAN)
, pt2_hat(NAN)
, sHat(NAN) {
}

EventPythia::~EventPythia() { }

bool EventPythia::Parse(const std::string& line) {
  static std::stringstream ss;
  ss.str("");
  ss.clear();
  ss << line;
  ss >>
  number >> number >>  // Skip first int in the line
  genevent >> process >> nucleon >> tgtparton >> xtgtparton >>
  beamparton >> xbeamparton >> thetabeamparton >> trueY >> trueQ2 >>
  trueX >> trueW2 >> trueNu >> leptonphi >> sHat >> t_hat >> u_hat >>
  pt2_hat >> Q2_hat >> F2 >> F1 >> R >> sigma_rad >> SigRadCor >> EBrems >>
  photonflux >> nTracks;
  // Protect against errors in the input file or the stream
  return !ss.fail();
}

// Look for the scattered lepton in the event record.
// This is the first (only?) particle that matches the following:
//  1) pdg code equals that of incident lepton beam.
//  2) status code is 1 i.e. it's a stable/final-state particle.
//  3) the parent is track three (counting from 1).
const ParticleMC* EventPythia::ScatteredLepton() const {
  // Look for the lepton beam to get the species.
  // If we don't get it we can't find the scattered
  // lepton so return NULL.
  const VirtualParticle* beam = BeamLepton();
  if (!beam) {
    return NULL;
  }  // if
  const int species = beam->Id().Code();
  // Get the final state particles and search them for
  // the scattered lepton.
  std::vector<const VirtualParticle*> final;
  FinalState(final);
  std::vector<const VirtualParticle*>::const_iterator iter;
  for (iter = final.begin(); iter != final.end(); ++iter) {
    // We already know the particle is final state, so
    // check its species and parent index.
    if ((*iter)->GetParentIndex() == 3 &&
       (*iter)->Id().Code() == species) {
      // Found it, cast to required particle type and return.
      return static_cast<const ParticleMC*>(*iter);
    }  // if
  }  // for
  // No luck, couldn't find the scattered lepton.
  return NULL;
}

  EventBeagle::EventBeagle(const std::string& str/* Unused */): 
    EventPythia(str)
    , lepton(std::numeric_limits<Int_t>::max()) 
    , Atarg(std::numeric_limits<Int_t>::max())
    , Ztarg(std::numeric_limits<Int_t>::max())
    , pzlep(NAN)
    , pztarg(NAN)
    , pznucl(NAN)
    , crang(NAN)
    , crori(NAN)
    , b(NAN)
    , Phib(NAN)
    , Thickness(NAN)
    , ThickScl(NAN)
    , Ncollt(std::numeric_limits<Int_t>::max())
    , Ncolli(std::numeric_limits<Int_t>::max())
    , Nwound(std::numeric_limits<Int_t>::max())
    , Nwdch(std::numeric_limits<Int_t>::max())
    , Nnevap(std::numeric_limits<Int_t>::max())
    , Npevap(std::numeric_limits<Int_t>::max())
    , Aremn(std::numeric_limits<Int_t>::max())
    , NINC(std::numeric_limits<Int_t>::max())
    , NINCch(std::numeric_limits<Int_t>::max())
    , d1st(NAN)
    , davg(NAN)
    , pxf(NAN)
    , pyf(NAN)
    , pzf(NAN)
    , Eexc(NAN)
    , RAevt(NAN)
    , User1(NAN)
    , User2(NAN)
    , User3(NAN)
  {
   
  } // EventBeagle::EventBeagle()

  EventBeagle::~EventBeagle() { }

bool EventBeagle::Parse(const std::string& line) {
  static std::stringstream ss;
  ss.str("");
  ss.clear();
  ss << line;
  ss >>
//  number >> number >>  // Skip first int in the line
//  genevent >> process >> nucleon >> tgtparton >> xtgtparton >>
//  beamparton >> xbeamparton >> thetabeamparton >> trueY >> trueQ2 >>
//  trueX >> trueW2 >> trueNu >> leptonphi >> sHat >> t_hat >> u_hat >>
//  pt2_hat >> Q2_hat >> F2 >> F1 >> R >> sigma_rad >> SigRadCor >> EBrems >>
//  photonflux >> nTracks;
  number >> number >>  // Skip first int in the line
  genevent >> 
  lepton >> Atarg >> Ztarg >> pzlep >> pztarg >> pznucl >> crang >> crori >>//added variables
  process >> nucleon >> tgtparton >> xtgtparton >>
  beamparton >> xbeamparton >> thetabeamparton >> trueY >> trueQ2 >>
  trueX >> trueW2 >> trueNu >> leptonphi >> sHat >> t_hat >> u_hat >>
  pt2_hat >> Q2_hat >> F2 >> F1 >> R >> sigma_rad >> SigRadCor >> EBrems >>
  photonflux >> 
  b >> Phib >> Thickness >> ThickScl >> Ncollt >> Ncolli >> Nwound >> Nwdch >> Nnevap >> Npevap >> Aremn >> //added variables
  NINC >> NINCch >> d1st >> davg >> pxf >> pyf >> pzf >> Eexc >> RAevt >> User1 >> User2 >> User3 >> //added variables
  nTracks;
  // Protect against errors in the input file or the stream
  return !ss.fail();
}
}  // namespace erhic
