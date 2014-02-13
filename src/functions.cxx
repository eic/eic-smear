/**
 \file
 Global function implementations.
 
 \author    Thomas Burton
 \date      2011-07-07
 \copyright 2011 Brookhaven National Lab
 */

#include "eicsmear/functions.h"

#include <cmath>
#include <fstream>
#include <set>
#include <sstream>
#include <string>

#include <TLorentzVector.h>
#include <TMath.h>
#include <TParticlePDG.h>
#include <TVector2.h>
#include <TVector3.h>

#include "eicsmear/erhic/EventMC.h"
#include "eicsmear/erhic/ParticleMC.h"

/**
 Returns the first non-blank character in a line.
 Returns \0 if there are no non-blank characters in the line.
 */
char getFirstNonBlank(const std::string& line) {
  char first('\0');  // NULL character by default
                      // Find the index first non-blank character.
  size_t index = line.find_first_not_of(" \t");

  if (index != std::string::npos) {
    first = line.at(index);
  }  // if

  return first;
}

/**
 Calculate the hadron azimuthal angle around the virtual photon direction with
 respect to the lepton scattering plane in the proton rest frame.
 We use the HERMES convention, returning an angle in the range [0,2pi].
 The vectors passed as arguments should already be boosted to the proton rest
 frame. Incident and scattered leptons, incident protons and virtual photons
 all return -999.
 */
double
computeHermesPhiH(const TLorentzVector& hadronInPrf,
                  const TLorentzVector& leptonInPrf,
                  const TLorentzVector& photonInPrf) {
  const TVector3& q = photonInPrf.Vect();
  const TVector3& k = leptonInPrf.Vect();
  const TVector3& Ph = hadronInPrf.Vect();
  const TVector3& qCrossK = q.Cross(k);
  const TVector3& qCrossPh = q.Cross(Ph);
  // Used to determine spin direction:
  double qCrossKDotPh = qCrossK.Dot(Ph);
  if (::fabs(qCrossKDotPh) < 1.0e-5) return -999.0;

  // First get phi in the range [0,pi]:
  double phi = TMath::ACos(qCrossK.Dot(qCrossPh)
                           / qCrossK.Mag()
                           / qCrossPh.Mag());

  // Handle spin direction to get result in range [-pi,pi] by
  // scaling by [-1,+1], depending on hemisphere.
  phi *= qCrossKDotPh / ::fabs(qCrossKDotPh);

  return TVector2::Phi_0_2pi(phi);
}

//
// class EventToDot.
//

struct Pair {
  int a;
  int b;
  Pair(int x, int y) : a(x), b(y) { }
  // Less-than operator required for entry in sorted container.
  bool operator<(const Pair& other) const {
    if (other.a == a) return b < other.b;
    return a < other.a;
  }
};


void EventToDot::Generate(const erhic::EventMC& event,
                          const std::string& name) const {
  std::ofstream file(name.c_str());
  std::ostringstream oss;

  file << "digraph G {" << std::endl;
  oss << "   label = \"Event " << event.GetN() << "\"';";
  file << oss.str() << std::endl;

  std::set<Pair> pairs;

  // Keep track of which particles we use in the diagram
  // so we can set node attributes at the end.
  // Use a set to avoid worrying about adding duplicates.
  std::set<const erhic::ParticleMC*> used;

  // Loop over all tracks in the event and accumulate indices giving
  // parent/child relationships.
  // We look from parent->child and child->parent because they
  // aren't always fully indexed both ways.
  for (unsigned i(0); i < event.GetNTracks(); ++i) {
    // Check parent of particle
    const erhic::ParticleMC* parent = event.GetTrack(i)->GetParent();
    if (parent) {
      pairs.insert(Pair(parent->GetIndex(), event.GetTrack(i)->GetIndex()));
      used.insert(parent);
      used.insert(event.GetTrack(i));
    }  // if
    // Check children of particle
    for (unsigned j(0); j < event.GetTrack(i)->GetNChildren(); ++j) {
      pairs.insert(Pair(event.GetTrack(i)->GetIndex(),
                        event.GetTrack(i)->GetChild(j)->GetIndex()));
      used.insert(event.GetTrack(i));
      used.insert(event.GetTrack(i)->GetChild(j));
    }  // for
  }  // for
  // Insert relationships between particles.
  for (std::set<Pair>::const_iterator i = pairs.begin();
      i != pairs.end();
      ++i) {
    const erhic::ParticleMC* a = event.GetTrack(i->a - 1);
    const erhic::ParticleMC* b = event.GetTrack(i->b - 1);
    oss.str("");
    oss << "   " << a->GetIndex() << " -> " <<
    b->GetIndex();
    file << oss.str() << std::endl;
  }  // for
  file << "#   Node attributes" << std::endl;
  // Apply labels, formatting to used particles.
  for (std::set<const erhic::ParticleMC*>::const_iterator i = used.begin();
      i != used.end();
      ++i) {
    std::string shape("ellipse");
    if ((*i)->GetStatus() == 1) {
      shape = "box";
    }  // if
    if ((*i)->GetIndex() < 3) {
      shape = "diamond";
    }  // if
    oss.str("");
    oss << "   " << (*i)->GetIndex() << " [label=\"" << (*i)->GetIndex() << " ";
    if ((*i)->Id().Info()) {
      oss << (*i)->Id().Info()->GetName();
    } else {
      oss << "unknown PDG " << (*i)->Id().Code();
    }  // if
    oss << "\", shape=" << shape << "];";
    file << oss.str() << std::endl;
  }  // for
  file << "}";
}
