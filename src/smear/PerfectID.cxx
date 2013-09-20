/**
 \file
 Implementation of class Smear::PerfectId.
 
 \author    Michael Savastio
 \date      2011-08-19
 \copyright 2011 Brookhaven National Lab
 */

#include "eicsmear/smear/PerfectID.h"

#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <TDatabasePDG.h>
#include <TParticlePDG.h>

namespace {

// Returns the name of the particle species with the PDG code,
// or the integer converted to a string if the PDG information
// cannot be found.
std::string particleName(int pdg) {
  std::stringstream stream;
  TParticlePDG* p = TDatabasePDG::Instance()->GetParticle(pdg);
  if (p) {
    stream << p->GetName();
  } else {
    stream << pdg;
  }  // if
  return stream.str();
}

}  // anonymous namespace

namespace Smear {

PerfectID::PerfectID(const std::vector<Int_t>& pdg)
: mPdg(pdg.begin(), pdg.end()) {
}

PerfectID::~PerfectID() {
}

void PerfectID::Smear(const erhic::VirtualParticle& in, ParticleMCS& out) {
  // If the PDG list is empty, always copy PDG code.
  // Otherwise, copy the PDG code if it is in the list.
  bool copy = mPdg.empty() || mPdg.find(in.Id()) != mPdg.end();
  if (copy) {
    out.SetId(in.Id());
  }  // if
}

PerfectID* PerfectID::Clone(const char*) const {
  return new PerfectID(*this);
}

void PerfectID::Print(Option_t* /* option */) const {
  std::stringstream stream;
  stream << "Copies PDG ID for ";
  if (mPdg.empty()) {
    stream << "all particles";
  } else {
    // List the PDG codes. Insert the first one into the stream,
    // then loop from the second (if it exists) onward so we
    // can delimit with commas.
    std::set<Int_t>::const_iterator iter = mPdg.begin();
    stream << particleName(*iter);
    for (++iter; iter != mPdg.end(); ++iter) {
      stream << ", " << particleName(*iter);
    }  // for
  }  // if
  std::cout << stream.str() << std::endl;
}

void PerfectID::Insert(Int_t i) {
  mPdg.insert(i);
}

}  // namespace Smear
