/**
 \file
 Implementation of class Smear::ParticleID.
 
 \author    Michael Savastio
 \date      2011-08-12
 \copyright 2011 BNL. All rights reserved.
 */

#include "eicsmear/smear/ParticleID.h"

#include <sstream>
#include <string>
#include <vector>

namespace Smear {

ParticleID::ParticleID()
: Ran(0)
, PMatPath("PIDMatrix.dat")
, bUseMC(false) {
  ReadP(PMatPath);
}

ParticleID::ParticleID(TString filename)
: Ran(0)
, PMatPath(filename)
, bUseMC(false) {
  ReadP(PMatPath);
}

ParticleID::~ParticleID() {
}

void ParticleID::Speak() {
  std::cout.setf(std::ios::fixed);
  std::cout.precision(6);
  for (unsigned i(0); i < PMax.size(); i++) {
    for (unsigned k(0); k < FalseIdent.size(); k++) {
      std::cout << 1 << '\t' << i + 1 << '\t' << PMin.at(i) << '\t' <<
      PMax.at(i) << '\t' << k + 1;
      for (unsigned j(0); j < TrueIdent.size(); j++) {
        std::cout << '\t' << PMatrix.at(i).at(j).at(k);
      }  // for
      std::cout << std::endl;
    }  // for
  }  // for
}

void ParticleID::SetPMatrixSize() {
  PMatrix.resize(PMin.size());
  for (unsigned i(0); i < PMatrix.size(); i++) {
    PMatrix.at(i).resize(TrueIdent.size());
    for (unsigned j(0); j < PMatrix.at(i).size(); j++) {
      PMatrix.at(i).at(j).resize(FalseIdent.size());
    }  // for
  }  // for
}

void ParticleID::SetupProbabilityArray() {
  double t = 0.;
  Range.resize(PMatrix.size());
  for (unsigned i(0); i < Range.size(); i++) {
    Range.at(i).resize(PMatrix.at(i).size());
    for (unsigned j(0); j < Range.at(i).size(); j++) {
      Range.at(i).at(j).resize(PMatrix.at(i).at(j).size());
      for (unsigned k = 0; k < Range.at(i).at(j).size(); k++) {
        t += PMatrix.at(i).at(j).at(k);
        Range.at(i).at(j).at(k) = t;
      }  // for
      t = 0.;
    }  // for
  }  // for
}

int ParticleID::Wild(int pbin, int trueID) {
  const double r = Ran.Rndm();
  const int k = InListOfTrue(trueID);
  int falseid(-999999999);
  // Get the cumulative probability values for this momentum bin
  // and true ID
  const std::vector<double>& values = Range.at(pbin).at(k);
  for (unsigned i(0); i < values.size(); i++) {
    if (0 == i) {
      if (r < values.at(i)) {
        falseid = FalseIdent.at(i);
      }  // if
    } else if (r > values.at(i-1) && r < values.at(i)) {
      falseid = FalseIdent.at(i);
    }  // if
  }  // for
  return falseid;
}

int ParticleID::InListOfTrue(int ID) {
  for (unsigned i(0); i < TrueIdent.size(); i++) {
    if (TrueIdent.at(i) == abs(ID)) {
      return i;
    }  // if
  }  // for
  return -1;
}

int ParticleID::InListOfFalse(int ID) {
  for (unsigned i(0); i < FalseIdent.size(); i++) {
    if (FalseIdent.at(i) == abs(ID)) {
      return i;
    }  // if
  }  // for
     // This is to accomodate the old HERMES format
  if (ID < 5 && ID > 0) {
    return ID - 1;
  }  // for
  return -1;
}

void ParticleID::Clear(Option_t* /* option */) {
  TrueIdent.clear();
  FalseIdent.clear();
  PMin.clear();
  PMax.clear();
  PMatrix.clear();
  Range.clear();
}

void ParticleID::ReadP(TString filename) {
  Clear("");
  // Open the file and check for errors
  std::ifstream Qfile;
  Qfile.open(filename);
  if (!Qfile.good()) {
    std::cerr << "Error in ParticleID: unable to open file "
    << filename << std::endl;
    return;
  }  // if
     // Returns true if s begins with pattern.
  struct StartsWith {
    bool operator()(const std::string& s,
                    const std::string& pattern) const {
      return s.find(pattern, 0) == 0;
    }
  }
  starts;
  std::stringstream ss;
  std::string line, dummy;
  bool gotTrue(false), gotFalse(false), gotBins(false);
  while (std::getline(Qfile, line).good()) {
    // Strip leading whitespace
    line.erase(0, line.find_first_not_of(" \t"));
    // Feed the line to the stringstream, clearing existing contents first
    ss.str("");  // Remove contents
    ss.clear();  // Clear flags
    ss << line;
    int tmpint(0);
    // Read true particle IDs
    if (starts(line, "!T")) {
      ss >> dummy;
      while ((ss >> tmpint).good()) {
        TrueIdent.push_back(tmpint);
      }  // while
      gotTrue = !TrueIdent.empty();
    } else if (starts(line, "!F")) {
      // Read misidentified particle IDs
      ss >> dummy;
      while ((ss >> tmpint).good()) {
        FalseIdent.push_back(tmpint);
      }  // while
      gotFalse = !FalseIdent.empty();
    } else if (starts(line, "!P")) {
      // Read number of momentum bins
      int pbin;
      ss >> dummy >> pbin;
      PMin.resize(pbin);
      PMax.resize(pbin);
      gotBins = !PMin.empty();
    } else if (starts(line, "1") || starts(line, "2") || starts(line, "3")) {
      // Read the probabilities for this momentum bin/true ID/false ID
      if (!(gotTrue && gotFalse && gotBins)) {
        std::cerr << "Error in ParticleID: " <<
        "P matrix input file has bad or missing format lines.\n";
        return;
      }  // if
      int table, pbin, pid;
      double pmin, pmax, p1, p2, p3;
      ss >> table >> pbin >> pmin >> pmax >> pid >> p1 >> p2 >> p3;
      if ((unsigned)pbin > PMin.size()) {
        std::cerr << "Error in ParticleID: " <<
        "Out of bounds momentum bin listing.\n";
        return;
      }  // if
      pbin -= 1;  // Go from [1, N] to [0, N) for vector index range
      PMin.at(pbin) = pmin;
      PMax.at(pbin) = pmax;
      pid = InListOfFalse(pid);
      if ((unsigned)pid > FalseIdent.size()) {
        std::cerr << "Error in ParticleID: " <<
        "P matrix has bad particle listing.\n";
        return;
      }  // if
      if (PMatrix.empty()) {
        SetPMatrixSize();
      }  // if
      // Set identification probabilities
      if (1 == table) {
        PMatrix.at(pbin).at(0).at(pid) = p1;
        PMatrix.at(pbin).at(1).at(pid) = p2;
        PMatrix.at(pbin).at(2).at(pid) = p3;
      }  // if
    }  // if
  }  // while
  SetupProbabilityArray();
}

void ParticleID::Smear(const erhic::VirtualParticle& prt,
                       ParticleMCS& prtOut) {
  double momentum(0.);
  if (bUseMC) {
    momentum = prt.GetP();
  } else {
    momentum = prtOut.p;
  }  // if
  const int pid = prt.Id();
  if (InListOfTrue(pid) != -1 && Accept.Is(prt)) {
    for (unsigned i(0); i < PMin.size(); i++) {
      if (momentum > PMin.at(i) && momentum < PMax.at(i)) {
        // Generated ID is always positive.
        // Keep same sign as input PID i.e. no error in charge sign
        if (pid > 0) {
          prtOut.id = Wild(i, pid);
        } else {
          prtOut.id = -Wild(i, pid);
        }  // if
      }  // if
    }  // for
  }  // if
}

void ParticleID::Print(Option_t* /* option */) const {
  std::cout << "ParticleID using " << PMatPath << std::endl;
}

}  // namespace Smear
