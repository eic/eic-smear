/**
 \file
 Declaration of class Smear::ParticleID.
 
 \author    Michael Savastio
 \date      2011-08-12
 \copyright 2011 BNL. All rights reserved.
 */

#ifndef INCLUDE_EICSMEAR_SMEAR_PARTICLEID_H_
#define INCLUDE_EICSMEAR_SMEAR_PARTICLEID_H_

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include <TF1.h>
#include <TF2.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>

#include "eicsmear/smear/Device.h"
#include "eicsmear/erhic/Kinematics.h"
#include "eicsmear/erhic/Particle.h"
#include "eicsmear/erhic/ParticleIdentifier.h"
#include "eicsmear/smear/Smear.h"
#include "eicsmear/smear/Smearer.h"
#include "eicsmear/erhic/VirtualParticle.h"

namespace Smear {

/**
 This structure is used to generate particle ID.
 
 The input file containing the P matrix must begin with format
 lines beginning with "!T", "!F" and "!P".  For example
 
 !T 211 321 2212
 !F 211 321 2212 0
 !P 15
 
 The first line tells the PID generator which particles false ID
 will be generated for. The second line tells it which particles the
 particles in the first line can be identified as. The third line tells
 it how many momentum bins there are.
 
 Lines of data should appear as:
 
 1 pbinNumber  pmin pmax  FalseID  P1  P2  P3
 
 The line must begin with a 1 to be read.
 pbinNumber is the number of the momentum bin.
 pmin and pmax are the bounds of the momentum bin.
 FalseID is the ID that the particle will be mis-identified as,
 this should be a PDG particle code (unless it is an old Hermes file).
 P1 is the probability that the first particle appearing in the "!T"
 line will be misidentified as FalseID. Likewise for P2 and P3.
 
 \todo Implement data hiding
 \remark Why is the bitwise AND assignment operator&= overloaded?!
 */
struct ParticleID : public Smearer {
  /**
   Default constructor.
   Default path for probability matrix is "PIDMatrix.dat".
   Also, by default, the particle ID uses smeared kinematic
   variables as input, rather than Monte Carlo values.
   */
  ParticleID();

  /**
   Constructor.
   Initialise the particle misidentification matrix from the named file.
   */
  explicit ParticleID(TString filename);

  /**
   Destructor.
   */
  virtual ~ParticleID();

  /** 
   Set the path to the file containing the particle ID probability matrix,
   and read it in.
   */
  void SetPMatrixPath(TString);

  /** 
   If true, the ParticleID will use the original Monte Carlo values
   to generate ID's, rather than smeared values.
   Default is false.
   */
  void SetPIDUseMC(bool useMc);

  /**
   Return the TRandom3 instance used by the ParticleID.
   */
  TRandom3& GetRandomGenerator();

  /**
   Set the seed of the random number generator instance used by the 
   ParticleID.
   See ROOT TRandom3 documentation for more information.
   */
  void SetRanSeed(int seed);

  /**
   Set the ParticleID instance to have the same acceptance as the device.
   This includes all acceptance zones as well as the list of specific
   particles to be detected (if there is one).
   */
  void GetAcceptanceFromDevice(const Device&);

  /**
   Resize the PMatrix, (# of pbins, # of true, # of false).
   */
  void SetPMatrixSize();

  /**
   Setup the range array used by wild.
   
   Range is used to set up zones in [0,1],
   if a random number falls in one of the zones, the associated ID is used.
   */
  void SetupProbabilityArray();

  /**
   Returns a copy of this ParticleID.
   Inherited from TObject.
   The const char* argument has no effect.
   */
  virtual ParticleID* Clone(const char* = "") const;

  /**
   Randomly generates a false ID based on the current, a momentum bin
   and a true particle ID.
   */
  int Wild(int pbin, int trueID);

  int InListOfTrue(int ID);

  int InListOfFalse(int ID);

  /**
   Read in a P matrix and set up the ParticleID to be ready to generate.
   See the documentation for the detector function ReadPIDMatrix.
   */
  void ReadP(TString filename);

  /**
   Generates particle ID if the particle is in the list of particles
   to be identified, falls within the momentum range of the PMatrix
   and falls within acceptence.
   New id is stored in prtOut.id.
   By default, this will use the momentum prtOut.p to make its
   determination, but you can set it to use the values stored in prt
   instead using SetPIDUseMC(true).
   */
  void Smear(const erhic::VirtualParticle&, ParticleMCS&);

  /**
   Dump the contents of the table to the screen.
   */
  void Speak();

  /**
   Clears existing table contents.
   */
  virtual void Clear(Option_t* = "");

  /**
   Print information about this device to standard output.
   */
  virtual void Print(Option_t* = "") const;

  TRandom3 Ran;
  TString PMatPath;
  std::vector<int> TrueIdent;
  std::vector<int> FalseIdent;
  std::vector<double> PMin;
  std::vector<double> PMax;
  // Indices are for momentum bin, true ID, and false ID
  // PMatrix[i][j][k] stores the probability for a particle of type
  // TrueIdent[j] in the momentum bin (PMin[i], PMax[i]) to be
  // identified as a particle of type FalseIdent[k].
  std::vector< std::vector<std::vector<double> > > PMatrix;
  // Range is the cumulative probability distribution of PMatrix.
  std::vector< std::vector<std::vector<double> > > Range;
  bool bUseMC;

  ClassDef(Smear::ParticleID, 1)
};

inline void ParticleID::SetPMatrixPath(TString str) {
  PMatPath = str;
  ReadP(PMatPath);
}

inline void ParticleID::SetPIDUseMC(bool d) {
  bUseMC = d;
}

inline TRandom3& ParticleID::GetRandomGenerator() {
  return Ran;
}

inline void ParticleID::SetRanSeed(int n) {
  Ran.SetSeed(n);
}

inline void ParticleID::GetAcceptanceFromDevice(const Device& dev) {
  Accept = dev.Accept;
}

inline ParticleID* ParticleID::Clone(const char*) const {
  // const char* argument comes from TObject::Clone(), usused here.
  return new ParticleID(*this);
}

}  // namespace Smear

#endif  // INCLUDE_EICSMEAR_SMEAR_PARTICLEID_H_
