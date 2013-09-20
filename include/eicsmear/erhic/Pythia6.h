/**
 \file
 Declaration of class erhic::Pythia6.

 \author    Thomas Burton
 \date      2012-01-18
 \copyright 2012 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_ERHIC_PYTHIA6_H_
#define INCLUDE_EICSMEAR_ERHIC_PYTHIA6_H_

#include <memory>
#include <string>

#include <Rtypes.h>  // For ClassDef macro

class TFile;
class TTree;

namespace erhic {

class EventPythia;
class EventMCFilterABC;
class VirtualEvent;
class VirtualEventFactory;

/**
 Runs PYTHIA 6 and builds an output tree file in the ROOT environment.
 The output tree is populated with EventPythia objects.
 PYTHIA options (including initialisation) should be provided via the
 ROOT interface class TPythia6 before running the tree-building.
 */
class Pythia6 {
 public:
  /**
   Constructor.
   
   Associates the output of this PYTHIA run with a file.
   Generate events using the provided event factory.
   Both the file and the factory should be dynamically allocated.
   Pythia6 takes ownership and deletes them.
   Define the number of events to produce, and optionally
   provide a tree and branch name.
   */
  Pythia6(TFile* file,
          VirtualEventFactory*,
          int nEvents,
          const std::string& treeName = "EICTree",
          const std::string& branchName = "event",
          int printInterval = 1000);

  /**
   Destructor
   */
  virtual ~Pythia6();

  /**
   Runs PYTHIA and writes output
   \todo Implement selection of correct event factory (either
   erhic::Pythia6EventBuilder or hadronic::Pythia6EventFactory)
   depending on requested beam types.
   */
  virtual bool Run();

 protected:
  int mPrintInterval;
  TFile* mFile;           //!< Pointer to the output file
  TTree* mTree;           ///< Pointer to the generated tree
  const int mNEvents;     ///< Number of events to produce
  int mNGenerated;        ///< Number of events passing PYTHIA selection
  int mNTrials;           ///< Number of events thrown by PYTHIA
  std::auto_ptr<VirtualEventFactory> mFactory;  //!< Event factory.

  ClassDef(erhic::Pythia6, 1)
};

}  // namespace erhic

#endif  // INCLUDE_EICSMEAR_ERHIC_PYTHIA6_H_
