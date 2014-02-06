/**
 \file
 Here lies the function for smearing trees
 
 \author    Michael Savastio
 \date      2011-08-12
 \copyright 2011 Brookhaven National Lab
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <memory>

#include <TClass.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TChain.h>
#include <TString.h>
#include <TRandom2.h>
#include <TTree.h>
#include <TFile.h>
#include <TStopwatch.h>
#include <TH1D.h>

#include "eicsmear/erhic/VirtualParticle.h"
#include "eicsmear/smear/Detector.h"
#include "eicsmear/smear/EventDisFactory.h"
#include "eicsmear/smear/ParticleMCS.h"
#include "eicsmear/smear/Smear.h"

#ifdef WITH_PYTHIA6
#include "eicsmear/hadronic/EventSmear.h"
#endif

/**
 Smear nEvents events from the TTree named EICTree in the named input file,
 using the smearing definitions in the Detector.
 Write the resulting Smeared TTree to a file named outFileName.
 If nEvents <= 0 smear all events in the tree.
 Returns 0 upon success, 1 upon failure.
 */
int SmearTree(const Smear::Detector& detector, const TString& inFileName,
              const TString& outFileName, Long64_t nEvents) {
  // Open the input file and get the Monte Carlo tree from it.
  // Complain and quit if we don't find the file or the tree.
  TFile inFile(inFileName, "READ");
  if (!inFile.IsOpen()) {
    std::cerr << "Unable to open " << inFileName << std::endl;
  }  // if
  TTree* mcTree(NULL);
  inFile.GetObject("EICTree", mcTree);
  if (!mcTree) {
    std::cerr << "Unable to find EICTree in " << inFileName << std::endl;
    return 1;
  }  // if
  std::auto_ptr<erhic::VirtualEventFactory> builder;
  // Need to determine the type of object in the tree to choose
  // the correct smeared event builder.
  TClass branchClass(mcTree->GetBranch("event")->GetClassName());
  if (branchClass.InheritsFrom("erhic::EventDis")) {
    builder.reset(new Smear::EventDisFactory(detector,
                                             *(mcTree->GetBranch("event"))));
#ifdef WITH_PYTHIA6
  } else if (branchClass.InheritsFrom("erhic::hadronic::EventMC")) {
    builder.reset(new Smear::HadronicEventBuilder(detector,
                                             *(mcTree->GetBranch("event"))));
#endif
  } else {
    std::cerr << branchClass.GetName() << " is not supported for smearing" <<
    std::endl;
  }  // if
  // Open the output file.
  // Complain and quit if something goes wrong.
  TString outName(outFileName);
  if (outName.IsNull()) {
    outName = TString(inFileName).ReplaceAll(".root", ".smear.root");
  }  // if
  TFile outFile(outName, "RECREATE");
  if (!outFile.IsOpen()) {
    std::cerr << "Unable to create " << outName << std::endl;
    return 1;
  }  // if
  TTree smearedTree("Smeared", "A tree of smeared Monte Carlo events");
  TBranch* eventbranch = builder->Branch(smearedTree, "eventS");
  if (mcTree->GetEntries() < nEvents || nEvents < 1) {
    nEvents = mcTree->GetEntries();
  }  // if
  std::cout <<
  "/-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-/"
  << std::endl;
  std::cout <<
  "/  Commencing Smearing of " << nEvents << " events."
  << std::endl;
  std::cout <<
  "/-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-/"
  << std::endl;
  for (Long64_t i(0); i < nEvents; i++) {
    if (i % 10000 == 0 && i != 0) {
      std::cout << "Processing event " << i << std::endl;
    }  // if
    mcTree->GetEntry(i);
    builder->Fill(*eventbranch);
  }  // for
  smearedTree.Write();
  detector.Write("detector");
  outFile.Purge();
  std::cout <<
  "|~~~~~~~~~~~~~~~~~~ Completed Successfully ~~~~~~~~~~~~~~~~~~~|"
  << std::endl;
  return 0;
}
