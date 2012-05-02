/**
 SmearTree.cxx
 
 \file
 Here lies the function for smearing trees
 
 \author Michael Savastio
 \date 08/12/2011
 \copyright 2011 BNL. All rights reserved.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

#include <TROOT.h>
#include <TSystem.h>
#include <TChain.h>
#include <TString.h>
#include <TRandom2.h>
#include <TTree.h>
#include <TFile.h>
#include <TStopwatch.h>
#include <TH1D.h>

#include "EventBase.h"
#include "Detector.h"
#include "ParticleMC.h"
#include "ParticleMCS.h"
#include "Smear.h"

// Return a ParticleMCS with the same properties as the Monte Carlo particle
Smear::ParticleMCS* mcToSmear(const erhic::ParticleMC& mc) {
   return new Smear::ParticleMCS(mc.Get4Vector(), mc.Id(), mc.GetStatus());
}

/**
 Smear nEvents events from the TTree named EICTree in the named input file,
 using the smearing definitions in the Detector.
 Write the resulting Smeared TTree to a file named outFileName.
 If nEvents <= 0 smear all events in the tree.
 Returns 0 upon success, 1 upon failure.
 */
int SmearTree(Smear::Detector detector, TString inFileName,
              TString outFileName = "", Long64_t nEvents= -1) {
   // Open the input file and get the Monte Carlo tree from it.
   // Complain and quit if we don't find the file or the tree.
   TFile inFile(inFileName, "READ");
   if(not inFile.IsOpen()) {
      std::cerr << "Unable to open " << inFileName << std::endl;
   } // if
   TTree* mcTree(NULL);
   inFile.GetObject("EICTree", mcTree);
   if(not mcTree) {
      std::cerr << "Unable to find EICTree in " << inFileName << std::endl;
      return 1;
   } // if
   erhic::EventMC* event(NULL);
   mcTree->SetBranchAddress("event", &event);
   // Open the output file.
   // Complain and quit if something goes wrong.
   if(outFileName.IsNull()) {
      outFileName = inFileName.ReplaceAll(".root", ".smear.root");
   } // if
   TFile outFile(outFileName, "RECREATE");
   if(not outFile.IsOpen()) {
      std::cerr << "Unable to create " << outFileName << std::endl;
      return 1;
   } // if
   TTree smearedTree("Smeared", "A tree of smeared Monte Carlo events");

   Smear::Event* eventS(NULL);

   smearedTree.Branch("eventS",&eventS);

   if(mcTree->GetEntries() < nEvents or nEvents < 1) {
      nEvents = mcTree->GetEntries();
   } // if
   
   std::cout <<
   "/-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-/"
   << std::endl;
   std::cout <<
   "/  Commencing Smearing of " << nEvents << " events."
   << std::endl;
   std::cout <<
   "/-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-/"
   << std::endl;
   
   for(Long64_t i(0); i < nEvents; i++) {
      if(i % 10000 == 0 and i not_eq 0) {
         std::cout << "Processing event " << i << std::endl;
      } // if
      mcTree->GetEntry(i);
      // We need the scattered lepton to compute event kinematics.
      // There isn't an easy way to determine which particle is the scattered
      // lepton from the smeared record after the fact, so we record its
      // index while we are looping through the particle list.
      ParticleIdentifier pid(event->BeamLepton()->Id());
      for(unsigned j(0); j < event->GetNTracks(); j++) {
         const erhic::ParticleMC* ptr = event->GetTrack(j);
         if(not ptr) {
            continue;
         } // if
         // It's convenient to keep the initial beams, unsmeared, in the
         // smeared event record, so copy their properties exactly
         if(event->BeamLepton() == ptr or event->BeamHadron() == ptr) {
            eventS->particles.push_back(mcToSmear(*ptr));
         } // if
         else {
            eventS->particles.push_back(detector.Smear(*ptr));
            // If this is the scattered lepton, record the index.
            // Check that the scattered lepton hasn't been set yet so we
            // don't replace it with a subsequent match.
            if(pid.isScatteredLepton(*ptr) and not eventS->ScatteredLepton()) {
               eventS->mScatteredIndex = ptr->GetIndex() - 1;
            } // if
         } // else
      } // for
      // Fill the event-wise kinematic variables.
      detector.FillEventKinematics(*event, eventS);
      smearedTree.Fill();
      eventS->Reset();
   } // for
   smearedTree.Write();
   detector.Write("detector");
   outFile.Purge();
   std::cout <<
   "|~~~~~~~~~~~~~~~~~~ Completed Successfully ~~~~~~~~~~~~~~~~~~~|"
   << std::endl;
   return 0;
}
