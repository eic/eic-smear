//
// Pythia6.cxx
//
// Created by TB on 1/18/12.
// Copyright 2012 BNL. All rights reserved.
//

#include <iostream>
#include <stdexcept>
#include <sstream>

#include <TFile.h>
#include <TObjString.h>
#include <TPythia6.h>
#include <TStopwatch.h>
#include <TTree.h>

#include "eicsmear/erhic/EventMCFilterABC.h"
#include "eicsmear/erhic/ParticleMC.h"
#include "eicsmear/erhic/Pythia6.h"
#include "eicsmear/erhic/EventFactory.h"

namespace erhic {
   
   // Need some special attention to get TTree::Branch to work with
   // these different types
   
   Pythia6::Pythia6(TFile* file,
                    VirtualEventFactory* factory,
                    int nEvents,
                    const std::string& treeName,
                    const std::string& branchName,
                    int printInterval)
   : mPrintInterval(printInterval)
   , mFile(file)
   , mTree(NULL)
   , mEvent(NULL)
   , mNEvents(nEvents)
   , mNGenerated(0)
   , mNTrials(0)
   , mFilter(NULL)
   , mFactory(factory) {
      try {
         if(not file) {
            throw std::runtime_error("No file provided");
         } // if
         if(not file->IsWritable()) {
            throw std::runtime_error("File is not writable");
         } // if
         file->cd();
         mTree = new TTree(treeName.c_str(), "PYTHIA 6 events");
         std::cout << factory->EventName().c_str() << std::endl;
//         mTree->Branch(branchName.c_str(), factory->EventName().c_str(),
//                       &mEvent, 32000, 99);
         mFactory->Branch(*mTree, branchName);
      } // try
      catch(std::exception& e) {
         std::cout << "Caught exception in erhic::Pythia6::Pythia6() - " <<
         e.what() << std::endl;
         exit(1);
      } // catch
   }
   
   Pythia6::~Pythia6() {
      // Don't delete the file as that was externally provided
      // Don't delete the tree as that is associated with the file
      // and will be deleted when the file closes.
   }
   
   bool Pythia6::Run() {
      TPythia6* pythia = TPythia6::Instance();
      
      TStopwatch timer;
      double lastTime(0.);
      
      while(mTree->GetEntries() < mNEvents) {
         
         const int initialNGenerated = pythia->GetMSTI(5);
         const int initialNTrials = pythia->GetPyint5()->NGEN[2][0];
         
         pythia->GenerateEvent();
//         mEvent = mFactory->Create();
         TBranch* branch = mTree->GetBranch("event");
         mFactory->Fill(*branch);
         
         // Count the number of generations and trials for this event
         mNGenerated += pythia->GetMSTI(5) - initialNGenerated;
         const int trials = pythia->GetPyint5()->NGEN[2][0] - initialNTrials;
//         mEvent->SetGenEvent(trials);
         mNTrials += trials;
         
         if(mFilter and not mFilter->Accept(*mEvent)) {
            if(mEvent) {
               delete mEvent;
            } // if
            mEvent = NULL;
            continue;
         } // if
         
         // Event numbers count from 1, not zero, so add 1
//         mEvent->SetN(mTree->GetEntries() + 1);
//         mTree->SetBranchAddress("event", &mEvent);
//         mTree->Fill();
         if((mTree->GetEntries() % mPrintInterval) == 0) {
            double time = timer.RealTime();
            
            std::cout << mTree->GetEntries() << " events in " <<
            time << " seconds (+" << time - lastTime << ")" << std::endl;
            lastTime = time;
            timer.Start(false);
         } // if
//         mTree->ResetBranchAddress(mTree->GetBranch("event"));
//         delete mEvent;
//         mEvent = NULL;
      } // while
      
      // Write the tree
      mFile->cd();
      mTree->Write();
      mFile->Purge();
      
      // Write run information
      
      std::stringstream ss;
      std::string s;
      
      // Write the total cross section
      ss << pythia->GetPARI(1) * 1000.; // * 1000 --> microbarn
      ss >> s;
      TObjString(s.c_str()).Write("crossSection");
      
      // Write the number of generated events
      ss.str("");
      ss.clear();
      ss << mNGenerated;
      ss >> s;
      TObjString(s.c_str()).Write("nEvents");
      
      // Write the number of trials required to make the events
      ss.str("");
      ss.clear();
      ss << mNTrials;
      ss >> s;
      TObjString(s.c_str()).Write("nTrials");
      
      return true;
   }
   
} // namespace erhic
