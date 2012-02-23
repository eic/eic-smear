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

#include "EventMCFilterABC.h"
#include "ParticleMC.h"
#include "Pythia6.h"
#include "Pythia6EventBuilder.h"

namespace erhic {
   
   
   Pythia6::Pythia6(TFile* file,
                    int nEvents,
                    const std::string& treeName,
                    const std::string& branchName)
   : mPrintInterval(100)
   , mFile(file)
   , mTree(NULL)
   , mEvent(NULL)
   , mNEvents(nEvents)
   , mNGenerated(0)
   , mNTrials(0)
   , mFilter(NULL) {
      
      try {
         if(not file)
            throw std::runtime_error("No file provided");
         if(not file->IsWritable())
            throw std::runtime_error("File is not writable");
         
         file->cd();
         mTree = new TTree(treeName.c_str(),
                           "PYTHIA 6 events");
         
         mTree->Bronch(branchName.c_str(),
                       "EventPythia", &mEvent,
                       32000,
                       99);
      } // try
      catch(std::exception& e) {
         std::cout << "Caught exception in erhic::Pythia6::Pythia6() - " <<
         e.what() << std::endl;
         exit(1);
      }
   }
   
   Pythia6::~Pythia6() {
      // Don't delete the file as that was externally provided
      // Don't delete the tree as that is associated with the file
      // and will be deleted when the file closes.
   }
   
   bool Pythia6::Run() {
      
      Pythia6EventBuilder builder;
      
      TStopwatch timer;
      double lastTime(0.);
      
      while(mTree->GetEntries() < mNEvents) {
         
         const int initialNGenerated =
            TPythia6::Instance()->GetMSTI(5);
         const int initialNTrials =
            TPythia6::Instance()->GetPyint5()->NGEN[2][0];
         
         TPythia6::Instance()->GenerateEvent();
         mEvent = builder.Create();
         
         // Count the number of generations and trials for this event
         mNGenerated +=
            TPythia6::Instance()->GetMSTI(5) - initialNGenerated;
         mNTrials +=
            TPythia6::Instance()->GetPyint5()->NGEN[2][0] - initialNTrials;
         
         if(mFilter and not mFilter->Accept(*mEvent)) {
            delete mEvent;
            mEvent = NULL;
            continue;
         } // if
         
         // Event numbers count from 1, not zero, so add 1
         mEvent->SetN(mTree->GetEntries() + 1);
         mTree->SetBranchAddress("event", &mEvent);
         mTree->Fill();
         if((mTree->GetEntries() % mPrintInterval) == 0) {
            double time = timer.RealTime();
            
            std::cout << mTree->GetEntries() << " events in " <<
            time << " seconds (+" << time - lastTime << ")" << std::endl;
            lastTime = time;
            timer.Start(false);
         } // if
         mTree->ResetBranchAddress(mTree->GetBranch("event"));
         delete mEvent;
         mEvent = NULL;
      } // while
      
      // Write the tree
      mFile->cd();
      mTree->Write();
      mFile->Purge();
      
      // Write run information
      
      std::stringstream ss;
      std::string s;
      
      // Write the total cross section
      ss << TPythia6::Instance()->GetPARI(1) * 1000.; // * 1000 --> microbarn
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
