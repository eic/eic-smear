//
// Forester.cxx
//
// Created by TB on 6/23/11.
// Copyright 2011 BNL. All rights reserved.
//

#include <iomanip>
#include <memory>
#include <stdexcept>

#include <TRefArray.h>

#include "File.h"
#include "Forester.h"
#include "functions.h" // For getFirstNonBlank()
#include "ParticleIdentifier.h"

ClassImp(erhic::Forester)

namespace erhic {
   
   
   Forester::Forester()
   : mQuit(false)
   , mVerbose(false)
   , mTree(NULL)
   , mEvent(NULL)
   , mMaxNEvents(0)
   , mInterval(1)
   , mInputName("default.txt")
   , mOutputName("default.root")
   , mTreeName("EICTree")
   , mBranchName("event") {
   }
   
   
   Forester::~Forester() {
      if(mFile) {
         delete mFile;
         mFile = NULL;
      } // if
      
      // We don't delete the output file as ROOT keeps track of all files.
      // We don't delete the mTree pointer because mRootFile
      // has ownership of it.
   }
   
   
   Long64_t
   Forester::Plant() {
      
      try {
         
         // Initialisation of the input and output files.
         OpenInput();
         SetupOutput();
         
         if(BeVerbose()) {
            std::cout << "\nProcessing " << GetInputFileName() << std::endl;
         } // if
         
         /** \todo Get rid of the static counter. Replace with a member
          that is reset every time Plant() is called.*/
         static int i(0);
         while(not MustQuit()) {
            ++i;
            if(BeVerbose() and i % mInterval == 0) {
               // Make the field just wide enough for the maximum
               // number of events.
               int width = static_cast<int>(::log10(GetMaxNEvents()) + 1);
               std::cout << "Processing event "<< std::setw(width) << i;
               if(GetMaxNEvents() > 0) {
                  std::cout << "/" << std::setw(width) << GetMaxNEvents();
               } // if
               std::cout << std::endl;
            } // if
            
            // Build the next event
            VirtualEvent<ParticleMC>* event = mFactory->Create();
            
            // Fill the tree
            if(event) {
               mTree->SetBranchAddress("event", &event);
               mTree->Fill();
               if(GetMaxNEvents() > 0 and i >= GetMaxNEvents()) {
                  SetMustQuit(true); // Hit max number of events, so quit
               } // if
               
               mStatus.ModifyEventCount(1);
               mStatus.ModifyParticleCount(event->GetNTracks());
               
               // We must ResetBranchAddress before deleting the event.
               mTree->ResetBranchAddress(mTree->GetBranch("event"));
               delete event;
            } // if
            else break;
         } // while

         Finish();
         
         return 0;
      } // try...
      
      catch(std::exception& e) {
         std::cerr << "Caught exception in Forester::Plant(): "
         << e.what() << std::endl;
         return -1;
      } // catch
   }
   
   
   bool
   Forester::OpenInput() {
      
      try {
         
         //	Open the input file for reading.
         mTextFile.open(GetInputFileName().c_str());
         
         // Throw a runtime_error if the file could not be opened.
         if(not mTextFile.good()) {
            std::string message("Unable to open file ");
            throw std::runtime_error(message.append(GetInputFileName()));
         }	//	if
         
         // Determine which Monte Carlo generator produced the file.
         mFile =
         erhic::FileFactory::GetInstance().GetFile(mTextFile);
         if(not mFile) {
            throw std::runtime_error(GetInputFileName() +
                                     " is not from a supported generator");
         } // for
         
         mFactory = mFile->CreateEventFactory(mTextFile);
         
         return true;
      } // try...
      
      // Pass the exception on to be dealt with higher up the food chain.
      catch(std::exception&) {
         throw;
      } // catch
   }
   
   
   bool
   Forester::SetupOutput() {
      
      try {
         
         // Open the ROOT file and check it opened OK
         mRootFile = new TFile(GetOutputFileName().c_str(), "RECREATE");
         if(not mRootFile->IsOpen()) {
            std::string message("Unable to open file ");
            throw std::runtime_error(message.append(GetOutputFileName()));
         } // if
         
         // Create the tree and check for errors
         mTree = new TTree(GetTreeName().c_str(), "my EIC tree");
         if(not mTree) {
            std::string message("Error allocating TTree ");
            throw std::runtime_error(message.append(GetTreeName()));
         }	//	if
         
         // Allocate memory for the branch buffer and
         // add the branch to the tree
         AllocateEvent();
         mTree->Branch(GetBranchName().c_str(), mEvent->ClassName(),
                       &mEvent, 32000, 99);
         
         // Auto-save every 500 MB
         mTree->SetAutoSave(500LL * 1024LL * 1024LL);
         
         // Align the input file at the start of the first event.
         FindFirstEvent();

         // Start timing after opening and creating files,
         // before looping over events
         mStatus.StartTimer();
         
         return true;
      } // try...
      
      catch(std::exception&) {
         throw;
      } // catch
   }
   
   
   void Forester::Finish() {
      
      if(BeVerbose()) {
         std::cout << "\nProcessed " << GetInputFileName() << std::endl;
         // If we received the quit signal, then we *should* have
         // completed an event.
         // If not, there will be particles in the event buffer:
         // the file must have terminated mid-event, which is a problem
         if(mEvent->GetNTracks() not_eq 0) {
            std::cerr <<
            "Warning: may have terminated mid-event - check input file"
            << std::endl;
         } // if
      } // if
      
      // Write the TTree to the file.
      mRootFile->cd();
      mTree->Write();
      mRootFile->ls();
      
      // Write the Forester itself to make it easier to reproduce the file
      // with the same settings.
      this->Write("forester");
      
      // Reset quit flag in case of further runs.
      SetMustQuit(false);
      
      // Stop timing the run.
      mStatus.StopTimer();
      
      if(BeVerbose()) {
         GetGetStatus().Print(std::cout); // Messages for the user
      } // if
      
      mRootFile->Close();
   }
   
   
   Int_t
   Forester::ProcessLine() {
      std::cout << "Processing line" << std::endl;
      int nTracksRead(0);
      
      if(AtEndOfEvent()) {
         // FinishEvent() will return <0 and set mQuit to true
         // if mMaxNEvents was reached
         std::cout << "At end of event" << std::endl;
         nTracksRead = FinishEvent();
         mStatus.ModifyEventCount(1);
         mStatus.ModifyParticleCount(nTracksRead);
      }	//	if...
      else if('0' == getFirstNonBlank(mLine)) {
         //	This is the first line of an event, so read the event data
         bool parsed = mEvent->Parse(mLine);
         if(BeVerbose() and not parsed) {
            std::cerr << "Warning: There was a problem reading event "
            << mTree->GetEntries() + 1 << " - check input file." << std::endl;
         } // if
      }	//	else if...
      else if('=' not_eq getFirstNonBlank(mLine)) {
         // This is a line of particle data, so add it to the current event
         AddParticle();
      }	//	else if
      
      return nTracksRead;
   }
   
   
   bool Forester::AddParticle() {
      
      std::auto_ptr<ParticleMC> particle(new ParticleMC(mLine));
      
      // That's a keeper!
      particle->SetEvent(mEvent);
      mEvent->AddLast(particle.release());
      
      return true;
   }
   
   
   Int_t
   Forester::FinishEvent() {
      Int_t number = mEvent->GetN();
      
      if(BeVerbose() and number % mInterval == 0) {
         // Make the field just wide enough for the maximum number of events.
         int width = static_cast<int>(::log10(GetMaxNEvents()) + 1);
         std::cout << "Processing event "<< std::setw(width) << number;
         if(GetMaxNEvents() > 0) {
            std::cout << "/" << std::setw(width) << GetMaxNEvents();
         } // if
         std::cout << std::endl;
      } // if
      
      mTree->Fill();
      
      Int_t nTracksRead = mEvent->GetNTracks();
      
      // Free memory before starting the next event.
      mEvent->Reset();
      
      // If an event limit was requested and we are at that limit, quit after
      // this event.
      const Long64_t nMax = GetMaxNEvents();
      if(nMax > 0 and number >= nMax) {
         SetMustQuit(true);
      } // if
      
      return nTracksRead;
   }
   
   
   char get1stNonBlank(std::istream& is) {
      char first = is.peek();
      while(not(first == ' ' or first == '\t')) {
         char c;
         is.get(c);
         first = is.peek();
      } // while
      return first;
   }
   
   
   bool
   Forester::AllocateEvent() {
      
      try {
         if(mEvent) {
            delete mEvent;
            mEvent = NULL;
         } // if
         mEvent = mFile->AllocateEvent();
         
         return mEvent;
      } // try...
      
      // Catch exceptions and pass up the food chain
      catch(std::exception&) {
         throw;
      } // catch
   }
   
   
   bool
   Forester::FindFirstEvent() {
      
      // Naughty kludge alert!
      // The first line was already read to determine the generator.
      // The header in the text files is six lines, so
      // read the remaining five lines of header.
      std::getline(mTextFile, mLine);
      std::getline(mTextFile, mLine);
      std::getline(mTextFile, mLine);
      std::getline(mTextFile, mLine);
      std::getline(mTextFile, mLine);
      
      return true;
   }
   
   void Forester::Print(std::ostream& os) const {
      os << "Input file: " << mInputName << std::endl;
      os << "Output file: " << mOutputName << std::endl;
      os << "Output tree: " << mTreeName << std::endl;
      os << "Output branch: " << mBranchName << std::endl;
      os << "Maximum number of events: " << mMaxNEvents << std::endl;
      if(mEvent) {
         os << "Event type: " << mEvent->ClassName() << std::endl;
      } // if
   }
   
} // namespace erhic
