//
// Forester.cxx
//
// Created by TB on 6/23/11.
// Copyright 2011 BNL. All rights reserved.
//

#include <iomanip>
#include <memory>
#include <stdexcept>

//#include <TLorentzRotation.h>
//#include <TLorentzVector.h>
//#include <TROOT.h>

/*#include "EventPythia.h"
#include "EventPepsi.h"
#include "EventDjangoh.h"
#include "EventMilou.h"
#include "EventRapgap.h"
 */
#include "File.h"
#include "Forester.h"
#include "functions.h" // For getFirstNonBlank()
#include "ParticleIdentifier.h"
#include <TRefArray.h>
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
      // We don't delete the mTree pointer because mRootFile has ownership of it.
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
               mStatus.ModifyParticleCount(event->NTracks());
               
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
         erhic::monte_carlo::FileFactory::GetInstance().GetFile(mTextFile);
         if(not mFile) {
            throw std::runtime_error(GetInputFileName() +
                                     " is not from a supported generator");
         } // for
         
         mFactory = mFile->CreateEventFactory(mTextFile);
         std::cout << "Factory at " << mFactory << std::endl;
         
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
         
         // Open the ROOT file.
         //      mRootFile.reset(new TFile(GetOutputFileName().c_str(), "RECREATE"));
         mRootFile = new TFile(GetOutputFileName().c_str(), "RECREATE");
         if(not mRootFile->IsOpen()) {
            std::string message("Unable to open file ");
            throw std::runtime_error(message.append(GetOutputFileName()));
         } // if
         
         mTree = new TTree(GetTreeName().c_str(), "my EIC tree");
         mRootFile->ls();
         
         if(not mTree) {
            std::string message("Error allocating TTree ");
            throw std::runtime_error(message.append(GetTreeName()));
         }	//	if
         
         // Allocate memory for the branch buffer and add the branch to the tree
         AllocateEvent();
         mTree->Bronch(GetBranchName().c_str(),
                       mEvent->ClassName(), 
                       &mEvent,
                       32000,
                       99);
//         mTree->BranchRef(); // For TRefArray in ParticleMC
         
         mTree->SetAutoSave(500LL * 1024LL * 1024LL); // Auto-save every 500 MB
         
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
         // If we received the quit signal, then we should have completed an event.
         // If not, there will be particles in the event buffer - the file must
         // have terminated mid-event.
         if(mEvent->NTracks() not_eq 0) {
            std::cerr <<
            "Warning: processing may have finished mid-event - check input file"
            << std::endl;
         } // if
      } // if
      
      // Write the TTree to the file.
//      mRootFile->Purge(); // Remove all but the latest tree version
//      mRootFile->Write();
      mRootFile->cd();
      mRootFile->ls();
      mTree->Write();
      std::cout << mTree->GetEntries() << std::endl;
      mRootFile->ls();
      
      // Write the Forester itself to make it easier to reproduce the file
      // with the same settings.
      this->Write("forester");
      
      // Reset quit flag for further runs.
      SetMustQuit(false);
      
      // Stop timing the run.
      mStatus.StopTimer();
      
      if(BeVerbose()) {
         GetStatus().Print(std::cout); // Messages for the user
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
      
//      ParticleMC particle(mLine);
      std::auto_ptr<ParticleMC> particle(new ParticleMC(mLine));
//      std::cout << "AddParticle " << std::endl;
//      particle.Dump();
      
//      if(particle.GetStatus() > 21) {
//         return false;
//      } // if
      
//      static ParticleIdentifier identifier;
      
      // The first particle in the event should be defined to be the beam lepton.
      // Read this particle and use its species to set the beam species type for
      // this event.
      // Set this before checking SkipParticle() as its result depends on the
      // lepton beam species.
      //      if(mEvent->particles.empty()) {
      //         identifier.SetLeptonBeamPdgCode(particle.id);
      //      } // if
      
      // I want to implement each Particle containing the parent particle's
      // PDG code. This is made easier by having no "jumps" in the particle
      // list, because I can then just look up the parent via its index.
      //      if(identifier.SkipParticle(particle)) {
      //         return false;
      //      } // if
      
      // Calculate properties based on those we read from the Monte Carlo record
#if 0      
//      particle.ComputeDerivedQuantities();
      
      // That's a keeper!
      mEvent->particles.push_back(particle);
      mEvent->particles.back().event = mEvent;
#endif
//      particle->ComputeDerivedQuantities();
      
      // That's a keeper!
      particle->SetEvent(mEvent);
      mEvent->AddLast(particle.release());
      
      return true;
   }
   
   Int_t
   Forester::FinishEvent() {
      Int_t number = mEvent->N();
      
      // DisplayProgress();
      // We have reached the end of an event.
      // Calculate any remaining quantities that require information about the
      // whole event, then fill the tree.
      if(BeVerbose() and number % mInterval == 0) {
         // Make the field just wide enough for the maximum number of events.
         int width = static_cast<int>(::log10(GetMaxNEvents()) + 1);
         std::cout << "Processing event "<< std::setw(width) << number;
         if(GetMaxNEvents() > 0) {
            std::cout << "/" << std::setw(width) << GetMaxNEvents();
         } // if
         std::cout << std::endl;
      } // if
#if 0
      // Find the beam particles/virtual photon and store in beams.
      // Calculate event and particle quantities dependent on the beam properties.
      ParticleIdentifier::IdentifyBeams(*mEvent, beams);
//      mEvent->Compute(beams);
      
      if(not mEvent->Compute()) {
         std::cerr << "Event " << number <<
         ": computation of event kinematics failed" << std::endl;
      } // if
#endif
#if 0
      // Set parent daughter relations:
      const unsigned nParticles = mEvent->particles.size();
      for(unsigned i(0); i < nParticles; ++i) {
//         std::cout << "Particle " << i << std::endl;
         UShort_t idaughter1 = mEvent->particles.at(i).daughter;
         UShort_t idaughterN = mEvent->particles.at(i).ldaughter;
//         std::cout << "daughter1 = " << idaughter1 << " daughterN = " << idaughterN << std::endl;
         if(idaughter1 == 0) continue; // No daughters
         if(idaughterN == 0) idaughterN = idaughter1; // 1 daughter
         
         // PYTHIA output gives daughters in range [1,N] while the particles
         // array is indexed [0, N-1], so subtract 1 from the daughter indices
         // to get the array indices (there are no gaps in the array).
         --idaughter1;
         --idaughterN;
         
         for(unsigned j(idaughter1); j <= idaughterN; ++j) {
            std::cout << "\tDaughter " << j << " of particle " << i << std::endl;
//            mEvent->particles.at(idaughter1).Print();
//            mEvent->particles.at(i).mChildrenAdd(&(mEvent->particles.at(j)));
            mEvent->particles.at(i).mChildren.Add(&(mEvent->particles.at(j)));
//            mEvent->particles.at(i).children.push_back(&(mEvent->particles.at(j)));
         } // for
         std::cout << "Particle " << i << " has " << mEvent->particles.at(i).mChildren.GetLast() + 1 << " daughters" << std::endl;
      } // for
#endif
      
      mTree->Fill();
      
      Int_t nTracksRead = mEvent->NTracks();
      
      // Free memory before starting the next event.
      mEvent->Reset();
      
//      std::cout << ParticleMC::smNInstances << " particles" << std::endl;
      
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
         
#if 0
         // We now need to 'wind' the file forward until we reach the start
         // of the first event (indicated by a '0' as the first non-blank
         // character of a line).
         char firstNonBlank = get1stNonBlank(mTextFile);
         while((firstNonBlank not_eq '0')) {
            std::string line;
            std::getline(mTextFile, line);
         } // while
#endif
         
         return mEvent;
      } // try...
      
      // Catch exceptions and pass up the food chain
      catch(std::exception&) {
         throw;
      } // catch
   }
   
   bool
   Forester::FindFirstEvent() {
      // Counter for the number of lines read without finding the first event.
//      int failCount(0);
      
      // Set this flag to true if we find the start of the first event.
      bool found(false);
      
      // Read lines until the fail counter exceeds the maximum or the input stream
      // is no longer good for I/O operations.
#if 0
      while(failCount < 10 and
            std::getline(mTextFile, mLine).good()) {
         
         // If we found the first event, set the flag to true and break out of the
         // loop.
         if('0' == getFirstNonBlank(mLine)) {
            found = true;
            break;
         } // if
         
         // This line isn't the first event. Increment the 'fail' counter.
         ++failCount;
      } // while
#endif
      // The first line was already read to determine the generator.
      // Read the remaining 5 lines of header.
      std::getline(mTextFile, mLine);
//      std::cout << mLine << std::endl;
      std::getline(mTextFile, mLine);
//      std::cout << mLine << std::endl;
      std::getline(mTextFile, mLine);
//      std::cout << mLine << std::endl;
      std::getline(mTextFile, mLine);
//      std::cout << mLine << std::endl;
      std::getline(mTextFile, mLine);
//      std::cout << mLine << std::endl;
      
      if(not found) mLine.clear();
      
      return found;
   }
   
} // namespace erhic
