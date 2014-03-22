/**
 \file
 Implementation of class erhic::Forester.
 
 \author    Thomas Burton
 \date      2011-06-23
 \copyright 2011 Brookhaven National Lab
 */

#include "eicsmear/erhic/Forester.h"

#include <iomanip>
#include <memory>
#include <stdexcept>
#include <string>

#include <TRefArray.h>

#include "eicsmear/erhic/EventFactory.h"
#include "eicsmear/erhic/File.h"
#include "eicsmear/erhic/ParticleIdentifier.h"

namespace erhic {

Forester::Forester()
: mQuit(false)
, mVerbose(false)
, mTree(NULL)
, mEvent(NULL)
, mFile(NULL)
, mRootFile(NULL)
, mMaxNEvents(0)
, mInterval(1)
, mTextFile(NULL)
, mInputName("default.txt")
, mOutputName("default.root")
, mTreeName("EICTree")
, mBranchName("event")
, mFactory(NULL) {
}

Forester::~Forester() {
  if (mFile) {
    delete mFile;
    mFile = NULL;
  }  // if
  if (mEvent) {
    delete mEvent;
    mEvent = NULL;
  }  // if
  if (mFactory) {
    delete mFactory;
    mFactory = NULL;
  }  // if
  if (mRootFile) {
    delete mRootFile;
    mRootFile = NULL;
  }  // if
  if (mTextFile) {
    delete mTextFile;
    mTextFile = NULL;
  }  // if
  // We don't delete the mTree pointer because mRootFile
  // has ownership of it.
}

Long64_t Forester::Plant() {
  try {
    // Initialisation of the input and output files.
    OpenInput();
    SetupOutput();
    if (BeVerbose()) {
      std::cout << "\nProcessing " << GetInputFileName() << std::endl;
    }  // if
    /**
     \todo Get rid of the static counter. Replace with a member
     that is reset every time Plant() is called.
     */
    static int i(0);
    while (!MustQuit()) {
      ++i;
      if (BeVerbose() && i % mInterval == 0) {
        // Make the field just wide enough for the maximum
        // number of events.
        int width = static_cast<int>(::log10(GetMaxNEvents()) + 1);
        std::cout << "Processing event "<< std::setw(width) << i;
        if (GetMaxNEvents() > 0) {
          std::cout << "/" << std::setw(width) << GetMaxNEvents();
        }  // if
        std::cout << std::endl;
      }  // if
      // Build the next event
      if (mEvent) {
        delete mEvent;
        mEvent = NULL;
      }  // if
      // Catch exceptions from event builder here so we don't break
      // out of the whole tree building loop for a single bad event.
      try {
        mEvent = mFactory->Create();
        // Fill the tree
        if (mEvent) {
          mTree->Fill();
          if (GetMaxNEvents() > 0 && i >= GetMaxNEvents()) {
            SetMustQuit(true);  // Hit max number of events, so quit
          }  // if
          mStatus.ModifyEventCount(1);
          mStatus.ModifyParticleCount(mEvent->GetNTracks());
          // We must ResetBranchAddress before deleting the event.
        } else {
          break;
        }  // if
      }  // try
      catch(std::exception& e) {
        std::cerr << "Caught exception in Forester::Plant(): "
        << e.what() << std::endl;
        std::cerr << "Event will be skipped..." << std::endl;
      }  // catch
    }  // while
    Finish();
    return 0;
  }  // try
  catch(std::exception& e) {
    std::cerr << "Caught exception in Forester::Plant(): "
    << e.what() << std::endl;
    return -1;
  }  // catch
}

bool Forester::OpenInput() {
  try {
    // Open the input file for reading.
    if (!mTextFile) {
      mTextFile = new std::ifstream;
    }  // if
    mTextFile->open(GetInputFileName().c_str());
    // Throw a runtime_error if the file could not be opened.
    if (!mTextFile->good()) {
      std::string message("Unable to open file ");
      throw std::runtime_error(message.append(GetInputFileName()));
    }  // if
    // Determine which Monte Carlo generator produced the file.
    mFile =
    erhic::FileFactory::GetInstance().GetFile(*mTextFile);
    if (!mFile) {
      throw std::runtime_error(GetInputFileName() +
                               " is not from a supported generator");
    }  // for
    mFactory = mFile->CreateEventFactory(*mTextFile);
    return true;
  }  // try...
  // Pass the exception on to be dealt with higher up the food chain.
  catch(std::exception&) {
    throw;
  }  // catch
}

bool Forester::SetupOutput() {
  try {
    // Open the ROOT file and check it opened OK
    mRootFile = new TFile(GetOutputFileName().c_str(), "RECREATE");
    if (!mRootFile->IsOpen()) {
      std::string message("Unable to open file ");
      throw std::runtime_error(message.append(GetOutputFileName()));
    }  // if
    // Create the tree and check for errors
    mTree = new TTree(GetTreeName().c_str(), "my EIC tree");
    if (!mTree) {
      std::string message("Error allocating TTree ");
      throw std::runtime_error(message.append(GetTreeName()));
    }  // if
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
  }  // try...
  catch(std::exception&) {
    throw;
  }  // catch
}

void Forester::Finish() {
  if (BeVerbose()) {
    std::cout << "\nProcessed " << GetInputFileName() << std::endl;
  }  // if
  // Write the TTree to the file.
  mRootFile = mTree->GetCurrentFile();
  mRootFile->cd();
  mTree->Write();
  mRootFile->ls();
  // Write the Forester itself to make it easier to reproduce the file
  // with the same settings.
  Write("forester");
  // Reset quit flag in case of further runs.
  SetMustQuit(false);
  // Stop timing the run.
  mStatus.StopTimer();
  if (BeVerbose()) {
    GetGetStatus().Print(std::cout);  // Messages for the user
  }  // if
  mRootFile->Close();
}

bool Forester::AllocateEvent() {
  try {
    if (mEvent) {
      delete mEvent;
      mEvent = NULL;
    }  // if
    mEvent = mFile->AllocateEvent();
    return mEvent;
  }  // try...
  // Catch exceptions and pass up the food chain
  catch(std::exception&) {
    throw;
  }  // catch
}

bool Forester::FindFirstEvent() {
  // Naughty kludge alert!
  // The first line was already read to determine the generator.
  // The header in the text files is six lines, so
  // read the remaining five lines of header.
  std::getline(*mTextFile, mLine);
  std::getline(*mTextFile, mLine);
  std::getline(*mTextFile, mLine);
  std::getline(*mTextFile, mLine);
  std::getline(*mTextFile, mLine);
  return true;
}

void Forester::Print(std::ostream& os) const {
  os << "Input file: " << mInputName << std::endl;
  os << "Output file: " << mOutputName << std::endl;
  os << "Output tree: " << mTreeName << std::endl;
  os << "Output branch: " << mBranchName << std::endl;
  os << "Maximum number of events: " << mMaxNEvents << std::endl;
  if (mEvent) {
    os << "Event type: " << mEvent->ClassName() << std::endl;
  }  // if
}

void Forester::Print(Option_t* /* not used */) const {
  Print(std::cout);
}
}  // namespace erhic
