/**
 \file
 Declaration of class erhic::Forester.
 
 \author    Thomas Burton
 \date      2011-06-23
 \copyright 2011 Brookhaven National Lab
*/

#ifndef INCLUDE_EICSMEAR_ERHIC_FORESTER_H_
#define INCLUDE_EICSMEAR_ERHIC_FORESTER_H_

// C(++) headers
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

// ROOT headers
#include <Rtypes.h>
#include <TFile.h>
#include <TStopwatch.h>
#include <TTree.h>

// Other headers
#include "eicsmear/erhic/EventMC.h"

namespace erhic {

class FileType;
class VirtualEventFactory;

/**
 Manages the creation of trees from plain-text Monte Carlo files.
 Bad pun, I know, but the ROOT guys started it.
 */
class Forester : public TObject {
 public:
  /**
   Default constructor.
   */
  Forester();

  /**
   Destructor.
   */
  virtual ~Forester();

  /**
   Processes a text file into a ROOT file.
   Returns the number of events processed.
   */
  Long64_t Plant();

  /**
   Sets the name of the input text file containing Monte Carlo data.
   */
  void SetInputFileName(const std::string&);

  /**
   Sets the name of the ROOT tree file to create.
   */
  void SetOutputFileName(const std::string&);

  /**
   Sets the name of the TTree to write to the file named by 
   SetOutputFileName().
   */
  void SetTreeName(const std::string& = "EICTree");

  /**
   Sets the name of the TBranch containing event objects.
   This is the only branch written to the TTree named by SetTreeName().
   */
  void SetBranchName(const std::string& = "event");

  /**
   Returns the name of the input text file containing Monte Carlo data.
   */
  std::string GetInputFileName() const;

  /**
   Returns the name of the ROOT tree file to create.
   */
  std::string GetOutputFileName() const;

  /**
   Returns the name of the TTree to write to the file named by 
   SetOutputFileName().
   */
  std::string GetTreeName() const;

  /**
   Returns the name of the TBranch containing event objects.
   */
  std::string GetBranchName() const;

  /**
   Sets the maximum number of events to process. Processing will terminate
   when this number of events or the end of the input file is reached, 
   whichever occurs first. A value <= 0 indicates to process all events in the
   input file (this is the default).
   */
  void SetMaxNEvents(Long64_t = 0);

  /**
   Returns the maximum number of events to process.
   */
  Long64_t GetMaxNEvents() const;

  /**
   Sets the event count interval at which to print a status message.
   A value <= 0 suppresses messages.
   */
  void SetMessageInterval(Long64_t = 10000);

  /**
   Prints the current configuration to the requested output stream.
   */
  void Print(std::ostream& stream) const;

  /**
   \overload

   To avoid hiding TObject::Print()
   */
  void Print(Option_t* = "not used") const;

  /**
   If set to true, prints messages during running.
   If set to false, runs silently except in the case of critical errors.
   */
  void SetBeVerbose(bool = false);

  /**
   Returns the verbosity i.e. whether to print status messages.
   */
  bool BeVerbose() const;

  /**
   Returns the file type information for the last processed file.
   Returns NULL if no input has been processed.
   Do not delete the returned object.
   */
  const erhic::FileType* GetFileType() const;

 protected:
  /**
   Stores summary information about the last call to Forester::Plant().
   */
  class Status {
   public:
    Status()
    : mNEvents(0)
    , mNParticles(0) {
      // Initialise the start and end time to the creation time and reset
      // the timer to ensure it is at zero.
      std::time(&mStartTime);
      mEndTime = mStartTime;
      mTimer.Reset();
    }

    virtual ~Status() { }

    virtual std::ostream& Print(std::ostream& os = std::cout) const {
      // Put start and end times in different os <<... otherwise I get
      // the same time for each...
      os << "Began on " << std::ctime(&mStartTime);
      os << "Ended on " << std::ctime(&mEndTime);
      os << "Processed " << mNEvents << " events containing "
      << mNParticles << " particles in "
      << mTimer.RealTime() << " seconds "
      << '(' << mTimer.RealTime()/mNEvents <<" sec/event)" << std::endl;
      return os;
    }

  protected:
    virtual void StartTimer() {
      std::time(&mStartTime);
      mTimer.Start();
    }

    virtual void StopTimer() {
      std::time(&mEndTime);
      mTimer.Stop();
    }

    virtual void ModifyEventCount(Long64_t count) {
      mNEvents += count;
    }

    virtual void ModifyParticleCount(Long64_t count) {
      mNParticles += count;
    }

    time_t mStartTime;
    time_t mEndTime;
    Long64_t mNEvents;
    Long64_t mNParticles;

    // The TStopwatch is mutable as "GetRealTime()" is non-const.
    mutable TStopwatch mTimer;

    friend class Forester;

    ClassDef(Status, 1)
  };

  // End of class Status

  /**
   Prints a summary of the last call to Plant()
   to the requested output stream.
   */
  const Status& GetGetStatus() const {
    return mStatus;
  }

  /**
   Opens the input file and checks that it was produced by a
   supported Monte Carlo generator.
   Returns true upon success or false upon and I/O error or
   if the Monte Carlo generator is unsupported or cannot be determined.
   */
  bool OpenInput();

  /**
   Opens the output ROOT file and creates the TTree ready for filling.
   */
  bool SetupOutput();

  /**
   Writes output and takes end-of-file actions.
   */
  void Finish();

  /**
   Allocate an event buffer for the TTree based on the generator type.
   Returns false if the generator has not yet been successfully
   determined.
   */
  bool AllocateEvent();

  /**
   Aligns the input text file on the first line of the first event.
   After a successful call, mLine stores the first line of the event
   and true is returned.
   If unsuccessful, mLine is blank and false is returned.
   */
  bool FindFirstEvent();

  /** Prints the status of the current Plant() call to the standard output. */
  void PrintStatus() const;

  /**
   Prints the quit flag status. A return value of true indicates
   that the input file has ended or the mMaxNEvents has been reached.
   Processing the file will
   quit after the end of the next call to FinishEvent().
   */
  bool MustQuit() const;

  /** Set the quit flag. */
  void SetMustQuit(bool quit);

  // Member variables.
  // Those with comment //!< will be treated as transient members by ROOT
  // and won't be written to a file.

  Bool_t mQuit;  ///< Quit status. Set to true once EoF or max events reached
  Bool_t mVerbose;  ///< Verbosity flag
  TTree* mTree;  //!< Output TTree, owned by mRootFile
  VirtualEvent* mEvent;  //!< Stores event branch address
  const erhic::FileType* mFile;  //!< File type information
  TFile* mRootFile;  //!< Pointer to output ROOT file
  Long64_t mMaxNEvents;  ///< Maximum number of events to process
  Long64_t mInterval;  ///< Event interval between printing status messages
  std::ifstream* mTextFile;  //!< Input text file
  std::string mInputName;  ///< Name of the input text file
  std::string mOutputName;  ///< Name of the output ROOT file
  std::string mTreeName;  ///< Name of the output TTree
  std::string mBranchName;  ///< Name of the event TBranch
  std::string mLine;  ///< Stores the latest text line read from the input file
  Status mStatus;  ///< Forester status information
  VirtualEventFactory* mFactory;  //!< Pointer to the event-builder object

  ClassDef(Forester, 1)
};

inline void Forester::SetInputFileName(const std::string& name) {
  mInputName = name;
}

inline void Forester::SetOutputFileName(const std::string& name) {
  mOutputName = name;
}

inline void Forester::SetTreeName(const std::string& name) {
  mTreeName = name;
}

inline void Forester::SetBranchName(const std::string& name) {
  mBranchName = name;
}

inline std::string Forester::GetInputFileName() const {
  return mInputName;
}

inline std::string Forester::GetOutputFileName() const {
  return mOutputName;
}

inline std::string Forester::GetTreeName() const {
  return mTreeName;
}

inline std::string Forester::GetBranchName() const {
  return mBranchName;
}

inline void Forester::SetMaxNEvents(Long64_t number) {
  mMaxNEvents = number;
}

inline Long64_t Forester::GetMaxNEvents() const {
  return mMaxNEvents;
}

inline void Forester::SetMessageInterval(Long64_t number) {
  mInterval = number;
}

inline bool Forester::MustQuit() const {
  return mQuit;
}

inline void Forester::SetMustQuit(bool flag) {
  mQuit = flag;
}

inline void Forester::SetBeVerbose(bool flag) {
  mVerbose = flag;
}

inline bool Forester::BeVerbose() const {
  return mVerbose;
}

inline const FileType* Forester::GetFileType() const {
  return mFile;
}

}  // namespace erhic

#endif  // INCLUDE_EICSMEAR_ERHIC_FORESTER_H_
