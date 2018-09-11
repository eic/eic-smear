/**
 \file
 Declarations of file information and log-reading classes.
 
 \author    Thomas Burton
 \date      2011-07-29
 \copyright 2011 Brookhaven National Lab
 */

#ifndef INCLUDE_EICSMEAR_ERHIC_FILE_H_
#define INCLUDE_EICSMEAR_ERHIC_FILE_H_

#include <iostream>
#include <map>
#include <string>

#include <TObject.h>
#include <TObjString.h>
#include <TString.h>

#include "eicsmear/erhic/EventBase.h"
#include "eicsmear/erhic/EventMilou.h"
#include "eicsmear/erhic/EventPythia.h"
#include "eicsmear/erhic/EventFactory.h"

namespace erhic {

/**
 Base class for log file processors.

 Reads a log file from a Monte Carlo generator and extracts information.
 Inherited classes for each generator type implement the Extract()
 method to gather the required information for that generator and the
 Save() method to store it to file.
 */
class LogReader : public TObject {
 public:
  /**
   Constructor.
   */
  LogReader() { }

  /**
   Destructor.
   */
  virtual ~LogReader() { }

  /**
   Return a new LogReader instance.
   */
  virtual LogReader* Create() const = 0;

  /**
   Extract data from the named log file.
   */
  virtual bool Extract(const std::string& file) = 0;

  /**
   Saves the extracted data to the current file, if one is open and is
   writeable.
   Returns -1 if the data cannot be saved.
   To write the LogReader itself, use LogReader::Write().
   */
  virtual Int_t Save() const = 0;

  ClassDef(erhic::LogReader, 1)
};

/**
 Processes PYTHIA log files.

 Reads a log file and finds the total cross section and the number
 of generated events when Extract() is called.
 Writes those values as TObjStrings to the current directory when
 Save() is called, assuming that directory is writeable.
 */
class LogReaderPythia : public LogReader {
 public:
  /**
   Constructor.
   */
  LogReaderPythia();

  /**
   Destructor.
   */
  virtual ~LogReaderPythia();

  /**
   Return a new LogReaderPythia instance.
   */
  LogReaderPythia* Create() const;

  /**
   Extract data from the named log file.
   */
  bool Extract(const std::string& file);

  /**
   Write the extracted information to the current file, if it is
   writeable. If you want to write the LogReaderPythia itself, use
   LogReaderPythia::Write().
   */
  Int_t Save() const;

 protected:
  TObjString crossSection_;  ///> Total cross section in microbarns
  TObjString nEvents_;  ///> Total number of events generated

  ClassDef(erhic::LogReaderPythia, 1)
};

inline LogReaderPythia* LogReaderPythia::Create() const {
  return new LogReaderPythia;
}

/**
 Processes PEPSI log files.

 Reads a log file and finds the total cross section and the number
 of generated events when Extract() is called.
 Writes those values as TObjStrings to the current directory when
 Save() is called, assuming that directory is writeable.
 */
class LogReaderPepsi : public LogReader {
 public:
  /**
   Constructor.
   */
  LogReaderPepsi();

  /**
   Destructor.
   */
  virtual ~LogReaderPepsi();

  /**
   Return a new LogReaderPepsi instance.
   */
  LogReaderPepsi* Create() const;

  /**
   Extract data from the named log file.
   */
  bool Extract(const std::string& file);

  /**
   Write the extracted information to the current file, if it is
   writeable. If you want to write the LogReaderPepsi itself, use
   LogReaderPepsi::Write().
   */
  Int_t Save() const;

 protected:
  TObjString crossSection_;  ///> Total cross section in microbarns
  TObjString nEvents_;  ///> Total number of events generated

  ClassDef(erhic::LogReaderPepsi, 1)
};

inline LogReaderPepsi* LogReaderPepsi::Create() const {
  return new LogReaderPepsi;
}

/**
 Processes DJANGOH log files.

 Reads a log file and finds the total cross section and the number
 of generated events when Extract() is called.
 Writes those values as TObjStrings to the current directory when
 Save() is called, assuming that directory is writeable.
 */
class LogReaderDjangoh : public LogReader {
 public:
  /**
   Constructor.
   */
  LogReaderDjangoh();

  /**
   Destructor.
   */
  virtual ~LogReaderDjangoh();

  /**
   Return a new LogReaderDjangoh instance.
   */
  LogReaderDjangoh* Create() const;

  /**
   Extract data from the named log file.
   */
  bool Extract(const std::string& file);

  /**
   Write the extracted information to the current file, if it is
   writeable. If you want to write the LogReaderDjangoh itself, use
   LogReaderDjangoh::Write().
   */
  Int_t Save() const;

 protected:
  TObjString crossSection_;  ///> Total cross section in microbarns
  TObjString nEvents_;  ///> Total number of events generated

  ClassDef(erhic::LogReaderDjangoh, 1)
};

inline LogReaderDjangoh* LogReaderDjangoh::Create() const {
  return new LogReaderDjangoh;
}

/**
 Processes PYTHIA log files.

 Reads a log file and finds the total cross section and the number
 of generated events when Extract() is called.
 Writes those values as TObjStrings to the current directory when
 Save() is called, assuming that directory is writeable.
 */
class LogReaderMilou : public LogReader {
 public:
  /**
   Constructor.
   */
  LogReaderMilou() { }

  /**
   Destructor.
   */
  virtual ~LogReaderMilou() { }

  /**
   Return a new LogReaderMilou instance.
   */
  LogReaderMilou* Create() const;

  /**
   Extract data from the named log file.
   */
  bool Extract(const std::string& file);

  /**
   Write the extracted information to the current file, if it is
   writeable. If you want to write the LogReaderMilou itself, use
   LogReaderMilou::Write().
   */
  Int_t Save() const;

  /**
   Returns the number of events reported by the log file.
   Extract() should be called first.
   */
  Int_t GetNEvents() const;

  /**
   Returns the total cross section reported by the log file.
   Extract() should be called first.
   */
  Double_t GetCrossSection() const;

  /**
   Returns the error on total cross section reported by the log file.
   Extract() should be called first.
   */
  Double_t GetCrossSectionError() const;

 protected:
  TObjString crossSection_;  ///> Total cross section in nb
  TObjString crossSectionError_;  ///> Cross section error in nb
  TObjString nEvents_;  ///> Total number of events generated

  ClassDef(erhic::LogReaderMilou, 1)
};

inline LogReaderMilou* LogReaderMilou::Create() const {
  return new LogReaderMilou;
}

inline Int_t LogReaderMilou::GetNEvents() const {
  return nEvents_.GetString().Atoi();
}

inline Double_t LogReaderMilou::GetCrossSection() const {
  return crossSection_.GetString().Atof();
}

/**
 Returns the error on total cross section reported by the log file.
 Extract() should be called first.
 */
inline Double_t LogReaderMilou::GetCrossSectionError() const {
  return crossSectionError_.GetString().Atof();
}

/**
 Processes gmc_trans log files.

 Reads a log file and finds the total cross section and the number
 of generated events.
 */
class LogReaderGmcTrans : public LogReader {
 public:
  /**
   Constructor.
   */
  LogReaderGmcTrans();

  /**
   Destructor.
   */
  virtual ~LogReaderGmcTrans();

  /**
   Return a new LogReaderGmcTrans instance.
   */
  LogReaderGmcTrans* Create() const;

  /**
   Search the named file.
   Store the total cross section and number of generated events
   if they can be found.
   */
  bool Extract(const std::string& filename);

  /**
   Write the stored cross section and number of events to
   the active ROOT directory.
   Returns the total number of bytes written, or a value <= 0
   in the case of an error.
   */
  Int_t Save() const;

  /**
   Returns the number of events reported by the log file.
   Extract() should be called first.
   */
  Int_t GetNEvents() const;

  /**
   Returns the total cross section reported by the log file.
   Extract() should be called first.
   */
  Double_t GetCrossSection() const;

 protected:
  TObjString mNEvents;  ///< Number of generated events.
  TObjString mCrossSection;  ///< Total cross section in microbarns.

  ClassDef(erhic::LogReaderGmcTrans, 1)
};

/**
 Factory class for LogReaders.

 Singleton class.
 Creates a LogReader instance corresponding to a Monte Carlo
 generator type.
 */
class LogReaderFactory {
 public:
  /**
   Returns the single instance of LogReaderFactory.
   */
  static LogReaderFactory& GetInstance();

  /**
   Returns a LogReader instance of the type for reading log files
   from the Monte Carlo generator event type 'event'.
   Returns NULL in the case of an unsupported generator.
   The LogReader must be deleted by the user.
   */
  LogReader* CreateReader(const EventBase& event) const;

  /**
   Returns a LogReader instance of the type for reading log files
   from the Monte Carlo generator named 'name'.
   Returns NULL in the case of an unsupported generator.
   The LogReader must be deleted by the user.
   */
  LogReader* CreateReader(const std::string& name) const;

  /**
   Returns a LogReader instance of the type for reading log files
   from the Monte Carlo generator which produced the content in an
   istream, by reading the first line of that stream.
   Returns NULL in the case of an unsupported generator.
   The LogReader must be deleted by the user.
   */
  LogReader* CreateReader(std::istream&) const;

  /**
   Attempts to locate a log file corresponding to the named Monte
   Carlo file.
   Searches for a file with the same base name and extension '.log'.
   Looks in
   the current directory and, if mcFile gives a path containing
   'TXTFILES', in the corresonding directory substituting 'LOGFILES'.
   */
  std::string Locate(const std::string& mcFile) const;

 protected:
  /**
   Constructor.
   */
  LogReaderFactory();

  /**
   Destructor.
   */
  ~LogReaderFactory();

  typedef std::map<std::string, LogReader*> Map;
  Map prototypes_;

  ClassDef(erhic::LogReaderFactory, 1)
};

/**
 Abstract base class for Monte Carlo file types.
 Describes a Monte Carlo file type and returns objects required for
 processing or analysis of that file type.
 */
class FileType : public TObject {
 public:
  /**
   Destructor.
   */
  virtual ~FileType() { }

  /**
   Returns a new FileType instance.
   */
  virtual FileType* Create() const = 0;

  /**
   Returns a new event object for the generator making this type of file.
   */
  virtual EventBase* AllocateEvent() const = 0;

  /**
   Returns the name of the generator making this type of file.
   */
  virtual std::string GetGeneratorName() const = 0;

  /**
   Returns a reader to process the log file corresponding to this type of file.
   */
  virtual LogReader* CreateLogReader() const = 0;

  /**
   Returns a new event object for the generator making this type of file.
   */
  virtual VirtualEventFactory* CreateEventFactory(std::istream&) const = 0;

  ClassDef(erhic::FileType, 1)
};

/*
 Templated file descriptor class, valid for Monte Carlo event classes.
 e.g. File<EventPythia> describes a Pythia event file.
 */
template<typename T>
class File : public FileType {
 public:
  /**
   Constructor.
   
   If the string argument is not empty, the File attempts to open
   a file with that name. If the file is opened it tries to extract
   */
  File();

  /**
   Destructor.
   */
  virtual ~File();

  /**
   Returns a new File object.
   */
  virtual File<T>* Create() const;

  /**
   Allocates an event of the type for this file.
   */
  virtual T* AllocateEvent() const;

  /**
   Returns the name of the generator.
   Entirely in lower case.
   */
  virtual std::string GetGeneratorName() const;

  /**
   Create a LogReader for this type of Monte Carlo file.
   Returns NULL if the file type is unsupported or has no LogReader
   class implemented.
   The LogReader must be deleted by the user.
   */
  virtual LogReader* CreateLogReader() const;

  /**
   Returns a new event factory instance.
   */
  virtual EventFromAsciiFactory<T>*
  CreateEventFactory(std::istream& is) const {
     return new EventFromAsciiFactory<T>(is);
  }

 protected:
  T* t_;

  // Warning: explicitly putting the erhic:: namespace before the class
  // name doesn't seen to work for template classes.
  ClassDef(File, 1)
};

template<typename T>
inline T* File<T>::AllocateEvent() const {
  return new T;
}

template<typename T>
inline File<T>* File<T>::Create() const {
  return new File<T>();
}

/**
 Factory class for Files.
 Singleton class.
 */
class FileFactory {
 public:
  /**
   Returns the single instance of FileFactory.
   */
  static FileFactory& GetInstance();

  /**
   Returns a FileType object for the named generator.
   */
  const FileType* GetFile(const std::string& generatorName) const;

  /**
   Returns a FileType object, determining the generator type from a stream.
   */
  const FileType* GetFile(std::istream&) const;

 protected:
  /**
   Constructor.
   */
  FileFactory();

  /**
   Destructor.
   */
  virtual ~FileFactory();

  typedef std::map<std::string, FileType*> Map;
  Map prototypes_;
};

}  // namespace erhic

#endif  // INCLUDE_EICSMEAR_ERHIC_FILE_H_
