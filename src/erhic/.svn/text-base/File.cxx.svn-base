/**
 \file
 Implementations of file information and log-reading classes.
 
 \author    Thomas Burton
 \date      2011-07-29
 \copyright 2011 Brookhaven National Lab
 */

#include "eicsmear/erhic/File.h"

#include <fstream>
#include <sstream>
#include <string>

#include <TSystem.h>

#include "eicsmear/erhic/EventPepsi.h"
#include "eicsmear/erhic/EventDjangoh.h"
#include "eicsmear/erhic/EventDpmjet.h"
#include "eicsmear/erhic/EventRapgap.h"
#include "eicsmear/erhic/EventGmcTrans.h"

namespace erhic {

LogReaderPythia::LogReaderPythia() { }

LogReaderPythia::~LogReaderPythia() { }

bool LogReaderPythia::Extract(const std::string& file) {
  // First we get the cross section from the log file.
  std::ifstream ifs(file.c_str(), std::ios::in);

  if (!ifs.is_open()) return false;

  std::string line;
  const std::string searchPattern("Pythia total cross section normalisation:");
  std::string normalisation;
  std::string nEvents;

  while (ifs.good()) {
    std::getline(ifs, line);
    size_t position = line.find(searchPattern);
    if (position != std::string::npos) {
      // We found the line.
      // Erase the text preceding the cross section.
      std::stringstream ss;
      line.erase(0, position + searchPattern.length());
      ss.str("");
      ss.clear();
      ss << line;
      ss >> normalisation;
    }  // if
    const std::string searchPattern2("Total Number of generated events");
    position = line.find(searchPattern2);
    if (position != std::string::npos) {
      // We found the line.
      // Erase the text preceding the cross section.
      std::stringstream ss;
      line.erase(0, position + searchPattern2.length());
      ss.str("");
      ss.clear();
      ss << line;
      ss >> nEvents;
    }  // if
  }  // while
  crossSection_.SetString(normalisation.c_str());
  std::cout << crossSection_.GetString().Atof() << std::endl;
  nEvents_.SetString(nEvents.c_str());
  std::cout << nEvents_.GetString().Atoi() << std::endl;
  std::cout << "Extracted information from " << file << std::endl;
  return true;
}

Int_t LogReaderPythia::Save() const {
  return
  crossSection_.Write("crossSection") +
  nEvents_.Write("nEvents");
  nEvents_.Write("nTrials");
}

//
// class LogReaderPepsi
//

LogReaderPepsi::LogReaderPepsi() { }

LogReaderPepsi::~LogReaderPepsi() { }

bool LogReaderPepsi::Extract(const std::string& file) {
  // First we get the cross section from the log file.
  std::ifstream ifs(file.c_str(), std::ios::in);
  if (!ifs.is_open()) return false;
  std::string line;
  const std::string searchPattern(
      "total cross section in pb from MC simulation");
  std::string normalisation;
  std::string nEvents;
  while (ifs.good()) {
    std::getline(ifs, line);
    size_t position = line.find(searchPattern);
    if (position != std::string::npos) {
      // We found the line.
      // Erase the text preceding the cross section.
      std::stringstream ss;
      line.erase(0, position + searchPattern.length());
      ss.str("");
      ss.clear();
      ss << line;
      double value;
      // Divide by 1,000,000 to go from pb to microbarn
      ss >> value;
      value /= 1.e6;
      ss.str("");
      ss.clear();
      ss << value;
      ss >> normalisation;
    }  // if
    const std::string searchPattern2("Total Number of trials");
    position = line.find(searchPattern2);
    if (position != std::string::npos) {
      // We found the line.
      // Erase the text preceding the cross section.
      std::stringstream ss;
      line.erase(0, position + searchPattern2.length());
      ss.str("");
      ss.clear();
      ss << line;
      ss >> nEvents;
    }  // if
  }  // while
  crossSection_.SetString(normalisation.c_str());
  std::cout << crossSection_.GetString().Atof() << std::endl;
  nEvents_.SetString(nEvents.c_str());
  std::cout << nEvents_.GetString().Atoi() << std::endl;
  std::cout << "Extracted information from " << file << std::endl;
  return true;
}

Int_t LogReaderPepsi::Save() const {
  return
  crossSection_.Write("crossSection") +
  nEvents_.Write("nEvents");
}

//
// class LogReaderDjangoh
//

LogReaderDjangoh::LogReaderDjangoh() { }

LogReaderDjangoh::~LogReaderDjangoh() { }

bool LogReaderDjangoh::Extract(const std::string& file) {
  // First we get the cross section from the log file.
  std::ifstream ifs(file.c_str(), std::ios::in);

  if (!ifs.is_open()) return false;

  std::string line;
  // There are two places we need to look for cross sections.
  // As standard, look for the line starting with
  //   Cross-section from HERACLES
  // Sometimes there will be an additional line starting
  //   Total cross section is now    SIGTOT =
  // *If* this line is present we take the cross section from there,
  // otherwise use the 'HERACLES' line, which is always present.
  // Both lines give cross section in pb.
  const std::string xsecPattern("Cross-section from HERACLES =");
  const std::string xsecAltPattern("Total cross section is now    SIGTOT =");
  std::string normalisation;
  std::string nEvents;

  while (ifs.good()) {
    std::getline(ifs, line);
    // Look for normal cross section line, unless the cross section
    // was already set (via the alternative line, see below).
    size_t position = line.find(xsecPattern);
    if (position != std::string::npos && normalisation.empty()) {
      // We found the line.
      // Erase the text preceding the cross section.
      std::stringstream ss;
      line.erase(0, position + xsecPattern.length());
      ss << line;
      double value;
      // Divide by 1,000,000 to go from pb to microbarn
      ss >> value;
      value /= 1.e6;
      ss.str("");
      ss.clear();
      ss << value;
      ss >> normalisation;
    }  // if
    position = line.find(xsecAltPattern);
    // Look for alternative cross section.
    if (position != std::string::npos) {
      // We found the line.
      // Erase the text preceding the cross section.
      std::stringstream ss;
      line.erase(0, position + xsecAltPattern.length());
      ss << line;
      double value;
      // Divide by 1,000,000 to go from pb to microbarn
      ss >> value;
      value /= 1.e6;
      ss.str("");
      ss.clear();
      ss << value;
      ss >> normalisation;
    }  // if
    const std::string searchPattern2("TOTAL EVENT NUMBER");
    position = line.find(searchPattern2);
    if (position != std::string::npos) {
      // We found the line.
      // Erase the text preceding the cross section.
      std::stringstream ss;
      line.erase(0, position + searchPattern2.length());
      ss.str("");
      ss.clear();
      ss << line;
      ss >> nEvents;
    }  // if
  }  // while

  crossSection_.SetString(normalisation.c_str());
  std::cout << crossSection_.GetString().Atof() << std::endl;
  nEvents_.SetString(nEvents.c_str());
  std::cout << nEvents_.GetString().Atoi() << std::endl;

  std::cout << "Extracted information from " << file << std::endl;
  return true;
}

Int_t LogReaderDjangoh::Save() const {
  return
  crossSection_.Write("crossSection") +
  nEvents_.Write("nEvents");
}


bool LogReaderMilou::Extract(const std::string& file) {
  // First we get the cross section from the log file.
  std::ifstream ifs(file.c_str(), std::ios::in);

  if (!ifs.is_open()) return false;

  std::string line;

  const std::string nEventPattern("Number of generated events    =");
  const std::string xsecPattern("Total cross-section (nb) :");
  const std::string errorPattern("Error                    :");

  std::string tmp("");

  std::stringstream ss;

  while (ifs.good()) {
    std::getline(ifs, line);
    ss.str("");
    ss.clear();
    // Look for one of the search patterns.
    if (line.find(nEventPattern) != std::string::npos) {
      line.erase(0,
                 line.find(nEventPattern) + nEventPattern.length());
      ss << line;
      ss >> tmp;
      nEvents_.SetString(tmp.c_str());
    } else if (line.find(xsecPattern) != std::string::npos) {
      line.erase(0,
                 line.find(xsecPattern) + xsecPattern.length());
      ss << line;
      ss >> tmp;
      crossSection_.SetString(tmp.c_str());
    } else if (line.find(errorPattern) != std::string::npos) {
      line.erase(0,
                 line.find(errorPattern) + errorPattern.length());
      ss << line;
      ss >> tmp;
      crossSectionError_.SetString(tmp.c_str());
    }  // if
  }  // while
  // Return true if all the strings are filled, or false if any
  // of them are empty.
  return !(nEvents_.GetString().IsNull() ||
              crossSection_.GetString().IsNull() ||
              crossSectionError_.GetString().IsNull());
}

Int_t LogReaderMilou::Save() const {
  Int_t bytes(0);
  bytes += nEvents_.Write("nEvents");
  bytes += crossSection_.Write("crossSection");
  bytes += crossSectionError_.Write("crossSectionError");
  return bytes;
}

LogReaderGmcTrans::LogReaderGmcTrans()
: mNEvents("")
, mCrossSection("") {
}

LogReaderGmcTrans::~LogReaderGmcTrans() {
}

LogReaderGmcTrans* LogReaderGmcTrans::Create() const {
  return new LogReaderGmcTrans;
}

bool LogReaderGmcTrans::Extract(const std::string& filename) {
  // Open the file, check for errors.
  std::ifstream file(filename.c_str(), std::ios::in);
  if (!file.is_open()) {
    return false;
  }  // if
  // The line with the number of generated events ends with this:
  const std::string eventPattern("generated events (IEVGEN)");
  // Thw line with the total cross section ends with this:
  const std::string xsecPattern("total xsec in microbarns after selector");
  // Read the file line-by-line, looking for the patterns.
  std::stringstream sstream;  // For splitting the string
  std::string line;
  while (file.good()) {
    std::getline(file, line);
    // Check for number-of-events pattern.
    if (line.find(eventPattern) != std::string::npos) {
      // Found it, the number of events is the first word in the line.
      std::string tmp;
      sstream.str("");
      sstream.clear();
      sstream << line;
      sstream >> tmp;
      mNEvents.SetString(tmp.c_str());
    }  // if
    // Check for total cross section pattern.
    if (line.find(xsecPattern) != std::string::npos) {
      std::string tmp;
      sstream.str("");
      sstream.clear();
      sstream << line;
      sstream >> tmp;
      mCrossSection.SetString(tmp.c_str());
    }  // if
  }  // while
  return !(mNEvents.GetString().IsNull() ||
           mCrossSection.GetString().IsNull());
}

Int_t LogReaderGmcTrans::Save() const {
  return mNEvents.Write("nEvents") + mCrossSection.Write("crossSection");
}

Int_t LogReaderGmcTrans::GetNEvents() const {
  return mNEvents.GetString().Atoi();
}

Double_t LogReaderGmcTrans::GetCrossSection() const {
  return mCrossSection.GetString().Atof();
}

LogReaderFactory& LogReaderFactory::GetInstance() {
  static LogReaderFactory theInstance;
  return theInstance;
}

/**
 Returns a LogReader instance of the type for reading log files
 from the Monte Carlo generator event type 'event'.
 Returns NULL in the case of an unsupported generator.
 The LogReader must be deleted by the user.
 */
LogReader* LogReaderFactory::CreateReader(const EventBase& event) const {
  // The event name will be "EventX" where "X" is the Monte Carlo
  // generator name.
  TString name = event.ClassName();
  name.ReplaceAll("erhic::", "");
  name.ReplaceAll("Event", "");
  name.ToLower();
  return CreateReader(name.Data());
}

/**
 Returns a LogReader instance of the type for reading log files
 from the Monte Carlo generator named 'name'.
 Returns NULL in the case of an unsupported generator.
 The LogReader must be deleted by the user.
 */
LogReader* LogReaderFactory::CreateReader(const std::string& name) const {
  // Use TString::ToLower() to convert the input name to all
  // lower case.
  TString str(name);
  str.ToLower();
  LogReader* reader(NULL);
  if (prototypes_.find(str.Data()) != prototypes_.end()) {
    reader = prototypes_.find(str.Data())->second->Create();
  }  // if
  return reader;
}

/**
 Returns a LogReader instance of the type for reading log files
 from the Monte Carlo generator named 'name'.
 Returns NULL in the case of an unsupported generator.
 The LogReader must be deleted by the user.
 */
LogReader* LogReaderFactory::CreateReader(std::istream& is) const {
  std::string line;
  std::getline(is, line);
  // Use TString::ToLower() to convert the input name to all
  // lower case.
  TString str(line);
  str.ToLower();
  LogReader* reader(NULL);
  if (str.Contains("pythia")) {
    reader = CreateReader("pythia");
  } else if (str.Contains("pepsi") || str.Contains("lepto")) {
    reader = CreateReader("pepsi");
  } else if (str.Contains("rapgap")) {
    reader = CreateReader("rapgap");
  } else if (str.Contains("djangoh")) {
    reader = CreateReader("djangoh");
  } else if (str.Contains("milou")) {
    reader = CreateReader("milou");
  }  // if
  return reader;
}

std::string LogReaderFactory::Locate(const std::string& mcFile) const {
  TString inFileName(mcFile);
  TString logFileName;
  std::string extension;
  if (mcFile.find_last_of('.') != std::string::npos) {
    extension = mcFile.substr(mcFile.find_last_of('.'));
  }  // if
  // If the input file is in the current directory, expand the full path:
  if (std::string(".") == gSystem->DirName(inFileName)) {
    inFileName.Prepend("/").Prepend(gSystem->pwd());
  }  // if

  // The standard data directory structure is to have a directory called
  // TXTFILES containing the Monte Carlo output, and a directory called
  // LOGFILES containing the corresponding log files. The sub-directory
  // structure of LOGFILES should match that of TXTFILES.
  // So, first we'll check if the input path contains TXTFILES, in which
  // case we just substitute LOGFILES:
  if (inFileName.Contains("TXTFILES")) {
    logFileName = inFileName;
    logFileName.ReplaceAll("TXTFILES", "LOGFILES");
    logFileName.ReplaceAll(extension.c_str(), ".log");
  }  // if...
  // Check if the file whose name we have constructed exists.
  // If not clear the string contents.
  if (gSystem->AccessPathName(logFileName, kFileExists)) {
    logFileName.Clear();
  }  // if
  // OK, that didn't work, so let's just look in the current directory
  if (logFileName.IsNull()) {
    logFileName = inFileName;
    if (extension.empty()) {
      logFileName.Append(".log");
    } else {
      logFileName.ReplaceAll(extension.c_str(), ".log");
    } // if
  }  // if
  // Check if the file whose name we have constructed exists.
  // If not clear the string contents.
  if (gSystem->AccessPathName(logFileName, kFileExists)) {
    logFileName.Clear();
  }  // if
  return logFileName.Data();
}

LogReaderFactory::LogReaderFactory() {
  prototypes_.insert(std::make_pair("pythia",
                                    new LogReaderPythia));
  prototypes_.insert(std::make_pair("milou",
                                    new LogReaderMilou));
  prototypes_.insert(std::make_pair("pepsi",
                                    new LogReaderPepsi));
  prototypes_.insert(std::make_pair("djangoh",
                                    new LogReaderDjangoh));
  prototypes_.insert(std::make_pair("gmctrans", new LogReaderGmcTrans));
}

LogReaderFactory::~LogReaderFactory() {
  Map::iterator i = prototypes_.begin();
  for (; i != prototypes_.end(); ++i) {
    delete i->second;
  }  // for
}

template<typename T>
File<T>::File() : t_(new T) { }

template<typename T>
File<T>::~File() {
  if (t_) {
    delete t_;
    t_ = NULL;
  }  // if
}

template<typename T>
std::string File<T>::GetGeneratorName() const {
  // The event class name is "EventX" where "X" is the generator
  // name.
  TString name = t_->ClassName();
  name.ReplaceAll("erhic::", "");
  name.ReplaceAll("Event", "");
  name.ToLower();
  return name.Data();
}

template<typename T>
LogReader* File<T>::CreateLogReader() const {
  return LogReaderFactory::GetInstance().CreateReader(*t_);
}

FileFactory& FileFactory::GetInstance() {
  static FileFactory theInstance;
  return theInstance;
}

const FileType* FileFactory::GetFile(const std::string& name) const {
  const FileType* file(NULL);
  if (prototypes_.find(name) != prototypes_.end()) {
    file = prototypes_.find(name)->second->Create();
  }  // if
  return file;
}

const FileType* FileFactory::GetFile(std::istream& is) const {
  std::string line;
  std::getline(is, line);
  // Use TString::ToLower() to convert the input name to all
  // lower case.
  TString str(line);
  str.ToLower();
  const FileType* file(NULL);
  if (str.Contains("pythia")) {
    file = GetFile("pythia");
  } else if (str.Contains("pepsi") || str.Contains("lepto")) {
    file = GetFile("pepsi");
  } else if (str.Contains("rapgap")) {
    file = GetFile("rapgap");
  } else if (str.Contains("djangoh")) {
    file = GetFile("djangoh");
  } else if (str.Contains("milou")) {
    file = GetFile("milou");
  } else if (str.Contains("gmctrans")) {
    file = GetFile("gmctrans");
  }  // if
  return file;
}

FileFactory::FileFactory() {
  prototypes_.insert(std::make_pair("djangoh",
                                    new File<EventDjangoh>()));
  prototypes_.insert(std::make_pair("dpmjet",
                                    new File<EventDpmjet>()));
  prototypes_.insert(std::make_pair("milou",
                                    new File<EventMilou>()));
  prototypes_.insert(std::make_pair("pepsi",
                                    new File<EventPepsi>()));
  prototypes_.insert(std::make_pair("pythia",
                                    new File<EventPythia>()));
  prototypes_.insert(std::make_pair("rapgap",
                                    new File<EventRapgap>()));
  prototypes_.insert(std::make_pair("gmctrans",
                                    new File<EventGmcTrans>()));
}

FileFactory::~FileFactory() {
  Map::iterator i = prototypes_.begin();
  for (; i != prototypes_.end(); ++i) {
    if (i->second) delete i->second;
  }  // for
}

}  // namespace erhic
