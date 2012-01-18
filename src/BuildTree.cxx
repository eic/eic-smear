#include <string>

#include <TString.h>
#include <TSystem.h>

#include "Forester.h"
#include "File.h"

/**
 This is an example function to generate ROOT files.
 It can be used "out of the box".
 If more control over the output is desired, then the settings of the
 Forester can be tweaked to do so.
 */
Long64_t
BuildTree(const TString& inputFileName,
          const TString& outputDirName,
          const Long64_t maxEvent,
          std::string logFile) {
   
	//	Set the maximum size of the tree on disk.
	// Once this size is reached a new file is opened for continued writing.
   // Set 10 Gb. Us LL to force long integer.
	TTree::SetMaxTreeSize(10LL * 1024LL * 1024LL * 1024LL);
   
   // Get the input file name, stripping any leading directory path via
   // use of the BaseName() method from TSystem.
   TString outName = gSystem->BaseName(inputFileName);
   
   // Remove the existing extension, if there is one.
   if(outName.Last('.') > -1) {
      outName.Replace(outName.Last('.'), outName.Length(), "");
   } // if
   
   // If we are analysing a subset of events, include the number of events in
   // the file name before the extension.
   if(maxEvent > 0) {
      outName.Append(".");
      outName += maxEvent;
      outName.Append("event");
   } // if
   
   outName.Append(".root");
   
   TString outDir(outputDirName);
   if(not outDir.EndsWith("/")) outDir.Append('/');
   outName.Prepend(outDir);
   
   // Configure an object of class Forester, which handles processing the text
   // file into a tree.
   erhic::Forester forester;
   forester.SetInputFileName(std::string(inputFileName));
   forester.SetOutputFileName(std::string(outName));
   forester.SetMaxNEvents(maxEvent);
   forester.SetMessageInterval(10000);
   forester.SetBeVerbose(true);
   
   Long64_t result = forester.Plant(); // Plant that tree!
   
   // Search the log file for information.
   // Use the provided log file name if there is one, otherwise attempt
   // automated procedure to locate it.
   if(logFile.empty()) {
#ifdef USE_NAMESPACE_ERHIC
      using namespace erhic::monte_carlo;
#endif
      logFile =
      LogReaderFactory::GetInstance().Locate(inputFileName.Data());
   } // if
   
   // Use the FileType created by Forester when running to generate a
   // LogReader to process the log file, assuming a FileType was created and
   // the log file was located.
   if(forester.GetFileType() and not logFile.empty()) {
      TFile rootFile(outName, "UPDATE");
      erhic::monte_carlo::LogReader* reader =
         forester.GetFileType()->CreateLogReader();
      bool wasRead = (reader ? reader->Extract(logFile) : false);
      if(wasRead) {
//         std::cout <<"wrote "<<reader->Save()<<" bytes"<<std::endl;
         reader->Save();
      } // if
        //      else std::cout <<"Failed to read log file "<<logFile<<std::endl;
   } // if
     //   else std::cout << "Couldn't find the log file "<< forester.GetFileType() << " " << logFile <<std::endl;
   
   return result;
}