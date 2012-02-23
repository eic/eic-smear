/**
 Pythia6.h

 \file
 Declaration of class Pythia6.

 Created by TB on 1/18/12.
 Copyright 2012 BNL. All rights reserved.
*/

#ifndef _Pythia6_H_
#define _Pythia6_H_

#include <string>

#include <Rtypes.h> // For ClassDef macro

class TFile;
class TTree;

namespace erhic {
   
   class EventPythia;
   class EventMCFilterABC;
   
   /**
    Runs PYTHIA 6 and builds an output tree file in the ROOT environment.
    The output tree is populated with EventPythia objects.
    PYTHIA options (including initialisation) should be provided via the
    ROOT interface class TPythia6 before running the tree-building.
    */
   class Pythia6 {
      
   public:
      
      /**
       Constructor.
       
       Associates the output of this PYTHIA run with a file.
       Define the number of events to produce, and optionally
       provide a tree and branch name.
       */
      Pythia6(TFile* file,
              int nEvents,
              const std::string& treeName = "EICTree",
              const std::string& branchName = "event");
      
      /** Destructor */
      virtual ~Pythia6();
      
      /** Optionally specify a filter to apply to events.
       Events failing this filter will not be stored.
       By default no filtering is applied.
       The passed filter object is not deleted by Pythia6. */
      virtual void SetFilter(EventMCFilterABC*);
      
      /** Runs PYTHIA and writes output */
      virtual bool Run();
      
   protected:
      
      int mPrintInterval;
      TFile* mFile;           ///< Pointer to the output file
      TTree* mTree;           ///< Pointer to the generated tree
      EventPythia* mEvent;    ///< Pointer to the event buffer
      const int mNEvents;     ///< Number of events to produce
      int mNGenerated;        ///< Number of events passing PYTHIA selection
      int mNTrials;           ///< Number of events thrown by PYTHIA
      EventMCFilterABC* mFilter; ///< Event filter
      
      ClassDef(Pythia6, 1)
   };
   
   inline void Pythia6::SetFilter(EventMCFilterABC* filter) {
      mFilter = filter;
   }
   
} // namespace erhic

#endif
