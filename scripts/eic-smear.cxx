// Authors: Wouter Deconinck, Kolja Kauder

// ROOT headers
#include <TROOT.h>
#include <TRint.h>
#include <TSystem.h>
#include <TThread.h>
#include <TString.h>

#include <iostream>
#include <string>

#include "eicsmear/functions.h"

R__EXTERN class SmearRint* gSmearRint;

class SmearRint : public TRint {

protected:

  static SmearRint* fExists;       ///< Check whether interface already existing

public:
  /// \brief Constructor
  SmearRint (const char* appClassName, int* argc, char** argv,
	     void* options = 0, int numOptions = 0, bool noLogo = kFALSE)
    : TRint (appClassName, argc, argv, options, numOptions, noLogo) {
    gSmearRint = this;

    // eic-smear command prompt
    SetPrompt("eic-smear [%d] ");

    // Pointer to self
    fExists = this;
  };

  /// \brief Destructor
  virtual ~SmearRint() {
    // Reset point to self
    if (fExists == this)
      fExists = NULL;
  };

}; // class smearRint

// Global pointers
SmearRint* gSmearRint = NULL;

// Pointer to self
SmearRint* SmearRint::fExists = NULL;

int main(int argc, char** argv)
{
  // Load library, assume in LD_LIBRARY_PATH
  gSystem->Load("libeicsmear");
  std::cout << "Using eic-smear version: " << erhic::EicSmearVersionString << std::endl;
  if ( argc >=2 ){
    std::string a1 = argv[1];
    if ( a1 == "-v" || a1 == "--version"){      
      return 0;
    }
  }    
  auto ErrorIgnoreLevel=gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal; 
  auto smeardets=gSystem->Load("libeicsmeardetectors");
  gErrorIgnoreLevel = ErrorIgnoreLevel;
  
  auto libeicsmearpath=gSystem->DynamicPathName("libeicsmear");
  // auto libdetectorpath="libeicsmeardetectors not found";
  auto libdetectorpath="";
  if ( smeardets>0 ) {
      libdetectorpath=gSystem->DynamicPathName("libeicsmeardetectors");
    }
  std::cout << "Using these eic-smear libraries : "<< std::endl
	    << libeicsmearpath << std::endl
	    << libdetectorpath << std::endl;

  // Start command prompt
  SmearRint* rint = new SmearRint("eic-smear ", &argc, argv);
  rint->Run();
  delete rint;
}
