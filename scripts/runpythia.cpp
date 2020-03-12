// Example for generating EventPythia objects directly from PYTHIA
// output, without intermediate text output.

#include<iostream>
#include<string>
#include<TString.h>
#include<TPythia6.h>
#include<TRandom3.h>
#include<TFile.h>
#include<TStopwatch.h>

#include <eicsmear/erhic/Pythia6.h>
#include <eicsmear/erhic/Pythia6EventBuilder.h>

using namespace std;


void runpythia(TString outFile,
	       int nEvents,
	       double pElectron,
	       double pProton,
	       double minQ2 = 1.,
	       double maxQ2 = -1.,
	       int messageInterval = 1000){
  // bool doFiltering = false,
  // double minz = 0.,
  // double maxcorrelation = 2.)
   
   cout << "Input arguments:" << endl;
   cout << "Output file: " << outFile << endl;
   cout << "Number of events: " << nEvents << endl;
   cout << "Electron momentum: " << pElectron << endl;
   cout << "Proton momentum: " << pProton << endl;
   cout << "Minimum Q2: " << minQ2 << endl;
   cout << "Maximum Q2: " << maxQ2 << endl;
   cout << "Message interval: " << messageInterval << endl;
   // cout << "Filtering: " << doFiltering << endl;
   // cout << "Minimum z: " << minz << endl;
   // cout << "Max pair correlation: " << maxcorrelation << endl;
   
   // The behaviour of PYTHIA 6 is configured via the ROOT interface
   // class TPythia6.
   TPythia6* pythia = TPythia6::Instance();
   
   // Configure for charm production
   // I copied these parameters from Elke's charm log files
   pythia->SetMSEL(4);
   pythia->SetMSTP(14, 30);
   pythia->SetMSTP(15, 0);
   pythia->SetMSTP(16, 1);
   pythia->SetMSTP(17, 4); // MSTP 17=6 is the R-rho measured as by hermes, =4 Default
   pythia->SetMSTP(18, 3);
   pythia->SetMSTP(19, 1); // Hermes MSTP-19=1 different Q2 suppression, default = 4
   pythia->SetMSTP(20, 0); // Hermes MSTP(20, 0 , default MSTP(20, 3
   pythia->SetMSTP(32, 8);
   pythia->SetMSTP(38, 4);
   pythia->SetMSTP(51, 10800); // PDFLIB/LHAPDF 10800 = CT10
   pythia->SetMSTP(52, 2); // MSTP(52) = 2 indicates to use PDFLIB
   pythia->SetMSTP(53, 3);
   pythia->SetMSTP(54, 1);
   pythia->SetMSTP(55, 5);
   pythia->SetMSTP(56, 1);
   pythia->SetMSTP(57, 1);
   pythia->SetMSTP(58, 5);
   pythia->SetMSTP(59, 1);
   pythia->SetMSTP(60, 7);
   pythia->SetMSTP(61, 2);
   pythia->SetMSTP(71, 1);
   pythia->SetMSTP(81, 0);
   pythia->SetMSTP(82, 1);
   pythia->SetMSTP(91, 1);
   pythia->SetMSTP(92, 3);      // hermes MSTP(92, 4
   pythia->SetMSTP(93, 1);
   pythia->SetMSTP(101, 3);
   pythia->SetMSTP(102, 1);
   pythia->SetMSTP(111, 1);
   pythia->SetMSTP(121, 0);
   pythia->SetPARP(13, 1);
   pythia->SetPARP(18, 0.4); // hermes pythia->SetPARP(18, 0.17
   pythia->SetPARP(81, 1.9);
   pythia->SetPARP(89, 1800);
   pythia->SetPARP(90, 0.16);
   pythia->SetPARP(91, 0.40);
   pythia->SetPARP(93, 5.);
   pythia->SetPARP(99, 0.40);
   pythia->SetPARP(100, 5);
   pythia->SetPARP(102, 0.28);
   pythia->SetPARP(103, 1.0);
   pythia->SetPARP(104, 0.8);
   pythia->SetPARP(111, 2.);
   pythia->SetPARP(161, 3.00);
   pythia->SetPARP(162, 24.6);
   pythia->SetPARP(163, 18.8);
   pythia->SetPARP(164, 11.5);
   pythia->SetPARP(165, 0.47679);
   pythia->SetPARP(166, 0.67597); // PARP165/166 are linked to MSTP17 as R_rho of HERMES is used
   pythia->SetPARJ(1, 0.100);
   pythia->SetPARJ(2, 0.300);
   pythia->SetPARJ(11, 0.5);
   pythia->SetPARJ(12, 0.6);
   pythia->SetPARJ(21,  0.40);
   pythia->SetPARJ(32, 1.0);
   pythia->SetPARJ(33,  0.80);
   pythia->SetPARJ(41,  0.30);
   pythia->SetPARJ(42,  0.58);
   pythia->SetPARJ(45,  0.5);
   pythia->SetMSTJ(1, 1);
   pythia->SetMSTJ(12, 1);
   pythia->SetMSTJ(45, 5);
   pythia->SetMSTU(112, 5);
   pythia->SetMSTU(113, 5);
   pythia->SetMSTU(114, 5);
   pythia->SetCKIN(1, 1.);
   pythia->SetCKIN(2, -1.);
   pythia->SetCKIN(3, 0.);
   pythia->SetCKIN(4, -1.);
   pythia->SetCKIN(5, 1.00);
   pythia->SetCKIN(6, 1.00);
   pythia->SetCKIN(7, -10.);
   pythia->SetCKIN(8, 10.);
   pythia->SetCKIN(9, -40.);
   pythia->SetCKIN(10, 40.);
   pythia->SetCKIN(11, -40.);
   pythia->SetCKIN(12, 40.);
   pythia->SetCKIN(13, -40.);
   pythia->SetCKIN(14, 40.);
   pythia->SetCKIN(15, -40.);
   pythia->SetCKIN(16, 40.);
   pythia->SetCKIN(17, -1.);
   pythia->SetCKIN(18, 1.);
   pythia->SetCKIN(19, -1.);
   pythia->SetCKIN(20, 1.);
   pythia->SetCKIN(21, 0.);
   pythia->SetCKIN(22, 1.);
   pythia->SetCKIN(23, 0.);
   pythia->SetCKIN(24, 1.);
   pythia->SetCKIN(25, -1.);
   pythia->SetCKIN(26, 1.);
   pythia->SetCKIN(27, -1.);
   pythia->SetCKIN(28, 1.);
   pythia->SetCKIN(31, 2.);
   pythia->SetCKIN(32, -1.);
   pythia->SetCKIN(35, 0.);
   pythia->SetCKIN(36, -1);
   pythia->SetCKIN(37, 0.);
   pythia->SetCKIN(38, -1.);
   pythia->SetCKIN(39, 4.);
   pythia->SetCKIN(40, -1.);
   pythia->SetCKIN(65, minQ2);        // Min for Q^2
   pythia->SetCKIN(66, maxQ2);       // Max for Q^2
   pythia->SetCKIN(67, 0.);
   pythia->SetCKIN(68, -1.);
   pythia->SetCKIN(77, 2.0);
   pythia->SetCKIN(78, -1.);
   
   // Set beam momenta.
   // The fourth argument is a dummy value when calling with frame "3MOM".
   // Set the electron 3-momenta via SetP(1, ...) (beam 1)
   // and the proton 3-momenta via SetP(2, ...) (beam 2)
   pythia->SetP(1, 1, 0.);  // px
   pythia->SetP(1, 2, 0.);  // py
   pythia->SetP(1, 3, -pElectron); // pz
   pythia->SetP(2, 1, 0.);   // px
   pythia->SetP(2, 2, 0.);   // py
   pythia->SetP(2, 3, pProton); // pz
   
   // Set the random number generator seed.
   // Seeds can be in the range [0, 900,000,000].
   // TRandom::Integer(N) returns number in the range [0, N).
   // Set the TRandom seed to zero so the seed is set to a random value.
   TRandom3 random(0);
   pythia->SetMRPY(1, random.Integer(900000001));
   std::cout << "PYTHIA random number seed MRPY(1) = " << pythia->GetMRPY(1)
   << std::endl;
   
   // Intialise PYTHIA for e-p collisions
   char* target = (char*) "p+";
   char* beam = (char*) "gamma/e-";
   char* frame = (char*) "3MOM";
   float WIN = pElectron;
   pythia->Pyinit(frame, beam, target, WIN);
   
   // Open the output   
   TFile* file = new TFile(outFile, "RECREATE");
   erhic::VirtualEventFactory* factory = new erhic::Pythia6EventBuilder();
   erhic::Pythia6 pythia6(file, factory, nEvents, "EICTree", "event", messageInterval);
   
   // The event filter code uses the particle filter
   // classes for the dihadron code.
   // The file is provided but untested.
   // gROOT->LoadMacro("filter.cpp+g"); // Defines class DiHadronEventFilter
   // DiHadronEventFilter filter(minz, maxcorrelation);
   // if(doFiltering) pythia6.SetFilter(&filter);
   
   // Run PYTHIA and write output
   TStopwatch timer;
   pythia6.Run();
      
   std::cout << "Completed in " << timer.RealTime() << " seconds" << std::endl;
   file->Close();

   return;   
}
