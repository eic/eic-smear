// wrapper to compile runpythia.cpp

#include<iostream>
#include<string>
#include<sstream>

#include<TString.h>

using namespace std;

void runpythia(TString outFile,
	       int nEvents,
	       double pElectron,
	       double pProton,
	       double minQ2,
	       double maxQ2,
	       int messageInterval );

	       // bool doFiltering,
	       // double minz,
	       // double maxcorrelation);

int main ( int argc, char** argv){
  TString outFile;
  int nEvents=0;
  double pElectron=10;
  double pProton=100;
  double minQ2 = 1.;
  double maxQ2 = -1.;
  int messageInterval = 1000;
  // bool doFiltering = false;
  // double minz = 0.;
  // double maxcorrelation = 2.;

  // string helpstring="Usage: compiled_runpythia outFile nEvents pElectron pProton [minQ2] [maxQ2] [messageInterval] [doFiltering] [minz] [maxcorrelation]";
  string helpstring="Usage: compiled_runpythia outFile nEvents pElectron pProton [minQ2] [maxQ2] [messageInterval]";
  if ( argc >1 && string(argv[1]) == "-h" ) {
    cout << helpstring << endl;
    return -1;
  }
  
  // for (int i=0; i<argc; ++i ){
  //   cout << i << "  " << argv[i] << endl;
  // }
  
  switch ( argc ){
  // case 12:
  //   maxcorrelation = atof(argv[11]);
  //   // Fallthrough
  // case 11:
  //  minz = atof(argv[10]);
  //   // Fallthrough
  // case 10:
  //    std::istringstream(argv[9]) >> doFiltering;
  //   // Fallthrough
  case 9:
    messageInterval = atoi(argv[8]);
    // Fallthrough
  case 8:
    maxQ2 = atof(argv[7]);
    // Fallthrough
  case 7:
   minQ2 = atof(argv[6]);
    // Fallthrough
  case 6:
    pElectron = atof(argv[5]);
    // Fallthrough
  case 5:
    pProton = atof(argv[4]);
    // Fallthrough    
  case 4:
    pElectron = atof(argv[3]);
    // Fallthrough
  case 3:
    nEvents = atoi(argv[2]);
    // Fallthrough
  case 2:
    outFile = argv[1];
    break;
  default:
    cout << helpstring << endl;
    return -1;
  }
  cout<<endl;

  
  
  // runpythia(outFile, nEvents, pElectron, pProton, minQ2, maxQ2,
  // 	  messageInterval, doFiltering, minz, maxcorrelation);
  runpythia(outFile, nEvents, pElectron, pProton, minQ2, maxQ2, messageInterval);


  return 0;
}

  
