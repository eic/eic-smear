#ifndef QAPLOTS_H
#define QAPLOTS_H

#include <TString.h>
#include <TH1.h>
#include <TH2.h>

#include <string>
#include <vector>

struct qaparameters{
  // std::string txtfilename="./tests/beagle_eD.txt";
  std::string txtfilename="./tests/ep_noradcorr.20x250.small.txt";
  TString outfilebase="./qaplots";
  std::string outpath="./";
  long nevents=-1;
  std::vector<int> pids = {}; // sign will be ignored. 0 for all. leave empty for e, pi, k, p. 
  std::string detstring = "BeAST"; // Capitalization does not matter
};

struct qacollection {
  TH2D* DelP_th;
  TH2D* DelE_th;
  TH2D* dTh_p;
  TH2D* dPhi_p;
};



qaparameters ParseArguments ( int argc, char* argv[] );

#endif // QAPLOTS_H
