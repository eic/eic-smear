#ifndef QAPLOTS_H
#define QAPLOTS_H

#include <TString.h>

#include <string>
#include <vector>

struct qaparameters{
  // std::string txtfilename="./tests/beagle_eD.txt";
  std::string txtfilename="./tests/ep_noradcorr.20x250.small.txt";
  TString outfilebase="./qaplots";
  long nevents=-1;
  std::vector<int> pids = {}; // sign will be ignored. 0 for all. leave empty for e, pi, k, p. 
  std::string detstring = "BeAST"; // Capitalization does not matter
};

qaparameters ParseArguments ( int argc, char* argv[] );

#endif // QAPLOTS_H
