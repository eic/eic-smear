#ifndef QAPLOTS_H
#define QAPLOTS_H

// Note: Nothing in here is necessary for eic-smear usage

#include <TString.h>
#include <TH1.h>
#include <TH2.h>

#include <string>
#include <vector>

struct qaparameters{
  std::string txtfilename="./tests/ep_hiQ2.20x250.small.txt";
  TString outfilebase="./qaplots";
  std::string outpath="./";
  long nevents=-1;
  std::vector<int> pids = {}; // sign will be ignored. 0 for all. leave empty for e, pi, k, p. 
  std::string detstring = "HANDBOOK"; // Capitalization does not matter

  long usedevents=-1;// pure convenience so I can access the true number when nevents=-1;
};

struct eventqacollection {
  TH2D* Q2_NM;
  TH2D* Q2_JB;
  TH2D* Q2_DA;
  TH2D* y_NM;
  TH2D* y_JB;
  TH2D* y_DA;
  TH2D* x_NM;
  TH2D* x_JB;
  TH2D* x_DA;

  long missedQ2_NM;
  long missedQ2_JB;
  long missedQ2_DA;
  long missedy_NM;
  long missedy_JB;
  long missedy_DA;
  long missedx_NM;
  long missedx_JB;
  long missedx_DA;

  TH2D* delQ2_NM;
  TH2D* delQ2_JB;
  TH2D* delQ2_DA;
  TH2D* dely_NM;
  TH2D* dely_JB;
  TH2D* dely_DA;
  TH2D* delx_NM;
  TH2D* delx_JB;
  TH2D* delx_DA;

};

struct pidqacollection {
  TH2D* DelE_E;
  TH2D* dPhi_p;
  TH2D* DelP_th;
  TH2D* DelE_th;
  TH2D* dTh_p;
  TH2D* DelP_eta;
  TH2D* DelE_eta;
  TH2D* dEta_p;
};



qaparameters ParseArguments ( int argc, char* argv[] );

#endif // QAPLOTS_H
