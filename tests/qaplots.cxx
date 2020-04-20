// Load and smear a file,
// Then create a variety of QA plots
// This implementation mirrors the "standard" work flow
// Build -> Smear -> Load and befriend -> Analyze

#include "qaplots.hh"

#include <eicsmear/functions.h>
#include <eicsmear/smear/functions.h>

#include "eicsmear/erhic/EventBase.h"
#include "eicsmear/erhic/EventPythia.h"
#include "eicsmear/erhic/Particle.h"
#include "eicsmear/smear/EventSmear.h"
#include "eicsmear/smear/EventS.h"
#include "eicsmear/smear/Smear.h"
#include "eicsmear/smear/ParticleMCS.h"

#include "ePHENIXDetector.h"
#include "eicsmear/smear/Detector.h"
Smear::Detector BuildHandBookDetector();
Smear::Detector BuildPerfectDetector();
Smear::Detector BuildNaiveHandBookDetector();
Smear::Detector BuildZeus();
Smear::Detector BuildBeAST();

// Note: The remaining includes are not necessary for eic-smear usage
#include <TSystem.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TLegend.h>
#include <TText.h>
#include <TStyle.h>

#include <string>
#include <iomanip>
#include <cctype>
#include <exception>
#include <vector>
#include <map>

// Convenience only
using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::map;

// some helpers
static const TString getrootname(const qaparameters& qapars );
static void initializepidqabook(const qaparameters& qapars, map<int,pidqacollection>& qabook );
static void initializeeventqa(const qaparameters& qapars, eventqacollection& eventqa );
static void FillEventQA( eventqacollection& eventqa , const erhic::EventMC* const inEvent, const Smear::Event* const inEventS );
static void FillParticleQA( map<int,pidqacollection>& qabook, const Particle* const inParticle, const Smear::ParticleMCS* const inParticleS  );
static void PlotQA ( const qaparameters& qapars, eventqacollection& eventqa, map<int,pidqacollection>& qabook );
  
int main(int argc, char* argv[]){
  
  // Parse arguments
  // ---------------
  // Defaults are set in qaplots.h
  qaparameters qapars = ParseArguments ( argc, argv );

  // Set up output name
  TString rootname = getrootname(qapars);
  
  // First, convert txt file to tree
  // -------------------------------
  BuildTree(qapars.txtfilename.c_str(), qapars.outpath.c_str(), qapars.nevents);

  // Smear the tree
  // --------------
  Smear::Detector detector;
  if ( qapars.detstring=="HANDBOOK" ) detector = BuildHandBookDetector();
  if ( qapars.detstring=="PERFECT" ) detector = BuildPerfectDetector();
  if ( qapars.detstring=="NAIVEHANDBOOK" ) detector = BuildNaiveHandBookDetector();
  if ( qapars.detstring=="BEAST" ) detector = BuildBeAST();
  if ( qapars.detstring=="ZEUS" ) detector = BuildZeus();
  if ( qapars.detstring=="EPHENIX" ) detector = BuildEphoenix();
  if ( detector.GetNDevices() == 0 ) {
    cerr << "Detector sepcified as " << qapars.detstring
	 << " not recognized or empty." << endl;
    return -1;
  }
  TString smearedname = rootname;
  smearedname.ReplaceAll (".root",".smeared.root" );
  // Can disable warnings here.
  // Many warnings are harmless ( x or y can be smeared to values >1)
  // But it is recommended to leave it on and follow up on "inf", "nan" etc. if you test a new detector
  // erhic::DisKinematics::BoundaryWarning=false;
  SmearTree( detector, rootname.Data(), smearedname.Data());  

  // -------------
  // Load the tree
  // -------------
  TChain* inTree = new TChain("EICTree");
  inTree->Add(rootname);
  inTree->AddFriend("Smeared",smearedname);
  
  // Setup Input Event Buffer
  erhic::EventMC* inEvent(NULL);
  Smear::Event* inEventS(NULL);
  inTree->SetBranchAddress("event",&inEvent);
  inTree->SetBranchAddress("eventS",&inEventS);

  // Open histo root file and book histograms
  // ----------------------------------------
  TFile * outfile = new TFile ( qapars.outfilebase + qapars.detstring + ".root", "RECREATE");

  // A collection of event-wise qa plots, like kinematics
  // eventqacollection is defined in the header
  eventqacollection eventqa = {}; // this syntax initializes everything to 0
  initializeeventqa ( qapars, eventqa );
  
  // We'll also have a collection of histos and maybe other info for every pid
  // pidqacollection is defined in the header
  map<int,pidqacollection> qabook;
  // By default, use the standard particles
  if ( qapars.pids.size() == 0 ) qapars.pids = { 11, 211, 321, 2212, 2112 }; // e, pi, K, p, n
  initializepidqabook ( qapars, qabook );
    
  // -------------
  // Analysis loop
  // -------------  
  for(long iEvent=0; iEvent<inTree->GetEntries(); iEvent++){
    
    //Read next Event
    if(inTree->GetEntry(iEvent) <=0) break;
    if(iEvent%10000 == 0) cout << "Event " << iEvent << endl;
    
    // -------------
    // event-wise QA
    // -------------
    // Following function just contains many statements of the form
    // if ( inEventS->GetQ2()>0 ) eventqa.Q2_NM->Fill ( std::log10(inEvent->GetQ2()), std::log10(inEventS->GetQ2()));
    FillEventQA( eventqa, inEvent, inEventS );

    // -------------------
    // Loop over Particles
    // -------------------
    for(int j=0; j<inEventS->GetNTracks(); j++){
      // Skip beam
      if ( j<3 ) continue;

      const Smear::ParticleMCS* inParticleS = inEventS->GetTrack(j); // Smeared Particle      
      const Particle* inParticle = inEvent->GetTrack(j); // Unsmeared Particle

      // Skip non-final particles. 
      if ( inParticle->GetStatus() != 1 ) continue;

      // Particle was not smeared
      if(inParticleS == NULL) continue; 

      // ----------------
      // Particle-wise QA
      // ----------------
      // Following function just contains many statements of the form
      // coll.dEta_p->Fill(inParticle->GetP(), inParticle->GetEta() - inParticleS->GetEta());
      FillParticleQA( qabook, inParticle, inParticleS );      
    }
  }
  
  // NOTE: The remainder of this long file is tedious and explicit creation and filling of histograms
  // The Fill*QA functions demonstrate how to access values in detail, but this is the end of the 
  // basic work flow to use eic-smear from start to finish! 

  qapars.usedevents = inTree->GetEntries();
  PlotQA( qapars, eventqa, qabook );  
  outfile->Write();
  
  return 0;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

void FillEventQA( eventqacollection& eventqa , const erhic::EventMC* const inEvent, const Smear::Event* const inEventS ){

    // Q2
    if ( eventqa.Q2_NM && inEventS->GetQ2()>0)                eventqa.Q2_NM->Fill ( std::log10(inEvent->GetQ2()), std::log10(inEventS->GetQ2()));
    else eventqa.missedQ2_NM++;
    if ( eventqa.Q2_JB && inEventS->GetQ2JacquetBlondel()>0 ) eventqa.Q2_JB->Fill ( std::log10(inEvent->GetQ2()), std::log10(inEventS->GetQ2JacquetBlondel()));
    else eventqa.missedQ2_JB++;
    if ( eventqa.Q2_DA && inEventS->GetQ2DoubleAngle()>0 )    eventqa.Q2_DA->Fill ( std::log10(inEvent->GetQ2()), std::log10(inEventS->GetQ2DoubleAngle()));
    else eventqa.missedQ2_DA++;
    
    if ( eventqa.delQ2_NM && inEventS->GetQ2()>0)                eventqa.delQ2_NM->Fill ( std::log10(inEvent->GetQ2()), std::abs(inEventS->GetQ2()-inEvent->GetQ2())/ inEvent->GetQ2());
    if ( eventqa.delQ2_JB && inEventS->GetQ2JacquetBlondel()>0 ) eventqa.delQ2_JB->Fill ( std::log10(inEvent->GetQ2()), std::abs(inEventS->GetQ2JacquetBlondel()-inEvent->GetQ2())/ inEvent->GetQ2());
    if ( eventqa.delQ2_DA && inEventS->GetQ2DoubleAngle()>0 )    eventqa.delQ2_DA->Fill ( std::log10(inEvent->GetQ2()), (inEventS->GetQ2DoubleAngle()-inEvent->GetQ2())/ inEvent->GetQ2());
    
    // y
    if ( eventqa.y_NM && inEventS->GetY()>0)                eventqa.y_NM->Fill ( inEvent->GetY(), inEventS->GetY());
    else eventqa.missedy_NM++;
    if ( eventqa.y_JB && inEventS->GetYJacquetBlondel()>0 ) eventqa.y_JB->Fill ( inEvent->GetY(), inEventS->GetYJacquetBlondel());
    else eventqa.missedy_JB++;
    if ( eventqa.y_DA && inEventS->GetYDoubleAngle()>0 )    eventqa.y_DA->Fill ( inEvent->GetY(), inEventS->GetYDoubleAngle());
    else eventqa.missedy_DA++;
    
    if ( eventqa.dely_NM && inEventS->GetY()>0)                eventqa.dely_NM->Fill ( inEvent->GetY(), std::abs(inEventS->GetY()-inEvent->GetY())/ inEvent->GetY());
    if ( eventqa.dely_JB && inEventS->GetYJacquetBlondel()>0 ) eventqa.dely_JB->Fill ( inEvent->GetY(), std::abs(inEventS->GetYJacquetBlondel()-inEvent->GetY())/ inEvent->GetY());
    if ( eventqa.dely_DA && inEventS->GetYDoubleAngle()>0 )    eventqa.dely_DA->Fill ( inEvent->GetY(), std::abs(inEventS->GetYDoubleAngle()-inEvent->GetY())/ inEvent->GetY());

    // if ( eventqa.y_NM && inEventS->GetY()>0)                eventqa.y_NM->Fill ( std::log10(inEvent->GetY()), std::log10(inEventS->GetY()));
    // else eventqa.missedy_NM++;
    // if ( eventqa.y_JB && inEventS->GetYJacquetBlondel()>0 ) eventqa.y_JB->Fill ( std::log10(inEvent->GetY()), std::log10(inEventS->GetYJacquetBlondel()));
    // else eventqa.missedy_JB++;
    // if ( eventqa.y_DA && inEventS->GetYDoubleAngle()>0 )    eventqa.y_DA->Fill ( std::log10(inEvent->GetY()), std::log10(inEventS->GetYDoubleAngle()));
    // else eventqa.missedy_DA++;
    
    // if ( eventqa.dely_NM && inEventS->GetY()>0)                eventqa.dely_NM->Fill ( std::log10(inEvent->GetY()), std::abs(inEventS->GetY()-inEvent->GetY())/ inEvent->GetY());
    // if ( eventqa.dely_JB && inEventS->GetYJacquetBlondel()>0 ) eventqa.dely_JB->Fill ( std::log10(inEvent->GetY()), std::abs(inEventS->GetYJacquetBlondel()-inEvent->GetY())/ inEvent->GetY());
    // if ( eventqa.dely_DA && inEventS->GetYDoubleAngle()>0 )    eventqa.dely_DA->Fill ( std::log10(inEvent->GetY()), std::abs(inEventS->GetYDoubleAngle()-inEvent->GetY())/ inEvent->GetY());

    // x
    if ( eventqa.x_NM && inEventS->GetX()>0)                eventqa.x_NM->Fill ( std::log10(inEvent->GetX()), std::log10(inEventS->GetX()));
    else eventqa.missedx_NM++;
    if ( eventqa.x_JB && inEventS->GetXJacquetBlondel()>0 ) eventqa.x_JB->Fill ( std::log10(inEvent->GetX()), std::log10(inEventS->GetXJacquetBlondel()));
    else eventqa.missedx_JB++;
    if ( eventqa.x_DA && inEventS->GetXDoubleAngle()>0 )    eventqa.x_DA->Fill ( std::log10(inEvent->GetX()), std::log10(inEventS->GetXDoubleAngle()));
    else eventqa.missedx_DA++;
    
    if ( eventqa.delx_NM && inEventS->GetX()>0)                eventqa.delx_NM->Fill ( std::log10(inEvent->GetX()), std::abs(inEventS->GetX()-inEvent->GetX())/ inEvent->GetX());
    if ( eventqa.delx_JB && inEventS->GetXJacquetBlondel()>0 ) eventqa.delx_JB->Fill ( std::log10(inEvent->GetX()), std::abs(inEventS->GetXJacquetBlondel()-inEvent->GetX())/ inEvent->GetX());
    if ( eventqa.delx_DA && inEventS->GetXDoubleAngle()>0 )    eventqa.delx_DA->Fill ( std::log10(inEvent->GetX()), std::abs(inEventS->GetXDoubleAngle()-inEvent->GetX())/ inEvent->GetX());

}
// ---------------------------------------------------------------
void FillParticleQA( map<int,pidqacollection>& qabook, const Particle* const inParticle, const Smear::ParticleMCS* const inParticleS  ){
    
  // If any component is smeared, all others are either smeared or 0 (meaning "not detected").
  // Could detect the latter case (with a given accuracy):
  // const double epsilon = 1e-9;
  
  // Fill histograms
  for ( auto& pidcoll : qabook ){
    auto& pid = pidcoll.first;
    auto& coll = pidcoll.second;
    if ( pid==0 || inParticle->GetPdgCode() == pid ){
      
      // if ( std::abs(inParticleS->GetP()) > epsilon ){
      //   auto delP = (inParticle->GetP() - inParticleS->GetP()) / inParticle->GetP();
      //   coll.DelP_th->Fill(inParticle->GetTheta(), delP);
      // }
      
      // if ( std::abs(inParticleS->GetE()) > epsilon ){
      //   auto delE = (inParticle->GetE() - inParticleS->GetE()) / inParticle->GetE();
      //   coll.DelE_th->Fill(inParticle->GetTheta(), delE);
      // }
      
      // if ( std::abs(inParticleS->GetTheta()) > epsilon ){
      //   coll.dTh_p->Fill(inParticle->GetP(), inParticle->GetTheta() - inParticleS->GetTheta());
      // }
      
      // if ( std::abs(inParticleS->GetPhi()) > epsilon ){
      //   coll.dPhi_p->Fill(inParticle->GetP(), inParticle->GetPhi() - inParticleS->GetPhi());
      // }
      
      auto delP = (inParticle->GetP() - inParticleS->GetP()) / inParticle->GetP();
      coll.DelP_th->Fill(inParticle->GetTheta(), delP);
      coll.DelP_eta->Fill(inParticle->GetEta(), delP);
      
      auto delE = (inParticle->GetE() - inParticleS->GetE()) / inParticle->GetE();
      coll.DelE_E->Fill(inParticle->GetE(), delE);
      coll.DelE_th->Fill(inParticle->GetTheta(), delE);
      coll.DelE_eta->Fill(inParticle->GetEta(), delE);
      
      coll.dTh_p->Fill(inParticle->GetP(), inParticle->GetTheta() - inParticleS->GetTheta());
      coll.dEta_p->Fill(inParticle->GetP(), inParticle->GetEta() - inParticleS->GetEta());
      coll.dPhi_p->Fill(inParticle->GetP(), inParticle->GetPhi() - inParticleS->GetPhi());
    }
  }
}

// -----------------------------------------------------------------------------

qaparameters ParseArguments ( int argc, char* argv[] ){
  vector<string> arguments(argv + 1, argv + argc);
  bool argsokay=true;
  qaparameters qapars;
  try{
    for ( auto parg = arguments.begin() ; parg!=arguments.end() ; ++parg){
      string arg=*parg;
      if ( arg == "-h" ){
	argsokay=false;
	break; 
      } else if ( arg == "-o" ){
	if (++parg == arguments.end() ){ argsokay=false; break; }
	qapars.outfilebase=*parg;
      } else if ( arg == "-i" ){
	if (++parg ==arguments.end() ){ argsokay=false; break; }
	qapars.txtfilename=*parg;
      } else if ( arg == "-N" ){      
	if (++parg==arguments.end() ){ argsokay=false; break; }
	qapars.nevents=std::stoi(parg->data());
      } else if ( arg == "-addpid" ){
	if ( ++parg == arguments.end() ){ argsokay=false; break; }
	qapars.pids.push_back(std::stoi(parg->data()));
      } else if ( arg == "-det" ){
	if (++parg == arguments.end() ){ argsokay=false; break; }
	qapars.detstring=*parg;
      } else {
	argsokay=false;
	break;
      }
    }
  } catch ( const std::exception& e){
    cerr << "Caught exception during argument parsing: "
	 << e.what() << endl;
    argsokay=false; 
  }  
  
  if ( !argsokay ) {
    cerr << "usage: " << argv[0] << endl
      	 << " [-i txtfilename] (Lund-style file)"  << endl
	 << " [-o OutFileBase] (extension will be added)"  << endl
      	 << " [-N Nevents] (<0 for all)" << endl
      	 << " [-addpid pid] (can be called multiple times)" << endl
	 << " [-det detstring] handbook, perfect, beast, ephenix, zeus, jleic, naivehandbook (capitalization does not matter.)" << endl
	 << endl;
    throw std::runtime_error("Not a valid list of options");
  }
  for (auto & c: qapars.detstring) c = toupper(c);
  
  return qapars;
}

// ---------------------------------------------------------------
const TString getrootname(const qaparameters& qapars ){
  // The root file name is created by replacing .txt by .root
  // Make sure input name is right and generate output name for loading
  if ( qapars.txtfilename.substr(qapars.txtfilename.length()-4,4) != ".txt" ){
    cerr << "Input file " << qapars.txtfilename << " doesn't end with .txt";
    throw std::runtime_error("Can't parse input file");
  }
  TString rootname = qapars.txtfilename.substr(0,qapars.txtfilename.length()-4);
  // BuildTree includes event number in partial transformation
  if ( qapars.nevents>=0 ) {
    rootname += ".";
    rootname += qapars.nevents;
    rootname += "event";
  }
  rootname += ".root";
  rootname = gSystem->BaseName( rootname );
  rootname.Prepend( qapars.outpath );
  cout << " ======================= " << endl;
  cout << " Transforming input file " << endl
       << qapars.txtfilename << endl
       << " into root file " << endl
       << rootname << endl;
  return rootname;
}
// ---------------------------------------------------------------
void initializepidqabook(const qaparameters& qapars, map<int,pidqacollection>& qabook ){
  gStyle->SetHistLineColor(kRed); // for Profiles

  TString s;
  float pmin = 0;
  float pmax = 20;
  int pbins = 80;

  float dpmin = -0.1;
  float dpmax = 0.1;
  int dpbins = 100;
  
  float emin = 0;
  float emax = 20;
  int ebins = 80;

  float demin = -1;
  float demax = 1;
  int debins = 100;

  float thmin = 0;
  float thmax = TMath::Pi();
  int thbins = 64;
  
  float dthmin = -0.1;
  float dthmax = 0.1;
  int dthbins = 100;

  float etamin = -5;
  float etamax = 5;
  int etabins = 100;
  
  float detamin = -0.1;
  float detamax = 0.1;
  int detabins = 100;

  float phimin = 0;
  float phimax = TMath::TwoPi();
  int phibins = 64;
  
  float dphimin = -0.1;
  float dphimax = 0.1;
  int dphibins = 100;

  for ( auto pid : qapars.pids ){
    pid = abs ( pid ); // ignoring charge

    s = qapars.detstring + "_DelE_E_"; s += pid;
    qabook[pid].DelE_E = new TH2D( s,s+";E;#Delta E/E", ebins, emin, emax, debins, demin, demax);

    s = qapars.detstring + "_dPhi_p_"; s += pid;
    qabook[pid].dPhi_p = new TH2D( s,s+";p;#Delta#phi", pbins, pmin, pmax, dphibins, dphimin, dphimax );    

    s = qapars.detstring + "_DelP_th_"; s += pid;
    qabook[pid].DelP_th = new TH2D( s,s+";#theta;#Delta p/p", thbins, thmin, thmax, dpbins, dpmin, dpmax);

    s = qapars.detstring + "_DelE_th_"; s += pid;
    qabook[pid].DelE_th = new TH2D( s,s+";#theta;#Delta E/E", thbins, thmin, thmax, debins, demin, demax);
    
    s = qapars.detstring + "_dTh_p_"; s += pid;
    qabook[pid].dTh_p = new TH2D( s,s+";p;#Delta#theta", pbins, pmin, pmax, dthbins, dthmin, dthmax );

    s = qapars.detstring + "_DelP_eta_"; s += pid;
    qabook[pid].DelP_eta = new TH2D( s,s+";#eta;#Delta p/p", etabins, etamin, etamax, dpbins, dpmin, dpmax);

    s = qapars.detstring + "_DelE_eta_"; s += pid;
    qabook[pid].DelE_eta = new TH2D( s,s+";#eta;#Delta E/E", etabins, etamin, etamax, debins, demin, demax);
    
    s = qapars.detstring + "_dEta_p_"; s += pid;
    qabook[pid].dEta_p = new TH2D( s,s+";p;#Delta#eta", pbins, pmin, pmax, detabins, detamin, detamax );

  }

}

// ---------------------------------------------------------------
void initializeeventqa(const qaparameters& qapars, eventqacollection& eventqa ){
  // recording log of x, y, Q2
  
  TString s;

  float Q2min = 1e-3;
  float Q2max = 1e5;
  int logQ2bins = 240;

  s = qapars.detstring + "_LogQ2_NM";
  eventqa.Q2_NM = new TH2D( s,s+";log Q^{2};log Q^{2}_{NM}", logQ2bins, std::log10(Q2min), std::log10(Q2max), logQ2bins, std::log10(Q2min), std::log10(Q2max));

  s = qapars.detstring + "_LogQ2_JB";
  eventqa.Q2_JB = new TH2D( s,s+";log Q^{2};log Q^{2}_{JB}", logQ2bins, std::log10(Q2min), std::log10(Q2max), logQ2bins, std::log10(Q2min), std::log10(Q2max));

  s = qapars.detstring + "_LogQ2_DA";
  eventqa.Q2_DA = new TH2D( s,s+";log Q^{2};log Q^{2}_{DA}", logQ2bins, std::log10(Q2min), std::log10(Q2max), logQ2bins, std::log10(Q2min), std::log10(Q2max));

  int delQ2bins = 100;
  s = qapars.detstring + "_delQ2_NM";
  eventqa.delQ2_NM = new TH2D( s,s+";log Q^{2};(Q^{2}_{NM}-Q^{2})/Q2", logQ2bins, std::log10(Q2min), std::log10(Q2max), delQ2bins, 0, 1.2);

  s = qapars.detstring + "_delQ2_JB";
  eventqa.delQ2_JB = new TH2D( s,s+";log Q^{2};(Q^{2}_{JB}-Q^{2})/Q2", logQ2bins, std::log10(Q2min), std::log10(Q2max), delQ2bins, 0, 1.2);

  s = qapars.detstring + "_delQ2_DA";
  eventqa.delQ2_DA = new TH2D( s,s+";log Q^{2};(Q^{2}_{DA}-Q^{2})/Q2", logQ2bins, std::log10(Q2min), std::log10(Q2max), delQ2bins, 0, 1.2);


  float ymin = 0;
  float ymax = 1.2; // see y>1 as well
  int ybins = 100;
  
  s = qapars.detstring + "_y_NM";
  eventqa.y_NM = new TH2D( s,s+";y;y_{NM}", ybins, ymin, ymax, ybins, ymin, ymax);

  s = qapars.detstring + "_y_JB";
  eventqa.y_JB = new TH2D( s,s+";y;y_{JB}", ybins, ymin, ymax, ybins, ymin, ymax);

  s = qapars.detstring + "_y_DA";
  eventqa.y_DA = new TH2D( s,s+";y;y_{DA}", ybins, ymin, ymax, ybins, ymin, ymax);

  int delybins = 100;
  s = qapars.detstring + "_dely_NM";
  eventqa.dely_NM = new TH2D( s,s+";y;(y_{NM}-y)/y", ybins, ymin, ymax, delybins, 0, 1.2);

  s = qapars.detstring + "_dely_JB";
  eventqa.dely_JB = new TH2D( s,s+";y;(y_{JB}-y)/y", ybins, ymin, ymax, delybins, 0, 1.2);

  s = qapars.detstring + "_dely_DA";
  eventqa.dely_DA = new TH2D( s,s+";y;(y_{DA}-y)/y", ybins, ymin, ymax, delybins, 0, 1.2);

  // float ymin = 1e-4;
  // float ymax = 1e1; // see y>1 as well
  // int logybins = 200;
  // s = qapars.detstring + "_Logy_NM";
  // eventqa.y_NM = new TH2D( s,s+";log y;log y_{NM}", logybins, std::log10(ymin), std::log10(ymax), logybins, std::log10(ymin), std::log10(ymax));

  // s = qapars.detstring + "_Logy_JB";
  // eventqa.y_JB = new TH2D( s,s+";log y;log y_{JB}", logybins, std::log10(ymin), std::log10(ymax), logybins, std::log10(ymin), std::log10(ymax));

  // s = qapars.detstring + "_Logy_DA";
  // eventqa.y_DA = new TH2D( s,s+";log y;log y_{DA}", logybins, std::log10(ymin), std::log10(ymax), logybins, std::log10(ymin), std::log10(ymax));

  // int delybins = 100;
  // s = qapars.detstring + "_dely_NM";
  // eventqa.dely_NM = new TH2D( s,s+";log y;(y_{NM}-y)/y", logybins, std::log10(ymin), std::log10(ymax), delybins, 0, 1.2);

  // s = qapars.detstring + "_dely_JB";
  // eventqa.dely_JB = new TH2D( s,s+";log y;(y_{JB}-y)/y", logybins, std::log10(ymin), std::log10(ymax), delybins, 0, 1.2);

  // s = qapars.detstring + "_dely_DA";
  // eventqa.dely_DA = new TH2D( s,s+";log y;(y_{DA}-y)/y", logybins, std::log10(ymin), std::log10(ymax), delybins, 0, 1.2);

  float xmin = 1e-4;
  float xmax = 1e1; // see x>1 as well
  int logxbins = 200;

  s = qapars.detstring + "_Logx_NM";
  eventqa.x_NM = new TH2D( s,s+";log x;log x_{NM}", logxbins, std::log10(xmin), std::log10(xmax), logxbins, std::log10(xmin), std::log10(xmax));

  s = qapars.detstring + "_Logx_JB";
  eventqa.x_JB = new TH2D( s,s+";log x;log x_{JB}", logxbins, std::log10(xmin), std::log10(xmax), logxbins, std::log10(xmin), std::log10(xmax));

  s = qapars.detstring + "_Logx_DA";
  eventqa.x_DA = new TH2D( s,s+";log x;log x_{DA}", logxbins, std::log10(xmin), std::log10(xmax), logxbins, std::log10(xmin), std::log10(xmax));

  int delxbins = 100;
  s = qapars.detstring + "_delx_NM";
  eventqa.delx_NM = new TH2D( s,s+";log x;(x_{NM}-x)/x", logxbins, std::log10(xmin), std::log10(xmax), delxbins, 0, 1.2);

  s = qapars.detstring + "_delx_JB";
  eventqa.delx_JB = new TH2D( s,s+";log x;(x_{JB}-x)/x", logxbins, std::log10(xmin), std::log10(xmax), delxbins, 0, 1.2);

  s = qapars.detstring + "_delx_DA";
  eventqa.delx_DA = new TH2D( s,s+";log x;(x_{DA}-x)/x", logxbins, std::log10(xmin), std::log10(xmax), delxbins, 0, 1.2);
}

// ---------------------------------------------------------------
void PlotQA ( const qaparameters& qapars, eventqacollection& eventqa, map<int,pidqacollection>& qabook ){

  // Stat  position and size
  // -----------------------
  gStyle->SetStatX(0.25);
  gStyle->SetStatW(0.15);
  gStyle->SetStatY(0.9);
  gStyle->SetStatH(0.15);

  // Position of the "Missed: " box
  float missx = 0.55;
  float missy = 0.2;
  float missy2 = 0.8;
  TText t;
  t.SetNDC();

  // prep a pdf collection
  new TCanvas;
  gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf[" );

  // event-wise qa

  // response-style
  // NM
  if ( eventqa.y_NM ) {
    eventqa.y_NM->Draw("colz");
    t.DrawText( missx,missy, Form("Missed: %ld / %ld",eventqa.missedy_NM, qapars.usedevents));
    gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
  }
  if ( eventqa.x_NM ) {
    eventqa.x_NM->Draw("colz");
    t.DrawText( missx,missy, Form("Missed: %ld / %ld",eventqa.missedx_NM, qapars.usedevents));
    gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
  }
  if ( eventqa.Q2_NM ) {
    eventqa.Q2_NM->Draw("colz");
    t.DrawText( missx,missy, Form("Missed: %ld / %ld",eventqa.missedQ2_NM, qapars.usedevents));
    gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
  }
  // DA
  if ( eventqa.y_DA ) {
    eventqa.y_DA->Draw("colz");
    t.DrawText( missx,missy, Form("Missed: %ld / %ld",eventqa.missedy_DA, qapars.usedevents));
    gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
  }
  if ( eventqa.x_DA ) {
    eventqa.x_DA->Draw("colz");
    t.DrawText( missx,missy, Form("Missed: %ld / %ld",eventqa.missedx_DA, qapars.usedevents));
    gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
  }
  if ( eventqa.Q2_DA ) {
    eventqa.Q2_DA->Draw("colz");
    t.DrawText( missx,missy, Form("Missed: %ld / %ld",eventqa.missedQ2_DA, qapars.usedevents));
    gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
  }
  // JB
  if ( eventqa.y_JB ) {
    eventqa.y_JB->Draw("colz");
    t.DrawText( missx,missy, Form("Missed: %ld / %ld",eventqa.missedy_JB, qapars.usedevents));
    gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
  }
  if ( eventqa.x_JB ) {
    eventqa.x_JB->Draw("colz");
    t.DrawText( missx,missy, Form("Missed: %ld / %ld",eventqa.missedx_JB, qapars.usedevents));
    gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
  }
  if ( eventqa.Q2_JB ) {
    eventqa.Q2_JB->Draw("colz");
    t.DrawText( missx,missy, Form("Missed: %ld / %ld",eventqa.missedQ2_JB, qapars.usedevents));
    gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
  }

  // resolution-style
  // NM
  if ( eventqa.dely_NM ) {
    eventqa.dely_NM->Draw("colz");
    t.DrawText( missx,missy2, Form("Missed: %ld / %ld",eventqa.missedy_NM, qapars.usedevents));
    gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
  }
  if ( eventqa.delx_NM ) {
    eventqa.delx_NM->Draw("colz");
    t.DrawText( missx,missy2, Form("Missed: %ld / %ld",eventqa.missedx_NM, qapars.usedevents));
    gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
  }
  if ( eventqa.delQ2_NM ) {
    eventqa.delQ2_NM->Draw("colz");
    t.DrawText( missx,missy2, Form("Missed: %ld / %ld",eventqa.missedQ2_NM, qapars.usedevents));
    gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
  }
  // DA
  if ( eventqa.dely_DA ) {
    eventqa.dely_DA->Draw("colz");
    t.DrawText( missx,missy2, Form("Missed: %ld / %ld",eventqa.missedy_DA, qapars.usedevents));
    gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
  }
  if ( eventqa.delx_DA ) {
    eventqa.delx_DA->Draw("colz");
    t.DrawText( missx,missy2, Form("Missed: %ld / %ld",eventqa.missedx_DA, qapars.usedevents));
    gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
  }
  if ( eventqa.delQ2_DA ) {
    eventqa.delQ2_DA->Draw("colz");
    t.DrawText( missx,missy2, Form("Missed: %ld / %ld",eventqa.missedQ2_DA, qapars.usedevents));
    gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
  }
  // JB
  if ( eventqa.dely_JB ) {
    eventqa.dely_JB->Draw("colz");
    t.DrawText( missx,missy2, Form("Missed: %ld / %ld",eventqa.missedy_JB, qapars.usedevents));
    gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
  }
  if ( eventqa.delx_JB ) {
    eventqa.delx_JB->Draw("colz");
    t.DrawText( missx,missy2, Form("Missed: %ld / %ld",eventqa.missedx_JB, qapars.usedevents));
    gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
  }
  if ( eventqa.delQ2_JB ) {
    eventqa.delQ2_JB->Draw("colz");
    t.DrawText( missx,missy2, Form("Missed: %ld / %ld",eventqa.missedQ2_JB, qapars.usedevents));
    gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
  }

  

  // particle QA
  // -----------
  gStyle->SetStatX(0.55); // reposition stat box
  for ( auto& pidcoll : qabook ){
    auto& pid = pidcoll.first;
    auto& coll = pidcoll.second;

    // option "s" in Profile shows rms
    
    // coll.DelP_th->Draw("colz");
    // gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );

    coll.DelP_eta->Draw("colz");
    coll.DelP_eta->ProfileX("_px",1,-1,"s")->Draw("same");
    gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );

    // coll.DelE_th->Draw("colz");
    // gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
    
    coll.DelE_eta->Draw("colz");
    coll.DelE_eta->ProfileX("_px",1,-1,"s")->Draw("same");
    gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );

    coll.DelE_E->Draw("colz");
    coll.DelE_E->ProfileX("_px",1,-1,"s")->Draw("same");
    gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );

    // coll.dTh_p->Draw("colz");
    // gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
    
    coll.dEta_p->Draw("colz");
    coll.dEta_p->ProfileX("_px",1,-1,"s")->Draw("same");
    gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );

    coll.dPhi_p->Draw("colz");
    coll.dPhi_p->ProfileX("_px",1,-1,"s")->Draw("same");
    gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
  }
  // close the pdf collection
  gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf]" );
}
// ---------------------------------------------------------------
