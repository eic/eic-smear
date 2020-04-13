// Load and smear a file,
// Then create a variety of QA plots
// This implementation mirrors the "standard" work flow
// Build -> Smear -> Load and befriend -> Analyze

#include "qaplots.hh"

#include <TSystem.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>

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
Smear::Detector BuildZeus();
Smear::Detector BuildBeAST();

#include <string>
#include <iomanip>
#include <cctype>
#include <exception>
#include <vector>
#include <map>

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::map;

// some helpers
const TString getrootname(const qaparameters& qapars );
void initializeqabook(const qaparameters& qapars, map<int,qacollection>& qabook );



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
  if ( qapars.detstring=="BEAST" ) detector = BuildBeAST();
  if ( qapars.detstring=="ZEUS" ) detector = BuildZeus();
  if ( qapars.detstring=="EPHENIX" ) detector = BuildEphoenix();
  // Did that work?
  if ( detector.GetNDevices() == 0 ) {
    cerr << "Detector sepcified as " << qapars.detstring
	 << " not recognized or empty." << endl;
    return -1;
  }
  TString smearedname = rootname;
  smearedname.ReplaceAll (".root",".smeared.root" );
  // Can disable warnings here.
  // erhic::DisKinematics::BoundaryWarning=false;
  SmearTree( detector, rootname.Data(), smearedname.Data());
  
  
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

  // We'll have a collection of histos and maybe other info for every pid
  // qacollection is defined in the header
  map<int,qacollection> qabook;
  // By default, use the five standard particles
  if ( qapars.pids.size() == 0 ) qapars.pids = { 11, 211, 321, 2212, 2112 }; // e, pi, K, p, n
  initializeqabook ( qapars, qabook );
    
  // Analysis loop
  // -------------  
  for(long iEvent=0; iEvent<inTree->GetEntries(); iEvent++){
    
    //Read next Event
    if(inTree->GetEntry(iEvent) <=0) break;
    if(iEvent%10000 == 0) cout << "Event " << iEvent << endl;

    // Loop over Particles
    for(int j=0; j<inEventS->GetNTracks(); j++){
      // Skip beam
      if ( j<3 ) continue;

      const Smear::ParticleMCS* inParticleS = inEventS->GetTrack(j); // Smeared Particle      
      const Particle* inParticle = inEvent->GetTrack(j); // Unsmeared Particle

      // Skip non-final particles. 
      if ( inParticle->GetStatus() != 1 ) continue;

      // Particle was not smeared
      if(inParticleS == NULL) continue; 

      // Fill histograms
      for ( auto& pidcoll : qabook ){
	auto& pid = pidcoll.first;
	auto& coll = pidcoll.second;
	if ( pid==0 || inParticle->GetPdgCode() == pid ){
	  auto delP = (inParticle->GetP() - inParticleS->GetP()) / inParticle->GetP();
	  coll.DelP_th->Fill(inParticle->GetTheta(), delP);
	  
	  auto delE = (inParticle->GetE() - inParticleS->GetE()) / inParticle->GetE();
	  coll.DelE_th->Fill(inParticle->GetTheta(), delE);

	  coll.dTh_p->Fill(inParticle->GetP(), inParticle->GetTheta() - inParticleS->GetTheta());
	  
	  coll.dPhi_p->Fill(inParticle->GetP(), inParticle->GetPhi() - inParticleS->GetPhi());	    
	}
      }
    }
  }


  // Make plots and save
  // -------------------
  new TCanvas;
  gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf[" );
  for ( auto& pidcoll : qabook ){
    auto& pid = pidcoll.first;
    auto& coll = pidcoll.second;

    coll.DelP_th->Draw("colz");
    gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );

    coll.DelE_th->Draw("colz");
    gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
    
    coll.dTh_p->Draw("colz");
    gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
    
    coll.dPhi_p->Draw("colz");
    gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
  }
  gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf]" );
  
  outfile->Write();
  return 0;
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
	 << " [-det detstring] beast, ephenix, zeus, handbook, jleic (capitalization does not matter.)" << endl
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
void initializeqabook(const qaparameters& qapars, map<int,qacollection>& qabook ){
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

  float demin = -0.1;
  float demax = 0.1;
  int debins = 100;

  float thmin = 0;
  float thmax = TMath::Pi();
  int thbins = 64;
  
  float dthmin = -0.1;
  float dthmax = 0.1;
  int dthbins = 100;

  float phimin = 0;
  float phimax = TMath::TwoPi();
  int phibins = 64;
  
  float dphimin = -0.1;
  float dphimax = 0.1;
  int dphibins = 100;

  for ( auto pid : qapars.pids ){
    pid = abs ( pid ); // ignoring charge

    s = qapars.detstring + "_DelP_th_"; s += pid;
    qabook[pid].DelP_th = new TH2D( s,s+";#theta;#Delta p/p", thbins, thmin, thmax, dpbins, dpmin, dpmax);

    s = qapars.detstring + "_DelE_th_"; s += pid;
    qabook[pid].DelE_th = new TH2D( s,s+";#theta;#Delta E/E", thbins, thmin, thmax, debins, demin, demax);
    
    s = qapars.detstring + "_dTh_p_"; s += pid;
    qabook[pid].dTh_p = new TH2D( s,s+";p;#Delta#theta", pbins, pmin, pmax, dthbins, dthmin, dthmax );

    s = qapars.detstring + "_dPhi_p_"; s += pid;
    qabook[pid].dPhi_p = new TH2D( s,s+";p;#Delta#phi", pbins, pmin, pmax, dphibins, dphimin, dphimax );    
  }

}
