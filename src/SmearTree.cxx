//here lies the function for smearing trees

#include <TROOT.h>
#include <TSystem.h>
#include <TChain.h>
#include <TString.h>
#include <TRandom2.h>
#include <TTree.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <TStopwatch.h>
#include <TH1D.h>

#include "EventBase.h"
#include "VirtualParticle.h"

#include "Smear.h"
#include "Detector.h"
#include "Device.h"
#include "ParticleID.h"
#include "SmearEvent.h"

//Long64_t ParticleS::count(0);
//Long64_t EventS::count(0);

using namespace std;

/**
 Smear a tree using the Detector det.  If you do not specify the number of events to be smeared, this function
 will automatically smear all events in the tree.
 */
int SmearTree(Smear::Detector det, TString inFileName, TString outFileName ="SAME", Long64_t nEvents=-1) {

	TFile inFile(inFileName,"READ");
	TTree *tree = NULL;
	
	inFile.GetObject("EICTree",tree);
	
	TString fileNames;
	
	if (outFileName.Contains("SAME")) {
		fileNames = inFileName.ReplaceAll(".root",".smear.root");
	} else {
		fileNames = outFileName;
	}
	
	TFile file(fileNames,"RECREATE");
	
	if (!(inFile.IsOpen()) || !(file.IsOpen())) {
		std::cout << "! Error ! File not found. \n";
		return 0;
	}
	
	TTree stre("Smeared","A tree with numerical fruit.");

	EventBase *event = NULL;
	EventS *eventS = NULL;

  tree->SetBranchAddress("event",&event);
	stre.Branch("eventS",&eventS);

   TH1D cpuTimePerEvent("cpuTimePerEvent", "", 200, -5., -4.);
   TH1D realTimePerEvent("realTimePerEvent", "", 200, -5., -4.);
   cpuTimePerEvent.SetBit(TH1::kCanRebin);
   realTimePerEvent.SetBit(TH1::kCanRebin);
   TStopwatch timer;
   
	if (tree->GetEntries() < nEvents || -1 == nEvents ) {
		nEvents = tree->GetEntries();
	}
	
	std::cout << "/-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-/" << std::endl;
	std::cout << "/  Commencing Smearing of " << nEvents << " events." << std::endl;
	std::cout << "/-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-/" << std::endl;
	
	for (Long64_t i(0); i<nEvents; i++) {
		timer.Start();
		if (i % 10000 == 0 && i != 0) {
			std::cout << "Processing event " << i << "." << std::endl;
//         std::cout << "# of ParticleS: " << ParticleS::count << std::endl;
//         std::cout << "# of EventS:    " << EventS::count << std::endl;
		}
		
		tree->GetEntry(i);
		
		int nParticles = event->NTracks();
		
		for (int j=0; j<nParticles; j++) {
         Particle* ptr = event->GetTrack(j);
         if(not ptr) {
            continue;
         } // if
         
         Particle& particle = *ptr;
			
			eventS->smearedParticles.push_back(det.DetSmear(particle) );
		} // for
		det.FillEventKinematics(event,eventS);
		
		stre.Fill();
		
		eventS->Reset();
      // RealTime() also stops the timer.
		cpuTimePerEvent.Fill(log10(timer.CpuTime()));
      realTimePerEvent.Fill(log10(timer.RealTime()));
	}
	
	file.cd();
   cpuTimePerEvent.Write();
   realTimePerEvent.Write();
	stre.Write();
   file.Purge();
	
	std::cout << "|~~~~~~~~~~~~~~~~~~ Completed Successfully ~~~~~~~~~~~~~~~~~~~|" << std::endl;
	return 0;
}