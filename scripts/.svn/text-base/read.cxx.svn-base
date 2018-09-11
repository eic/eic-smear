// read.cxx
//
// Created by TB on 6/13/11.
// Copyright 2011 BNL. All rights reserved.
//
// Example of how to read a file produced by BuildTree for a simple analysis.
// To run, in ROOT do:
// root [0] .L /path/to/read.cxx
// root [1] read("myInputFile.root", 10000 )

void read(TString inFileNames, int nEvents ) {
   
   // If the analysis solely uses TTree::Draw statements,
   // you don't need to load
   // the shared library. You will receive warnings such as
   // Warning in <TClass::TClass>: no dictionary for class Particle
   // is available
   // but these can be ignored. However if you wish to work with the event
   // objects themselves, the shared library must be loaded:
   // Load the shared library, if not done automaticlly:
   //   gSystem->Load("/path/to/libeicsmear.so" );
   
   // The TTrees are named EICTree.
   // Create a TChain for trees with this name.
   TChain tree("EICTree");
   
   // Add the file(s) we want to analyse to the chain.
   // We could add multiple files if we wanted.
   tree.Add(inFileNames); // Wild cards are allowed e.g. tree.Add("*.root" );
// tree.Add(/path/to/otherFileNames ); // etc... 
   
   // Create an object to store the current event from the tree.
   // This is how we access the values in the tree.
   // If you want to use generator-specific values, then
   // the event type should match the type in the TTree. Valid types are
   // EventPythia, EventPepsi, EventRapgap, EventDjangoh, EventMilou.
   // If you only need shared quantities like x, Q2 and the particle list
   // you can use EventBase and the macro will be general for any Monte Carlo.
   erhic::EventPythia* event(NULL);// = new EventPythia;
// EventBase* event(NULL);
   
   // Now associate the contents of the branch with the buffer.
   // The events are stored in a branch named event:
   tree.SetBranchAddress("event", &event ); // Note &event, not event.
   
   // Now we can do some analysis...
   
   // We record the largest particle pT we find here:
   double highestPt(-1. );
   
   // Histograms for our analysis.
   TH2D Q2VsX("Q2VsX",
              "Q^{2} vs. Bjorken x;log_{10}(x);log_{10}(Q^{2})",
               100, -5., 0., 100, -2., 3. );
   TH1D ptHist("ptHist",
               "pT of charged pions",
               50, 0.0, 2.0 );
   TH1D deltaPhi("deltaPhi",
                 "Delta-phi of hadrons",
                 40, 0.0, 3.1415 );
   
   // Loop over events:
   for(int i(0); i < nEvents; ++i ) {
      
      // Read the next entry from the tree.
      tree.GetEntry(i);
      
      // Fill the Q2 vs. x histogram:
      Q2VsX.Fill(TMath::Log10(event->GetX()),
                 TMath::Log10(event->GetQ2()));
      
      // The event contains a vector (array) of particles.
      int nParticles = event->GetNTracks();
      
      // We now know the number of particles in the event, so loop over
      // the particles:
      for(int j(0); j < nParticles; ++j ) {
         const Particle* particle = event->GetTrack(j);
         // Let's just select charged pions for this example:
         int pdg = particle->GetPdgCode();
         if(abs(pdg) != 211 ) continue;
         
         ptHist.Fill(particle->GetPt());
         
         // Update the highest pT:
         if(particle->GetPt() > highestPt ) {
            highestPt = particle->GetPt();
         } // if
      } // for
   } // for
   
   std::cout << "The highest pT was " << highestPt << " GeV/c" << std::endl;
   
   TCanvas canvas;
   Q2VsX.Draw("colz" );
   canvas.Print("Q2VsX.png" );
   ptHist.Draw();
   canvas.Print("pt.png" );
}
