//
// draw.cxx
//
// Created by TB on 6/13/11.
// Copyright 2011 BNL. All rights reserved.
//
// Example of a simple analysis using TTree::Draw() statements.

void draw(const TString inputFile ) {
   
   // When you use only TTree::Draw() you can ignore errors like these:
   //Warning in <TClass::TClass>: no dictionary for class EventBase is available
   //Warning in <TClass::TClass>: no dictionary for class Particle is available
   // ROOT doesn't need the listed dictionary files to access the data.
   TFile file(inputFile, "READ");
   
   // The TTree is named EICTree
   TTree* tree(NULL );
   file.GetObject("EICTree", tree );
   if(! tree) return; // Oops!
   
   TCanvas* canvas = new TCanvas;
   canvas->Divide(2, 1);
   
   canvas->cd(1);
   // For event-wise quantities you don't need to prepend "event." to the
   // Draw() statement string; ROOT will automatically resolve it:
   tree->Draw("QSquared"); // Equivalent to "event.QSquared"
   
   canvas->cd(2);
   // Similarly you don't need to prepend "event.particles." to access particle
   // variables.
   tree->Draw("pt"); // Equivalent to "event.particles.pt".
}
