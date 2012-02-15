//
// build.cxx
//
// Created by TB on 6/13/11.
// Copyright 2011 BNL. All rights reserved.
//
// Demonstration of how to run BuildTree
// You can specify a maximum number of events to read:
// <= 0 means read all events (default).

void build(const TString inputFile, TString outputDir, Long64_t maxEvents) {
   
   // Load the shared library, if not done automaticlly:
//   gSystem->Load("/path/to/libeicsmear.so" );
   
   
   // Call the BuildTree() function.
   // The output file name is the same as the input file name but with the
   // extension ".root" e.g. myfile.txt --> myfile.root
   BuildTree(inputFile, outputDir, maxEvents);
}
