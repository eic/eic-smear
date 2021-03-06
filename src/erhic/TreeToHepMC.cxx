/**
 \file
 Defines the main TreeToHepMC function.

 \author    Kolja Kauder
 \date      2021-07-07
 \copyright 2021 Brookhaven National Lab
 */

#include <string>
#include <iostream>

#include <TString.h>
#include <TSystem.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TDatabasePDG.h>

#include "eicsmear/erhic/Forester.h"
#include "eicsmear/erhic/File.h"

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenCrossSection.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/Attribute.h"

using std::cout;
using std::cerr;
using std::endl;

using HepMC3::FourVector;
using HepMC3::GenRunInfo;
using HepMC3::GenEvent;
using HepMC3::GenParticle;
using HepMC3::GenParticlePtr;
using HepMC3::GenVertex;
using HepMC3::GenVertexPtr;


/**
 This function converts our tree format to HepMC3
 It would be better to skip the ROOT step, but that
 would require a lot of duplication and/or refactorization
 */
Long64_t TreeToHepMC(const std::string& inputFileName,
		     const std::string& outputDirName,
		     Long64_t maxEvent) {

  // Get the input file name, stripping any leading directory path via
  // use of the BaseName() method from TSystem.
  TString outName = gSystem->BaseName(inputFileName.c_str());
  
  // Make sure this is a root file, 
  if ( !outName.EndsWith(".root", TString::kIgnoreCase) ){
    cerr << "Warning: " << inputFileName << " does not end with .root" << endl;
  }
  
  // Replace the extension
  outName.Replace(outName.Last('.'), outName.Length(), "");
  outName.Append(".hepmc");

  TString outDir(outputDirName);
  if (!outDir.EndsWith("/")) outDir.Append('/');
  outName.Prepend(outDir);

  // Open the input file and get the Monte Carlo tree from it.
  // Complain and quit if we don't find the file or the tree.
  TFile inFile(inputFileName.c_str(), "READ");
  if (!inFile.IsOpen()) {
    std::cerr << "Unable to open " << inputFileName << std::endl;
  }  // if
  TTree* mcTree(NULL);
  // TODO: Extend to smeared trees
  inFile.GetObject("EICTree", mcTree);
  if (!mcTree) {
    std::cerr << "Unable to find EICTree in " << inputFileName << std::endl;
    return 1;
  }
  erhic::EventMC* inEvent(NULL);
  mcTree->SetBranchAddress("event",&inEvent);

  // Get generator name
  TClass* branchClass = TClass::GetClass(mcTree->GetBranch("event")->GetClassName());
  TString generatorname = branchClass->GetName();
  if (branchClass->InheritsFrom("erhic::EventDis")) {
    generatorname.ReplaceAll("erhic::Event","");
  } else {
    cerr << branchClass->GetName() << " is not supported." << endl;
    return -1;
  }  // if

  if (branchClass->InheritsFrom("erhic::EventBeagle")) {
    cout << "BeAGLE input is currently not supported (can't fix mother-daughter structure yet)" << endl;
    return -1;    
  }

  // Run info
  std::shared_ptr<GenRunInfo> run = std::make_shared<GenRunInfo>();
  struct GenRunInfo::ToolInfo generator={
    std::string(generatorname),
    std::string("unknown version"),
    std::string("Used generator")
  };
  run->tools().push_back(generator);
  // Can be used to save the name of the control card if known (not usually the case)
  // struct HepMC3::GenRunInfo::ToolInfo config={cardname,"1.0",std::string("Control cards")};

  // can be used for customized weight names
  // currently only DEMP uses weights, named "weight"
  // We'll need to catch that later but for here just use the default
  std::vector<std::string> wnames;
  if (!wnames.size()) wnames.push_back("default");
  run->set_weight_names(wnames);
  
  // cross-section et al are stored as special strings
  // We don't have incremental information, so attach the full info to the header
  // instead of using HepMC3::GenCrossSection
  // crossSection is in mbarn!
  // The super set, not all generators supply all of these
  std::vector <string> RunAttributes = {"crossSection", "crossSectionError", "nEvents", "nTrials" };
  for ( auto att : RunAttributes ){
    TObjString* ObjString(nullptr);
    inFile.GetObject(att.c_str(), ObjString);
    if (ObjString) {
      double value = std::atof(ObjString->String());
      cout << " Adding to the header: " << att << "  " << value << endl;
      run->add_attribute( att, std::make_shared<HepMC3::DoubleAttribute>( value )) ;
    }
  }

  // Open the output file.
  HepMC3::WriterAscii file(outName.Data(),run);

  // Event Loop
  if (mcTree->GetEntries() < maxEvent || maxEvent < 1) {
    maxEvent = mcTree->GetEntries();
  } 
  std::cout <<
  "/-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-/"
  << std::endl;
  std::cout <<
  "/  Commencing conversion of " << maxEvent << " events."
  << std::endl;
  std::cout <<
  "/-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-/"
  << std::endl;
  for (Long64_t i(0); i < maxEvent; i++) {
    if (i % 10000 == 0 && i != 0) {
      std::cout << "Processing event " << i << std::endl;
    } 
    mcTree->GetEntry(i);
    
    // Construct new empty HepMC3 event and fill it.
    // Using GeV and cm (!)
    GenEvent hepmc3evt( HepMC3::Units::GEV, HepMC3::Units::CM );
    hepmc3evt.set_event_number(i);
    hepmc3evt.weights().clear();
    hepmc3evt.weights().push_back(1.0);

    // Go through event-wise variables
    // Leaves -> particles but also generator-specific variables
    auto leaves = mcTree->GetListOfLeaves();
    for ( int l = 0 ; l < leaves->GetEntries(); ++l ){
      TLeaf* leaf = (TLeaf*) leaves->At(l);
      TString lname = leaf->GetName();
      TString c = leaf->GetTypeName();
      if ( lname.BeginsWith("particles") ) continue;
      // cout << lname << "  " << leaf->GetValue() << endl;

      // Catch weight
      if ( lname == "weight"){
	hepmc3evt.weights().clear();
	hepmc3evt.weights().push_back( leaf->GetValue() );
	continue;
      }
      // Store generator variables - upconvert types
      if ( lname.Contains( "char", TString::kIgnoreCase) ) {
	// This can be a char type or a C string. I'm not aware
	// of either use case, so don't waste time to differentiate, just ignore
	continue;
      }
      if ( lname.Contains( "long", TString::kIgnoreCase) ) {
	hepmc3evt.add_attribute(lname.Data(),std::make_shared<HepMC3::LongAttribute>( leaf->GetValue() )) ;
	continue;
      }
      if ( lname.Contains( "int", TString::kIgnoreCase) 
	   || lname.Contains( "short", TString::kIgnoreCase) ) {
	hepmc3evt.add_attribute(lname.Data(),std::make_shared<HepMC3::IntAttribute>( leaf->GetValue() )) ;
	continue;
      }
      if ( lname.Contains( "float", TString::kIgnoreCase)
	   || lname.Contains( "double", TString::kIgnoreCase) ) {
	hepmc3evt.add_attribute(lname.Data(),std::make_shared<HepMC3::DoubleAttribute>( leaf->GetValue() )) ;
	continue;
      }
      // ignore everything else, e.g. bool      
    } // leaf list

    // First, fix sloppily implemented mother-daughter relations
    // Note: Multiple parents seem to only be in BeAGLE
    //       and somewtimes there seem to be exactly 2 parents, sometimes a range like for daughters.
    //       I cannot differentiate between the two, and for the exact case, it erroneously gives the impression
    //       of a large range, like 17 -- 254 which will wreak havoc on the graph.
    for( unsigned int t=0; t<inEvent->GetNTracks(); ++t) {
      const Particle* inParticle = inEvent->GetTrack(t);
      
      // Do my children know me?
      auto myindex = inParticle->GetIndex();
      // std::cout << "Processing track " << t << " with index " << myindex << std::endl;
      auto c1 = inParticle->GetChild1Index();
      auto cN = inParticle->GetChildNIndex();
      if ( cN==0 ) cN =c1;
      if ( c1>cN ) std::swap(c1,cN);
      if ( c1>0 ) {
	for ( UShort_t c = c1; c<=cN; ++c ){ // sigh. index starts at 1, tracks at 0;
	  Particle* child = inEvent->GetTrack(c-1);
	  // std::cout << "     Processing child with index " << child->GetIndex() << std::endl;
	  auto p1 = child->GetParentIndex();
	  auto pN = child->GetParentIndex1();
	  if ( p1>pN ) std::swap(p1,pN);
	  if ( p1==0 && pN==0 ){ // child erroneously believes to be motherless
	    child->SetParentIndex( myindex );
	  } else if ( p1==0 ) { // We are the only parent, is it correctly assigned?
	    if ( pN != myindex ){
	      // Nothing we can do, e.g. pythia allows multiple parenthood but lacks a way to describe that, see:
	      // 12     12       2101        5       18       31
	      // ...
	      // 16     11          2       10       18       31
	      // ...
	      // 26      1       -211       12        0        0
	      // 27      1        211       16        0        0
	      // cerr << "My child thinks its mother is " << pN << ", but it should be " << myindex << endl;
	      // return -1;
	    }
	  } else { // We have more than one parent, are they correct?
	    if ( myindex < p1 || myindex > pN ){
	      std::cout << "Processing event " << i << std::endl;
	      std::cout << "Processing track " << t << " with index " << myindex << std::endl;
	      std::cout << "     Processing child with index " << child->GetIndex() << std::endl;
	      cerr << "My child thinks its mothers range between " << p1 << " and " << pN
		   << ", but I am " << myindex << endl;
	      return -1;
	    }
	  }
	}
      }
      // Do my parents acknowledge me?
      auto p1 = inParticle->GetParentIndex();
      auto pN = inParticle->GetParentIndex1();
      if ( p1>pN ) std::swap(p1,pN);
      if ( p1==0 ) p1 =pN;
      // Do all my parents acknowledge me?
      if (pN > 0){
	for ( unsigned int p = p1; p<=pN ; ++p ){
	  Particle* parent = inEvent->GetTrack(p-1);
	  auto pc1 = parent->GetChild1Index();
	  auto pcN = parent->GetChildNIndex();
	  if ( pc1>pcN ) std::swap(pc1,pcN);
	  if ( pc1 > myindex ){
	    // cout << "hello1" << endl;
	    parent->SetChild1Index( myindex );
	  }
	  if ( pcN < myindex ){
	    // cout << "hello2" << endl;
	    parent->SetChildNIndex( myindex );
	  }
	}
      }      
    }

    // Perform consistency checks and collect particles
    std::vector<GenParticlePtr> hepevt_particles;
    hepevt_particles.reserve( inEvent->GetNTracks() );

    for( unsigned int t=0; t<inEvent->GetNTracks(); ++t) {
      const Particle* inParticle = inEvent->GetTrack(t);
      // Particles with status 1 cannot have children
      auto status = inParticle->GetStatus();
      if ( status==1 ){
	if (inParticle->GetNChildren() != 0 ){
	  cout << "Status is 1 but we have children?" << endl;
	  inParticle->Print();
	} 
      }

      // // All child-less particles should have a "safe" status (like 21), best would be 1
      // if (inParticle->GetNChildren() == 0 ){
      // 	if ( status !=1 && status !=21 ){
      // 	  cerr << "Status is " << status << " but we have no children?" << endl;
      // 	  inParticle->Print();
      // 	}
      // }

      // // Mother-less particles should be the beam only
      // Alas, that's not the case :-/
      // if ( t>1 && inParticle->GetParentIndex()==0 && inParticle->GetParentIndex1() ==0 ){
      // 	cout << "Event: " << i << " -- We have no mother" << endl;
      // 	inParticle->Print();
      // }

      FourVector pv = FourVector( inParticle->GetPx(), inParticle->GetPy(), inParticle->GetPz(),inParticle->GetE() );
      int statusHepMC = inParticle->GetStatus();
      // We should use only 1 (final), 2 (decayed hadron or lepton), 4 (beam), and >10, <=200 (anything else)
      // This may need to be decided on a generator-by-generator basis
      // We can assume final particles already have status 1, because that's
      // what EventMC::FinalState uses (and it's not overridden in existing classes)

      // Catch decayed leptons and hadrons
      if ( t>3 ){     // Ignore the beam
	if (inParticle->GetNChildren() != 0 ){ // ignore final particles
	  auto pdg = TDatabasePDG::Instance()->GetParticle( inParticle->Id() );
	  if ( pdg ){ // ignore unknown particles (e.g. pomerons, ions)
	    if ( TString(pdg->ParticleClass()).Contains("Lepton")
		 || TString(pdg->ParticleClass()).Contains("Baryon")
		 || TString(pdg->ParticleClass()).Contains("Meson")
		 ){
	      // Now our status should be 2!
	      // cout << statusHepMC << endl;
	      // inParticle->Print();
	      statusHepMC = 2;
	    }
	  }
	}
      }
      // force everything else to be legal
      if ( statusHepMC != 1 && statusHepMC != 2 && statusHepMC != 4 ){
	while ( statusHepMC <=10 ) statusHepMC+=10;
	while ( statusHepMC >200 ) statusHepMC-=10;
      }

      // Create GenParticle
      hepevt_particles.push_back( std::make_shared<GenParticle>( pv, inParticle->Id(), statusHepMC ));
      hepevt_particles.back()->set_generated_mass( inParticle->GetM() );

    }

    // Build the event
    // beam particles
    // --------------
    // Default is 1 = e-, 2 = hadron, 3 = scattered e-, 4 = exchange boson
    // But this can (and does) differ, especially for the scattered lepton
    // As always, be aware of Fortran starting to count at 1
    // Vertex: We don't keep track of time
    auto lepton=inEvent->BeamLepton();
    int index_lepton = lepton->GetIndex();
    if ( index_lepton !=1 ) std::cout << "Warning: Found BeamLepton at " << index_lepton << endl;
    auto hep_lepton = hepevt_particles.at( index_lepton-1);
    hep_lepton->set_status(4);

    auto boson=inEvent->ExchangeBoson();
    int index_boson = boson->GetIndex();
    // This happens in Sartre who puts the boson at 3
    // if ( index_boson !=4 ) std::cout << "Warning: Found ExchangeBoson at " << index_boson << endl;
    // if ( boson->GetParentIndex() != index_lepton && boson->GetParentIndex1() != index_lepton ){
    //   // This is common for Sartre, and any others that treat the boson like the beam
    //   std::cout << "Warning: ExchangeBoson doesn't recognize the beam as its mother " << endl;
    // }
    auto hep_boson = hepevt_particles.at( index_boson-1);
    // if needed / desired, could force hep_boson->set_status(4);

    GenVertexPtr v_lepton = std::make_shared<GenVertex>();
    v_lepton->add_particle_in  (hep_lepton);
    v_lepton->add_particle_out (hep_boson);
    hepmc3evt.add_vertex(v_lepton);
      
    auto hadron=inEvent->BeamHadron();
    int index_hadron = hadron->GetIndex();
    if ( index_hadron !=2 ) std::cout << "Warning: Found BeamHadron at " << index_hadron << endl;
    auto hep_hadron = hepevt_particles.at( index_hadron-1);
    GenVertexPtr v_hadron = std::make_shared<GenVertex>();
    hep_lepton->set_status(4);
    v_hadron->add_particle_in (hep_hadron);
    hepmc3evt.add_vertex(v_hadron);

    // Now work our way through the remaining particles
    // - Attach each particle that has a mother to their end vertex
    // ---> Create / overwrite production vertex in the process.
    //      If it's inconsistent, there's not much we can do
    // - attach motherless particles to the exchange boson
    // ---> In that case, leave the production vertex location in peace
    // Topological order should just translate to the fact that
    // children always have a higher index than their parents
    // Note: Multiple parents would wreak havoc here - need to think more about BeAGLE
    for( unsigned int t=0; t<inEvent->GetNTracks(); ++t) {
      const Particle* inParticle = inEvent->GetTrack(t);

      // Skip what we already have
      int index = inParticle->GetIndex();
      if ( index==index_lepton || index==index_boson || index==index_hadron) continue;
      auto hep_in = hepevt_particles.at( index-1);
      
      // Mother?
      auto hep_mom = hep_boson;
      int momindex = inParticle->GetParentIndex();
      if ( momindex > 1 ){
	hep_mom = hepevt_particles.at( momindex-1);
      }
      
      // Does mom have an end vertex yet?
      auto momend = hep_mom->end_vertex();
      if (!momend) {
	momend = std::make_shared<GenVertex>();
	momend->add_particle_in(hep_mom);
	hepmc3evt.add_vertex(momend);
      }

      momend->add_particle_out(hep_in);
	
      // update prod vertex?
      if ( momindex > 1){
	auto vnew = inParticle->GetVertex();
	momend->set_position( FourVector( vnew.x(), vnew.y(), vnew.z(), 0));
      }

    }
    // Done! Write the event.
    file.write_event(hepmc3evt);

    // There's a bunch of cleanup one should do now, with all the dynamical
    // vertices and particles. BUT shared_ptr should take care of that. Revisit if there are memory leaks.
			      
    
  }  // event loop

    
  
  Long64_t result = 0;
  
  return result;
}
