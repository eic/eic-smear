
void analysis()
{
  gSystem->Load(".libs/libeicsmear.so");

  // Input simulated & reconstructed files;
  TFile *ff = new TFile("head.10event.root");
  //TFile *ff = new TFile("pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.8..100k-lines.10event.root");
  TTree *tree = ff->Get("EICTree"); 

  // Branches of interest;
  erhic::EventBeagle *event = new erhic::EventBeagle();
  //erhic::EventPythia *event = new erhic::EventPythia();
  tree->SetBranchAddress("event", &event);

  // Loop through all events; NB: for box-generated events without secondaries 
  // could simply use cbmsim->Project() as well; in general EicEventAssembler 
  // should be used for "true" physics events;
  int nEvents = tree->GetEntries(); printf("%d event(s)\n", nEvents);
  for(unsigned ev=0; ev<nEvents; ev++) {
    tree->GetEntry(ev);

    printf(" -> %3d particle(s)\n", event->GetNTracks());

    for(unsigned pt=0; pt<event->GetNTracks(); pt++) {
      erhic::ParticleMC *particle = event->GetTrack(pt);
      printf("   pt#%3d (%6d)    -> %7.3f\n", pt, (int)particle->Id(), particle->GetPz());
      if (particle->eA) 
	//printf("(%3d .. %3d)!\n", (int)particle->eA->orig1, (int)particle->orig1);
	printf("(%3d .. %3d)!\n", particle->orig1, particle->orig);
	//printf("(%8.3f .. %8.3f)!\n", particle->eA->pz, particle->pz);
      else
	printf("(%3d)!\n", particle->orig);
    } //for pt
  } //for ev
} // analysis()
