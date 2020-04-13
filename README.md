# EIC-smear

## About

Documentation in the BNL Wiki
* https://wiki.bnl.gov/eic/index.php/Monte_Carlo_and_Smearing

Contacts
* Alexander Kiselev <ayk@bnl.gov>
* Kolja Kauder <kkauderl@bnl.gov>
* Maxim Potekhin <potekhin@bnl.gov>

## Overview

eic-smear is a Monte Carlo analysis package developed and used by the BNL EIC task force.

It contains classes and routines for:
1) Building events in a C++ object and writing them to a ROOT file in a tree data structure.
2) Performing fast detector smearing on those Monte Carlo events.

The tree-building portion processes plain text files, formatted according to
the EIC convention, into a ROOT TTree containing events.
The following Monte Carlo generators are supported:
* PYTHIA
* RAPGAP
* PEPSI
* DJANGOH
* MILOU
* BeAGLE
* Sartre
* DPMJet
* gmc_trans
Most of these are currently hosted at https://gitlab.com/eic/mceg.
Please see the associated documentation for further information on
individual generators.
Creation will typically be of the form
```
pythiaeRHIC < STEER_FILE > out.log
```
A few small example files are included for testing.


Each entry in the TTree named EicTree is a single C++ event object,
storing event-wise quantities common to all generators,
generator-specific variables and a list of particle objects.

The smearing portion of the code provides a collection of classes and routines
that allow simple parameterisations of detector performance, provided via a
ROOT script, to be applied to any of the above Monte Carlo
events. Detector acceptance can also be defined. When run, a ROOT file with
smeared events is produced.
To access truth level information, use the "friend" mechanism in ROOT,
see below.

The output tree, called Smeared, will mirror the behavior of a true
detector system, i.e. it will only contain entries for particles that
were smeared (=measured), and only partial information if only parts
were smeared. E. g., if only momentum is smeared, the energy field will be
zero reflecting information gathered only by a tracker. In such a
case, the analyzer can of course rely on the truth level, but a more
realistic approach would be to make same kind of assumptions that one
would have to make for a physical detector, such as assuming pion mass
for all tracks baring PID information etc.


Both portions of the code are included in the eic-smear shared
library.


## Building

### Prerequisites

* CMake version >3.1 is required.
* Compiler with C++11 support
* ROOT is required. ROOT6 strongly preferred, although ROOT5 may work

### Procedure
Create a directory in which to build eic-smear and navigate to that
```sh
cd eic-smear
mkdir build
cd build
```

Configure using cmake. Optionally, you can specify a location where to
install include files and libraries:
We'll assume that the installation path is in $EICDIRECTORY
```
setenv EICDIRECTORY=</path/to/install>
cmake ../ -DCMAKE_INSTALL_PREFIX=$EICDIRECTORY
```

Build and install (the -j flag specifies how many parallel compilation
threads to use)
```
make -j 2
make install
```

### Notes:

* If you see instances of things like
```
Error in cling::AutoloadingVisitor::InsertIntoAutoloadingState:
   Missing FileEntry for eicsmear/smear/Smear.h
   requested to autoload type erhic::VirtualParticle
```
please setenv or export the environment variable ROOT_INCLUDE_PATH to point to the include directory in your installation.

If building at BNL, get ROOT6 in the following manner
```source /afs/rhic.bnl.gov/eic/restructured/etc/eic_cshrc.csh
setenv EIC_LEVEL pro
# verify
which root
```

* If you want to build PYTHIA6-dependent components, pass the location
of libPythia6 to cmake:
```
cmake ../ -DCMAKE_INSTALL_PREFIX=</path/to/install> -DPYTHIA6_LIBDIR=/path/to/pythia6/lib
```
This will generate additional classes that allow creation of the
EicTree from within the framework, bypassing the text file generation.
More detailed documentation of this feature to follow.

* For some reason, tab completion inside ROOT currently only works after
explicitly loading the library
```
root [] gSystem->Load("libeicsmear");
```
even if that same command is in your rootlogon.C. 

## Tests and Examples
if you prepare building using the cmake option
```
-DBUILD_TESTS=ON
```
a variety of tests are generated from the tests directory:

```
./tests/test_simple/test_simple_buildtree
```
will read a (provided) e+D BeAGLE file.

```
./tests/particlegun
```
is a customizable particle gun that creates a few simplee histograms
and plots to see and test the acceptance dependence of smearing.

### A canonic example
A particularly useful example is 
```
./tests/qaplots
```
This starts from (provided) examples of text MC output, builds the
EicTree, smears it, and runs a simple analysis on the result. The
current set of QA plots will grow, but the source file
```
tests/qaplot.cxx
```
and associated cmake file
```
tests/CMakeLists.txt
```
is intended to be readable and give insight how to generate compiled
code.

### Deconstructed example
All steps can be performed interactively and individually in ROOT as
follows:

#### Generate EicTree
```
root [] gSystem->Load("libeicsmear");
root [] BuildTree ("tests/ep_noradcorr.20x250.small.txt",".",-1);
tests/ep_noradcorr.20x250.small.txt
```
BuildTree accepts the name of an input file, the output directory, and
the number of events to generate (-1 for all);
IMPORTANT: The file type is by default assumed to be pythia6. For
other files, please make sure to include the generator name in the
filename. Currently accepted are pythia, pepsi, lepto, rapgap, djangoh, beagle,milou, sartre, simple.


#### Smear the tree
```
root [] gSystem->Load("libeicsmear")
root [] .L smearBeAST.cxx // Assuming you copied this here
root [] SmearTree(BuildBeAST(), "ep_noradcorr.20x250.small.root", "smeared.root",-1)
/-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-/
/  Commencing Smearing of 10000 events.
/-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-/
```
BuildTree accepts the name of an input file, the output file, and
the number of events to smear (-1 for all);

* IMPORTANT: The file type is by default assumed to be pythia6. For
other files, please make sure to include the generator name in the filename.

#### Analyze the Tree
(Suppressing the root prompts for easier copy/paste):
```
root -l 
gSystem->Load("libeicsmear");
TFile mcf ("ep_noradcorr.20x250.small.root"); // truth
TTree* mc=(TTree*)mcf.Get("EICTree");
mc->AddFriend("Smeared","smeared.root"); // befriend
	
// Setup Input Event Buffer
erhic::EventMC* inEvent(NULL);
Smear::Event* inEventS(NULL);
mc->SetBranchAddress("event",&inEvent);
mc->SetBranchAddress("eventS",&inEventS);
	
// histo and event loop
TH2D* EEprime = new TH2D("EEprime", ";E;Eprime", 100, 0, 20, 100, 0, 20);
for(long iEvent=0; iEvent<mc->GetEntries(); iEvent++){    
    //Read next Event
    if(mc->GetEntry(iEvent) <=0) break;
    // Loop over Particles
    for(int j=0; j<inEventS->GetNTracks(); j++){
      if ( j<3 ) continue;       // Skip beam
      const Smear::ParticleMCS* inParticleS = inEventS->GetTrack(j); // Smeared Particle      
      const Particle* inParticle = inEvent->GetTrack(j); // Unsmeared Particle
      if(inParticleS == NULL) continue; 
      EEprime->Fill(inParticle->GetE(), inParticleS->GetE());
   }
}
EEprime->Draw("colz");
```

## Anatomy of a Smearer
A "detector" is constructed as follows. For details,
please also see examples included in the scripts/directory.
```c++
// ... omitted some includes and helpers
Smear::Detector BuildMyDetector() {
	
// Create the detector object to hold all devices
  Smear::Detector det;
  // The detector will calculate event kinematics from smeared values
  det.SetEventKinematicsCalculator("NM JB DA");
```
* Implement a function returning a Smear::Detector. Setup the detector
object and activate some additional calculation options.

```c++
  // Tracking
  // eta = -3.5 --  -2.5
  // sigma_p/p ~ 0.1% p+2.0%
  Smear::Acceptance::Zone TrackBack1Zone(ThetaFromEta ( -2.5 ),ThetaFromEta ( -3.5 ));
  Smear::Device TrackBack1P(Smear::kP, "sqrt( pow ( 0.001*P*P, 2) + pow ( 0.02*P, 2) )");
  TrackBack1P.Accept.AddZone(TrackBack1Zone);
  TrackBack1P.Accept.SetCharge(Smear::kCharged);
  det.AddDevice(TrackBack1P);
```
* Set up a device that smears. In this case, momentum is smeared
  * in eta = -3.5 --  -2.5,
  * with $`\sigma_p/p = 0.1 \% p+2.0 \%`$,
  * accepting all charged particles.

* IMPORTANT: For more realistic representation of electrons, it may
  make sense to only accept hadrons and create a separate device for
  electrons that represents the combined information from multiple
  detectors, e.g.:

```c++
  TrackBack1P.Accept.SetCharge(Smear::kCharged);
  TrackBack1P.Accept.SetGenre(Smear::kHadronic);
  det.AddDevice(TrackBack1P);
  
  Smear::Device EMTrackBack1P(Smear::kP, "sqrt(  0.01 * pow(p.Pt(), 2) + 0.005 * p.Pt())");
  EMTrackBack1P.Accept.AddZone(TrackBack1Zone);
  EMTrackBack1P.Accept.SetCharge(Smear::kCharged);
  TrackBack1P.Accept.SetGenre(Smear::kElectromagnetic);
  det.AddDevice(EMTrackBack1P);
```

* Continue on, adding $`\phi`$ and $`\theta`$ devices, calorimetry,
etc.

```
  return det;
}
```
* Finally, after adding all desired devices return the complete
detector.

* Formulas are based on ROOT::TFormula and accept kP, kPhi, kTheta, kE. In
principle, kPt and kPz is also supported but currently not working.
* IMPORTANT: Due to the limitations of TFormula only four different
  variables can be used at a time (because internally they get
  translated dynamically into the four available free variables
  x,y,z,t).
  



## Doxygen
A recent version of the detailed class documentation is temporarily
hosted at www4.rcf.bnl.gov/~eickolja/.
You can also create an up-to-date one by running
```
doxygen
```
in the top level directory and directing your web browser to the
created local file
```
doxygen/html/index.html
```
You can obtain doxygen at www.doxygen.nl.





## Historical note

These instructions were used in the context of an older configuration
not currently in use. The original Subversion reporitory (not up to date) is at:

```sh
 http://svn.racf.bnl.gov/svn/eic/Utilities/eic-smear/trunk eic-smear
```

It was configured using autoconf:
```sh
autoreconf --force --install
cd /path/to/eic-smear/
./configure  --prefix=/path/to/install-dir
make
make install
```

If you want to build PYTHIA6-dependent components, pass the location
of libPythia6 to configure via
```sh
/path/to/eic-smear/configure --with-pythia6-libdir=/path/to/pythia6/lib
```
