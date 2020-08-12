# EIC-smear

## About

Documentation in the BNL Wiki
* https://wiki.bnl.gov/eic/index.php/Monte_Carlo_and_Smearing

Contacts
* Alexander Kiselev <ayk@bnl.gov>
* Kolja Kauder <kkauder@bnl.gov>

#### Overview ####

eic-smear is a Monte Carlo analysis package originally developed by the BNL EIC task force.

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

#### Smearing scripts ####
Detector parameterizations have been moved to [their own repository](https://github.com/eic/eicsmeardetectors) and need to be installed separately.

#### Building ####

##### Prerequisites #####

* CMake version >= 3.1 is required.
* Compiler with C++11 support
* ROOT version >= 6.10 is required, testing was only done for at least 6.14.

##### Procedure #####

Create a directory in which to build eic-smear and navigate to that
```
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

Detector parameterizations are developed on a different timescale and
are more prone to changes by users; they have therefore been moved to
https://github.com/eic/eicsmeardetectors, together with tests and
examples You should head over there now and install them.

A wrapper allows to start ROOT with libraries loaded and displays
version information as well as the library locations.
```
$ eic-smear
Using eic-smear version: 1.1.0
Using these eic-smear libraries :
/Users/kkauder/software/lib/libeicsmear.dylib
/Users/kkauder/software/lib/libeicsmeardetectors.dylib
eic-smear [0]
```

It can also be used for simple one liners:
```
echo 'BuildTree ("ep_hiQ2.20x250.small.txt.gz");SmearTree(BuildMatrixDetector_0_1(),"ep_hiQ2.20x250.small.root")' | eic-smear
```

##### Notes: #####

* If you see instances of things like
```
Error in cling::AutoloadingVisitor::InsertIntoAutoloadingState:
   Missing FileEntry for eicsmear/smear/Smear.h
   requested to autoload type erhic::VirtualParticle
```
please ```setenv``` or ```export``` the environment variable ROOT_INCLUDE_PATH to point to the include directory in your installation. It should no longer be necessary with recent build improvements, but for ROOT versions above 6.20, the necessity returns.

* If building at BNL, you can get ROOT6 in the following manner
```
source /afs/rhic.bnl.gov/eic/restructured/etc/eic_cshrc.csh
setenv EIC_LEVEL pro
#verify
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

* ```BuildTree()``` supports HepMC2 and HepMC3 input, if a HepMC3 installation is found by cmake. One can also pass a specific installation directory:
```
cmake ../ -DCMAKE_INSTALL_PREFIX=</path/to/install> -DHepMC3=/path/to/HepMC3
```
The filename should contain `hepmc`, the reader determines the used
version automatically.

* For some reason, tab completion inside ROOT currently only works after
explicitly loading the library
```
root [] gSystem->Load("libeicsmear");
```
even if that same command is in your rootlogon.C.

* The above assumes you are using ```csh```. In ```bash```, simply replace ```setenv A B``` with ```export A=B```.


#### Doxygen ####

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
not currently in use. The original Subversion repository (not up to date) is at:
```
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
