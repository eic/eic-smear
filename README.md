# EIC-smear

## About

Additional documentation can be found on the
[EIC GitHub pages](https://eic.github.io/software/eicsmear.html).

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
* Additionally, HepMC2 and HepMC3 files are supported, allowing for example Pythia8 and eSTARlight output to be processed.

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
zero reflecting information gathered only by a tracker.
Particles further have methods of the form
```c++
bool psmeared = p->IsPSmeared();
```
to differentiate between "not measured" and "measured as 0".
In such a
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
setenv EICDIRECTORY </path/to/install>
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
Using eic-smear version: 1.1.9
Using these eic-smear libraries :
/Users/kkauder/software/lib/libeicsmear.dylib
/Users/kkauder/software/lib/libeicsmeardetectors.dylib
eic-smear [0]
```

It can also be used for simple one liners:
```
echo 'BuildTree ("ep_hiQ2.20x250.small.txt.gz");SmearTree(BuildMatrixDetector_0_1(),"ep_hiQ2.20x250.small.root")' | eic-smear
```

One some architectures and ROOT versions, ```TRint``` has an obscure
bug that will cause segmentation faults when using ```std::cout``` and
similar commands inside this interpreter. Use printf instead, or just
load the libraries directly in a generic root instance.


#### Creation of HepMC output ####
eic-smear's ROOT trees can be converted to HepMC output using the `TreeToHepMC`
macro, e.g.:
```
$ eic-smear
Using eic-smear version: 1.1.9
Using these eic-smear libraries :
/Users/kkauder/software/lib/libeicsmear.dylib
/Users/kkauder/software/lib/libeicsmeardetectors.dylib
eic-smear [0] BuildTree("pythia6.txt")
eic-smear [1] TreeToHepMC("pythia6.root")
```
or
```
$ echo 'BuildTree("pythia6.txt"); TreeToHepMC("pythia6.root")' | eic-smear

```
BeAGLE's structure is flattened. All intermediate particles are there,
but their parental struture may not be preserved.
For all generators, including BeAGLE, hadronic and leptonic decays
should be correct.

HepMC2 output is possible via a flag:
```
$ eic-smear
eic-smear [0] TreeToHepMC("beagle_eD.root",".",-1,true);
```


This macro has been tested to work with HepMC tools and rivet as well as eAST.
Any feedback for these or other generators is very welcome!

Note: The necessity to first generate ROOT trees as an intermediate step is
awkward but currently unavoidable.
BuildTree() is inextricably linked with many internal classes and functionalities.
Significant refactorization to eliminate it in the future is possible
but will take longer.


##### Notes: #####

* If you see instances of things like
```
Error in cling::AutoloadingVisitor::InsertIntoAutoloadingState:
   Missing FileEntry for eicsmear/smear/Smear.h
   requested to autoload type erhic::VirtualParticle
```
please ```setenv``` or ```export``` the environment variable ```ROOT_INCLUDE_PATH``` to point to the include directory in your installation. It should no longer be necessary with recent build improvements, but for ROOT versions above 6.20, the necessity returns.

* If building at BNL, you can get ROOT6 in the following manner
```
setenv EIC_LEVEL dev
source /afs/rhic.bnl.gov/eic/restructured/etc/eic_cshrc.csh
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

* A recently added script ```TreeToHepMC()``` can be used to transform
our ROOT trees to HepMC3 format.  TObjStrings for cross section etc. are saved in the RunInfo, all event generator-specific variables are saved in the event info. Parent-child relationships are repaired/reserved as much as possible, motherless particles get attached to the exchange boson.


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

## Developer Notes

### Code Conventions

There are clear style guidelines adhered to in the original code, for ease of maintenance and collaboration. Please adhere to all standards unless there are compelling reasons. Note that due to changing maintainers, rapid reactions to immediate issue requests, and things like replacement of deprecated and now removed features pre C++11, these guidelines aren't followed as strictly anymore. Nevertheless, please do your best.

#### Naming and file structure
Please observe the following:
* Header files should have a ```.h``` suffix and implementation files should have a ```.cxx``` suffix (following the ROOT convention).
* Stick to "one class, one file", or at most a few closely related classes in a single file. A file name should correspond to the class name, with declarations and definitions in separate header and implementation files. e.g. class ```MyAmazingClass``` should be declared in ```MyAmazingClass.h``` and be implemented in ```MyAmazingClass.cxx```.
* Keep all code in an appropriate namespace e.g. ```namespace erhic``` for general EIC and Monte Carlo code, ```namespace Smear``` for smearing-specific code.
* Files should be placed in the directory structure to match the namespace in which the code appears e.g. class ```erhic::EventDisplay``` would have its header in ```include/eicsmear/erhic/EventDisplay.h``` and its implementation in ```src/erhic/EventDisplay.cxx```.
* All names should follow **camelCase**:
  * Class names (and therefore filenames) and member function names are ```CapitalizedLikeThis```.
  * Non-class function names are ```capitalizedLikeThis```.
  * Member variable names are prefixed with a lower-case "m", and capitalized like ```mSomeMember```.

### Documentation
Make liberal use of [Doxygen](https://www.doxygen.nl/manual/docblocks.html) comments to document the code.
HTML documentation is generated automatically from these and is part of the central [EIC Doxygen](https://eic.github.io/doxygen/) page.
At a minimum give a brief description of each file and class, and preferably document all class methods (at least public ones).
Ask yourself whether a new user would be able to understand the basic purpose of each class or function you write, and how to use it, by looking at the provided comments; if not, write more!
Comments in the code at complicated or important points are encouraged to aid fellow developers.

#### Coding style
Beyond naming and file structure, coding conventions closely follow those of [Google](https://google.github.io/styleguide/cppguide.html), with a few exceptions:
* Streams are permitted, and encouraged in preference to functions such as ```printf()``` and ```scanf()```.
* Non-const reference function arguments are permitted.
There is a script, [cpplint](https://github.com/cpplint/cpplint), produced by Google to check files for compliance with style guidelines.
A (rather old) copy is included in the eicpy directory of the eic-smear distribution.
Run this to check for code compliance before committing changes.
It is not perfect, but catches most issues.
cpplint accepts a ```--filter``` argument to suppress warnings of different types.
To eliminate the exceptions to the Google style guide listed above, run cpplint as follows
```
 cpplint.py --filter=-runtime/references,-readability/streams MyFile.cxx
```

Some false positives that are known to occur include:
* Complaining about ROOT- or Doxygen-style inline comments e.g. ```//!``` (ROOT) and ```///<``` (Doxygen).
* ROOT headers will be misinterpreted as C system headers and cpplint will suggest to place them before C++ headers, but do **not** do so. The order of ```#include``` statements should be as follows:
  1. (Only in an implementation file) the name of the corresponding header file.
  1. C system headers.
  1. C++ system headers.
  1. 3rd-party package headers (but please think carefully before adding extra external dependencies to eic-smear).
  1. ROOT headers.
  1. eic-smear headers.

While cpplint can also be instructed to filter these types of warnings, this is not recommended as one may then miss genuine errors of that type.

#### Contributing
We follow the standard (GitHub Standard Fork & Pull Request Workflow )[https://gist.github.com/Chaser324/ce0505fbed06b947d962]. Alternatively to Fork'ing, you can also request to be added to the group of contributors and create a branch for a slightly more convenient work flow.

#### Versioning
Eic-smear is versioned according to (Semantic Versioning)[https://semver.org/], with some relaxation.
PATCH increases should never break backward compatibility. MINOR increases may make a few changes in existing smearing scripts necessary. MAJOR updates introduce significant new functionality and may seriously break backward-compatibility.




#### Historical, **deprecated** version control and installation instructions

Version control policy follows "best practices" from the SVN guide. To reiterate, as an example assume we have already released version 1.2, and are working on a new release 1.3:
1. New work is committed to /trunk: new features, bug fixes etc.
1. When the new release is (nearly) ready, copy /trunk to a release branch /branches/1.3 for final development version 1.3.
1. Testing of /branches/1.3 continues in parallel with new additions to /trunk. Port fixes for bugs between the two as they are found.
1. When testing of version 1.3 is complete, fix the new release by copying /branches/1.3 to /tags/1.3.0. *Do not modify * anything in /tags - these a fixed "snapshots" of the code.
1. Maintain /branches/1.3 with bug fixes ported from /trunk. When enough changes are made to warrant a new release, copy /branches/1.3 to /tags/1.3.1 etc.

To summarize: /trunk contains all newly added features and fixes. /branches/X.Y is the "maintenance" branch for version X.Y. /tags/X.Y.Z are fixed "snapshot" releases. /trunk and /branches/X.Y are modified, while /tags/X.Y.Z remain unchanged.

Feel free to create your own personal branches whenever you want to work on new features and fixes without interfering with /trunk. To make life easier, remember to frequently port changes from /trunk to your branch, to avoid problems when merging the branch back to /trunk. Once you are finished and have ported your new features to /trunk, the personal branch can be deleted.

After porting changes between /trunk and a branch with svn merge, always provide the following information in the message when you commit the change: the file or files modified; the revision from which the change came; a brief summary of the change; the source of the change. e.g.
```
  svn commit -m "AUTHORS: ported r3 (added list of names) from branches/1.3"
  svn commit -m "include/eicsmear/erhic/ParticleMC.h: ported r64 (fixed bug in calculation of Feynman x) from trunk"
```


These installation and version control instructions were used in the context of an older configuration
not currently in use. The original Subversion repository (not up to date) is at:
```
 http://svn.racf.bnl.gov/svn/eic/Utilities/eic-smear/trunk eic-smear
```

Old source tarballs are still available at
* [Version 1.0.3](https://wiki.bnl.gov/eic/upload/Eic-smear-1.0.3.tar.bz2)
* [Version 1.0.2](https://wiki.bnl.gov/eic/upload/Eic-smear-1.0.2.tar.bz2)
* [Version 1.0.1](https://wiki.bnl.gov/eic/upload/Eic-smear-1.0.1.tar.bz2)
* [Version 1.0.0](https://wiki.bnl.gov/eic/upload/Eic-smear-1.0.0.tar.bz2)


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
