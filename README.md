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
* DPMJet
* gmc_trans

Please see the EIC wiki pages for descriptions of the data formats and
the EIC AFS area for installations of each of the Monte Carlo generators.
Each entry in the TTree is a single C++ event object, storing event-wise 
quantities common to all generators, generator-specific variables and
a list of particle objects.

The smearing portion of the code provides a collection of classes and routines
that allow simple parameterisations of detector performance, provided via a
ROOT script, to be applied to any of the above Monte Carlo
events. Detector acceptance can also be defined. When run, a ROOT file with
smeared events is produced, which can be analysed in conjunction with the
input Monte Carlo events.

Both portions of the code are included in the eic-smear shared library.


## Building

## Prerequisites

* CMake version >3.1 is required.
* ROOT6 is required. This implies you need a compiler with c++11 support

If building at BNL, get ROOT6 in the following manner
```source /afs/rhic.bnl.gov/eic/restructured/etc/eic_cshrc.csh
setenv EIC_LEVEL dev
# verify
which root
```


### Procedure
Create adirectory in which to build eic-smear and navigate to that
```sh
cd eic-smear
mkdir build
cd build
```

Configure using cmake. Optionally, you can specify a location where to
install include files and libraties:
```
cmake ../ -DCMAKE_INSTALL_PREFIX=</path/to/install>
```

Build and install (the -j flag specifies how many parallel compilation
threads to use)
```
make -j 2
make install
```

### Pythia
If you want to build PYTHIA6-dependent components, pass the location
of libPythia6 to cmake:
```
cmake ../ -DCMAKE_INSTALL_PREFIX=</path/to/install> -DPYTHIA6_LIBDIR=/path/to/pythia6/lib
```


For those with access to the EIC AFS area, the stable versions of the
code are installed in the current pro environment.


## A historical note

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
