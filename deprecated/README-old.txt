eic-smear
https://wiki.bnl.gov/eic/index.php/Monte_Carlo_and_Smearing

Contacts:
elke@bnl.gov
tpb@bnl.gov
macl@bnl.gov

--------------------------------------------------------------------------------
Overview

eic-smear is the BNL EIC task force Monte Carlo analysis package.
It contains classes and routines for:
1) Building events in a C++ object and writing them to a ROOT file in a tree
   data structure.
2) Performing fast detector smearing on those Monte Carlo events.

The tree-building portion processes plain text files, formatted according to
the EIC convention, into a ROOT TTree containing events.
The following Monte Carlo generators are supported:
PYTHIA, RAPGAP, PEPSI, DJANGOH, MILOU, DPMJet, gmc_trans.
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

--------------------------------------------------------------------------------
Building

For those with access to the EIC AFS area, the stable versions of the
code are installed under
/afs/rhic.bnl.gov/eic/PACKAGES/eic-smear

Those who wish to build the code themselves can check it out via Subversion:

svn checkout http://svn.racf.bnl.gov/svn/eic/Utilities/eic-smear/trunk eic-smear

Navigate to /path/to/eic-smear, and run

autoreconf --force --install

Create a directory in which to build eic-smear, navigate to that and run

/path/to/eic-smear/configure
make

If you wish to install the generated libraries in a location install-dir,
instead do

/path/to/eic-smear/configure --prefix=/path/to/install-dir
make
make install

If you want to build PYTHIA6-dependent components, pass the location
of libPythia6 to configure via
/path/to/eic-smear/configure --with-pythia6-libdir=/path/to/pythia6/lib