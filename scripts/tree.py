# Created by tpb on 15th Feb 2012

"""Example routines for creating and reading a ROOT tree.

Python is simply a superior way to write many ROOT scripts than
using normal CINT macros. It is more stable, and the simpler
syntax of Python compared to C++ is better suited to common
tree and histogramming operations.
This module contains simple routines outlining how to create
and access tree files with EIC events via Python.
Note you need a version of ROOT compiled against Python and
the necessary environment variables set (see the ROOT installation
instructions).
"""

def load():
    """Import the PyROOT module and use ROOT.gSystem.Load
    to import the eic-smear class library.
    """
    import ROOT
    ROOT.gSystem.Load('.libs/libeicsmear') # Or whatever your path is


def build(inputname, outputdir = '.', nevents = -1):
    """Example of creating a tree.

    Once the eic-smear library, you can just use the
    BuildTree routine as usual.
    """
    import ROOT
    ROOT.gSystem.Load('.libs/libeicsmear') # Or whatever your path is
    
    ROOT.BuildTree(inputname, outputdir, nevents)


def write(outputname, nevents):
    """Example of creating a tree.

    Creating and writing events manually without the use of BuildTree.
    """
    import ROOT
    ROOT.gSystem.Load('.libs/libeicsmear') # Or whatever your path is

    file = ROOT.TFile(outputname, 'recreate')
    tree = ROOT.TTree('events', 'A ROOT tree of EIC events')

    # Create a branch buffer of whatever event type you fancy
    event = ROOT.erhic.EventPythia()
    tree.Branch('event', event)

    for i in xrange(0, nevents):
        event.SetN(i)
        # ...
        # Build your event
        # ...
        tree.Fill()

    file.Write()


def read(inputname, treename):
    """Example of reading back a tree.
    """
    import ROOT
    ROOT.gSystem.Load('.libs/libeicsmear') # Or whatever your path is

    file = ROOT.TFile(inputname, 'read')

    # Note that you can access the tree simply via file.treename
    file.Get(treename).Draw('GetN()')

