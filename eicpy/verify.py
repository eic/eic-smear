#!/usr/bin/env python

from array import array
import os

import ROOT
try:
    ROOT.PyConfig.IgnoreCommandLineOptions = True
except AttributeError:
    pass
ROOT.gSystem.Load('libeicsmear')
ROOT.gROOT.SetBatch(True)

def parse():
    """Parse command line arguments and return in an argparse.Namespace."""
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('file', nargs='+', help='input ASCII file')
    parser.add_argument('-o', '--outdir', default='.',
                        help='output directory')
    parser.add_argument('-f', '--format', default='pdf', choices=['ps', 'pdf'],
                        help='output file format')
    parser.add_argument('-n', '--nevents', default=-1, type=int,
                        help='number of events per file')
    return parser.parse_args()

def build(filename, outdir='.', nevents=-1):
    """Build a ROOT tree from the named text file.
    
    Returns the name of the resultant ROOT file.
    
    """
    base, extension = os.path.splitext(filename)
    ROOT.BuildTree(filename, outdir, nevents)
    return '.'.join([base, 'root'])

def histogram(typename, *args, **kwargs):
    """Return a histogram of the requested type.
    
    args -- normal histogram constructor arguments
    kwargs -- additional attributes to attach to the histogram
    
    e.g. to create a TH1D with a 'cut' attribute:
    
    histogram(ROOT.TH1D, 'myhist', '', 10, 0, 10, cut='x>10')
    
    """
    h = typename(*args)
    attributes = {'select': h.GetName(), 'cut': '', 'log': ''} # Defaults
    attributes.update(kwargs)
    for name, value in attributes.iteritems():
        setattr(h, name, value)
    h.SetOption(attributes.get('opt', ''))
    return h

def th1d(*args, **kwargs):
    """Return a TH1D"""
    return histogram(ROOT.TH1D, *args, **kwargs)

def th2d(*args, **kwargs):
    """Return a TH2D"""
    return histogram(ROOT.TH2D, *args, **kwargs)

def fill_histogram_from_tree(self, tree):
    """Fill a histogram from a tree.
    
    Uses the histogram's 'select' and 'cut' attributes to determine
    what is plotted.
    
    """
    try:
        counts = tree.Project(self.GetName(), self.select, self.cut)
    except:
        print 'Error filling histogram', self.GetName(),\
              'from tree, histogram will not be filled'

# Add as a new member function to TH1
ROOT.TH1.FillTree = fill_histogram_from_tree

class Histograms(object):
    """A collection of basic histograms."""

    def name(self, string):
        """Return the string prepended by this object's file basename."""
        return '_'.join([self.base, string])

    def __init__(self, rootfile):
        """Initialise histograms and associate them with a ROOT file."""
        rootfile.cd()
        self.file = rootfile
        self.base, extension = os.path.splitext(file.GetName())
        self.out = '.'.join([self.base, 'pdf'])
        name = self.name
        self.window = ROOT.TCanvas(name('canvas'), '', 1, 1, 800, 800)
        self.page(open=True)
        # erhic::EventDis variables
        self.histograms = [
            th1d(name('beamlepton'), 'Beam lepton ID', 300, -3000., 3000.,
                 select='BeamLepton().Id().Code()'),
            th1d(name('beamhadron'), 'Beam hadron ID', 300, -3000., 3000.,
                 select='BeamHadron().Id().Code()'),
            th1d(name('scattered'), 'Scattered lepton ID', 300, -3000., 3000.,
                 select='ScatteredLepton().Id().Code()'),
            th1d(name('boson'), 'Exchange boson ID', 300, -3000., 3000.,
                 select='ExchangeBoson().Id().Code()'),
            th1d(name('nu'), 'nu', 25, 1., 5., select='log10(nu)'),
            th2d(name('Q2x'), '$Q^{2}$ vs. $x_{Bj}$', 25, -5., 0., 25, -1., 4.,
                 select='log10(QSquared):log10(x)', opt='colz'),
            th1d(name('xf'), '$x_{F}$', 55, -1.1, 1.1, select='xFeynman',
                 cut='KS==1 && id!=ScatteredLepton().Id().Code()'),
            th1d(name('z'), '$z$', 60, -0.1, 1.1, select='z',
                 cut='KS==1 && id!=ScatteredLepton().Id().Code()', log='y'),
            th1d(name('pt'), '$p_{T}$ vs. virtual photon', 60, 0., 3.,
                 select='ptVsGamma', log='y'),
            th1d(name('ptot'), 'Total final-state momentum', 100, 0., 300.,
                 select='FinalStateMomentum().P()'),
            th1d(name('ctot'), 'Total final-state charge', 100, -2., 2.,
                 select='FinalStateCharge()'),
            th2d(name('pthb'), 'Get4VectorInHadronBosonFrame().Pt():ptVsGamma',
                 50, 0., 10., 50, 0., 10.,
                 select='Get4VectorInHadronBosonFrame().Pt():ptVsGamma',
                 opt='colz'),
            th2d(name('mass'), 'Photon mass$^{2}$ vs. $Q^{2}$',
                 25, -1., 4., 25, -1., 4.,
                 select='log10(-(ExchangeBoson().Get4Vector().M2())):'
                        'log10(QSquared)',
                 opt='colz')
        ]
        # x, Q2, W2 and y, for each calculation method
        for method in ['', 'JB', 'DA']:
            def var(string):
                """Insert method into a string."""
                if '{' in string:
                    return string.format(method)
                else:
                    return ''.join([string, method])
            self.histograms += [
                th1d(name(var('x')), var('$x_{{Bj}}{}$'), 25, -5., 0.,
                     select=var('log10(x{})')),
                th1d(name(var('Q2')), var('$Q^{{2}}{}$'), 25, -1., 4.,
                     select=var('log10(QSquared)'), log='y'),
                th1d(name(var('y')), var('$y{}$'), 25, -0.1, 1.1,
                     select=var('y')),
                th1d(name(var('W2')), var('$W^{{2}}{}$'), 25, 1., 6.,
                     select=var('log10(WSquared)'))
            ]
        # erhic::EventMC variables:
        self.histograms += [
            th1d(name('number'), 'Event number', 100, 0.,
                file.EICTree.GetEntries(), select='number'),
            th1d(name('process'), 'Process', 200, 0., 200., select='process'),
            th1d(name('nTracks'), 'Number of tracks', 200, 0., 200.,
                select='nTracks'),
            th2d(name('tracks'), 'particles.GetEntries() vs. nTracks',
                 200, 0., 200., 200, 0., 200.,
                 select='@particles.GetEntries():nTracks', opt='colz'),
            th1d(name('Ein'), 'ELeptonInNucl', 25, 0., 12000.,
                select='ELeptonInNucl'),
            th1d(name('Eout'), 'ELeptonOutNucl', 25, 0., 12000.,
                select='ELeptonOutNucl'),
            th2d(name('hadronic'),
                 'HadronicFinalStateMomentum().M2():WSquared',
                 25, 1., 6., 25, 1., 6.,
                 select='log10(HadronicFinalStateMomentum().M2()):'
                        'log10(WSquared)', opt='colz')
        ]

    def page(self, open=False, close=False, log=''):
        """Print the current contents of the window."""
        if open:
            self.window.Print(self.out + '[')
        elif close:
            self.window.Print(self.out + ']')
        else:
            self.window.SetLogx('x' in log)
            self.window.SetLogy('y' in log)
            self.window.SetLogz('z' in log)
            self.window.Print(self.out)

    def make(self):
        """Fill and print all histograms."""
        self.file.cd()
        for histogram in self.histograms:
            histogram.FillTree(self.file.EICTree)
            self.window.cd()
            histogram.Draw()
            self.page(log=histogram.log)

    def close(self):
        """Cloes the ROOT file and image file."""
        self.file.Write()
        self.page(close=True)


class PythiaHistograms(Histograms):
    """PYHTIA-specific histograms"""

    def __init__(self, file):
        super(PythiaHistograms, self).__init__(file)
        name = self.name
        self.histograms += [
            th2d(name('trueX'), '$x_{Bj}$', 25, -5., 0., 25, -5., 0.,
                 select='log10(x):log10(trueX)', opt='colz'),
            th2d(name('trueQ2'), '$Q^{2}$', 25, -1., 4., 25, -1., 4.,
                 select='log10(QSquared):log10(trueQ2)', opt='colz'),
            th2d(name('trueY'), '$y$', 25, -0.1, 1.1, 25, -0.1, 1.1,
                 select='y:trueY', opt='colz'),
            th2d(name('trueW2'), '$W^{2}$', 25, 1., 6., 25, 1., 6.,
                 select='log10(WSquared):log10(trueW2)', opt='colz')
        ]


class DjangohHistograms(Histograms):
    """DJANGOH-specific histograms"""

    def __init__(self, file):
        super(DjangohHistograms, self).__init__(file)
        name = self.name
        self.histograms += [
            th2d(name('dtrueX'), '$x_{Bj}$', 25, -5., 0., 25, -5., 0.,
                 select='log10(x):log10(dtrueX)', opt='colz'),
            th2d(name('dtrueQ2'), '$Q^{2}$', 25, -1., 4., 25, -1., 4.,
                 select='log10(QSquared):log10(dtrueQ2)', opt='colz'),
            th2d(name('dtrueY'), '$y$', 25, -0.1, 1.1, 25, -0.1, 1.1,
                 select='y:dtrueY', opt='colz'),
            th2d(name('dtrueW2'), '$W^{2}$', 25, 1., 6., 25, 1., 6.,
                 select='log10(WSquared):log10(dtrueW2)', opt='colz')
        ]


# Define Histograms to make for each generator
HISTOGRAMS = {'pythia': PythiaHistograms,
              'lepto': Histograms,
              'milou': Histograms,
              'djangoh': DjangohHistograms,
              'pepsi': Histograms,
              'gmc_trans': Histograms}

def event_type(file):
    """Return the event type for a ROOT file.
    
    e.g. 'pythia', 'lepto'.
    
    """
    file.EICTree.GetEntry(0)
    return file.EICTree.event.ClassName().lstrip('erhic::Event').lower()

if __name__ == '__main__':
    options = parse()
    # Build ROOT files
    names = [build(name, options.outdir, options.nevents)
             for name in options.file]
    # Open ROOT files for updating with histograms
    root_files = [ROOT.TFile(name, 'update') for name in names]
    # Initialise histograms
    histograms = [HISTOGRAMS[event_type(file)](file) for file in root_files]
    # Fill, draw and write histograms
    for histogram in histograms:
        histogram.make()
        histogram.close()