#!/usr/bin/env python

'''
Suite of standard QA plots for EIC trees to check they are filled properly.
It's fairly slow, but you don't need to run many events to check that all
values are being filled properly.
'''

import argparse
from collections import namedtuple

# Configure the environment for ROOT/eic-smear.
# Do this before importing ROOT, as initialise() will
# attempt to set the environment to find libPyROOT
# if it isn't already.
import eicpy
eicpy.initialise(False)

import ROOT

from binning import bin_log10

# Don't pass command line arguments to ROOT so it
# doesn't intercept use of --help: we want the help
# to be printed for this script, not ROOT.
# root.cern.ch/phpBB3/viewtopic.php?f=14&t=13582
ROOT.PyConfig.IgnoreCommandLineOptions = True

# Speed things up a bit
ROOT.SetSignalPolicy(ROOT.kSignalFast)

# No need to draw graphics to screen
ROOT.gROOT.SetBatch(True)

class LogAxis(namedtuple('LogAxis', 'x y z')):
    '''
    Utility for setting ROOT pad logarithmic axes.
    '''
    def apply(self, pad):
        '''
        Apply the log axis settings to a ROOT pad.
        '''
        pad.SetLogx(self.x != 0)
        pad.SetLogy(self.y != 0)
        pad.SetLogz(self.z != 0)

    def rebin(self, hist):
        if self.x:
            bin_log10(hist.GetXaxis())
        if self.y and hist.GetDimension() > 1:
            bin_log10(hist.GetYaxis())
        if self.z and hist.GetDimension() > 2:
            bin_log10(hist.GetZaxis())


class THWrapper(object):
    '''
    Simple wrapper around ROOT histogram adding some functionality.
    '''
    def __init__(self, h, log, option):
        '''
        Initialise with a histogram and a LogAxis giving whether
        each axis (x, y, z) should be log10 (1) or not (0).
        '''
        self.h = h
        self.log = log
        if not h.GetSumw2():
            self.h.Sumw2()
        self.h.SetOption(option)
        self.log.rebin(self.h)

    def draw(self, pad):
        '''
        Draw the histogram to a ROOT pad.
        '''
        initial = ROOT.gPad
        pad.cd()
        self.log.apply(pad)
        self.h.Draw()
        ROOT.gPad = initial


class EventHists(object):
    '''
    A collection of standard histograms for event-wise properties.
    '''
    def __init__(self):
        self.hists = {}
        self.hists['process'] = THWrapper(
            ROOT.TH1D('process', 'Event process number', 200, 0., 200.),
            LogAxis(0, 1, 0), '')
        self.hists['tracks'] = THWrapper(
            ROOT.TH1D('tracks', 'Number of tracks per event', 100, 0., 200),
            LogAxis(0, 0, 0), '')
#        self.hists['eInNucl'] = THWrapper(
#            ROOT.TH1D('eInNucl', 'E of incident lepton in nuclear rest frame',
#                      40, 0., 10000),
#            LogAxis(0, 0, 0), '')
#        self.hists['eScNucl'] = THWrapper(
#            ROOT.TH1D('eScNucl', 'E of scattered lepton in nuclear rest frame',
#                      40, 0., 10000),
#            LogAxis(0, 0, 0), '')
        self.hists['Q2x'] = THWrapper(
            ROOT.TH2D('Q2x', 'Q^{2} vs. x', 30, 1.e-5, 10., 16, 0.1, 1000.),
            LogAxis(1, 1, 1), 'colz')
        self.hists['y'] = THWrapper(
            ROOT.TH1D('y', 'y', 48, -0.1, 1.1), LogAxis(0, 1, 0), '')
        self.hists['W2'] = THWrapper(
            ROOT.TH1D('W2', 'W^{2} (GeV^{2})', 60, 0.1, 1.e5),
            LogAxis(1, 0, 0), '')
        self.hists['nu'] = THWrapper(
            ROOT.TH1D('nu', '#nu (GeV)', 60, 0.1, 1.e5), LogAxis(1, 0, 0), '')
        # Jacquet-Blondel and double-angle kinematics.
        # Clone the original histograms to maintain the same axis ranges.
        for i in ['y', 'Q2x', 'W2']:
            orig = self.hists[i]
            name = i + 'jb'
            self.hists[name] = THWrapper(orig.h.Clone(name), orig.log,
                                         orig.h.GetOption())
            self.hists[name].h.SetTitle(orig.h.GetTitle() +
                                        ' (Jacquet-Blondel)')
            name = i + 'da'
            self.hists[name] = THWrapper(orig.h.Clone(name), orig.log,
                                         orig.h.GetOption())
            self.hists[name].h.SetTitle(orig.h.GetTitle() +
                                        ' (double-angle)')

    def output(self, pad, filename):
        '''
        Print all histograms to an output PostScript file.
        '''
        initial = ROOT.gPad
        pad.cd()
        for j in sorted(self.hists):
            self.hists[j].draw(pad)
            pad.Print(filename)
        ROOT.gPad = initial

    def fill(self, event):
        '''
        '''
        self.hists['process'].h.Fill(event.GetProcess())
        self.hists['tracks'].h.Fill(event.GetNTracks())
#        self.hists['eInNucl'].h.Fill(event.ELeptonInNucl)
#        self.hists['eScNucl'].h.Fill(event.EScatteredInNucl)
        self.hists['Q2x'].h.Fill(event.GetX(), event.GetQ2())
        self.hists['y'].h.Fill(event.GetY())
        self.hists['W2'].h.Fill(event.GetW2())
        self.hists['nu'].h.Fill(event.GetNu())
        self.hists['Q2xjb'].h.Fill(event.GetXJacquetBlondel(),
                                   event.GetQ2JacquetBlondel())
        self.hists['yjb'].h.Fill(event.GetYJacquetBlondel())
        self.hists['W2jb'].h.Fill(event.GetW2JacquetBlondel())
        self.hists['Q2xda'].h.Fill(event.GetXDoubleAngle(),
                                   event.GetQ2DoubleAngle())
        self.hists['yda'].h.Fill(event.GetYDoubleAngle())
        self.hists['W2da'].h.Fill(event.GetW2DoubleAngle())


class TrackHists(object):
    '''
    A collection of standard histograms for track-wise properties.
    Optionally specify a list of PDG codes to accept.
    Only tracks with one of those codes will be used when filling histograms.
    '''
    def __init__(self):
        # Particle status
        self.KS = THWrapper(
            ROOT.TH1D('KS', 'Track status code', 110, -10., 100.),
            LogAxis(0, 0, 0), '')
        # Particle parent index
        self.orig = THWrapper(
            ROOT.TH1D('orig', 'Track parent index', 210, -10., 200.),
            LogAxis(0, 0, 0), '')
        # Child indicies
        self.child = THWrapper(
            ROOT.TH2D('child', 'Child track indices;first;last',
                      210, -10., 200., 210, -10., 200.),
            LogAxis(0, 0, 0), 'colz')
        # y vs. x momentum
        self.pypx = THWrapper(
            ROOT.TH2D('pypx', 'p_{y} vs. p_{x} (GeV/c);p_{x};p_{y}',
                      60, -3., 3., 60, -3., 3.),
            LogAxis(0, 0, 1), 'colz')
        # Transverse vs. longitudinal momentum
        self.ptpz = THWrapper(
            ROOT.TH2D('ptpz', 'p_{T} vs. p_{z} (GeV/c);p_{z};p_{T}',
                      60, -40., 260., 60, -1., 11.),
            LogAxis(0, 0, 1), 'colz')
        # Energy vs. mass
        self.em = THWrapper(
            ROOT.TH2D('em', 'Energy vs. mass (final-state only);m;E',
                      60, -1., 11., 60, -40., 260.),
            LogAxis(0, 0, 1), 'colz')
        # Vertex y vs. x
        self.vxy = THWrapper(
            ROOT.TH2D('vxy', 'Track vertex y vs. x',
                      51, -100., 100., 51, -100., 100.),
            LogAxis(0, 0, 1), 'colz')
        # Energy fraction (z) vs. total momentum
        self.zp = THWrapper(
            ROOT.TH2D('zp', 'Energy fraction vs. p (final-state only);p;z',
                      60, -40., 260., 48, -0.1, 1.1),
            LogAxis(0, 0, 1), 'colz')
        # Angles theta vs. phi
        self.thetaphi = THWrapper(
            ROOT.TH2D('thetaphi', '#theta vs. #phi (rad)',
                      35, -0.4, 6.6, 40, -0.5, 3.5),
            LogAxis(0, 0, 1), 'colz')
        # Pseudorapidity vs. rapidity
        self.etay = THWrapper(
            ROOT.TH2D('etay', '#eta vs. rapidity',
                      40, -20., 20., 40, -20., 20.),
            LogAxis(0, 0, 1), 'colz')
        # p_T vs. polar angle in photon-proton frame
        self.gamma = THWrapper(
            ROOT.TH2D('gamma', 'p_{T} vs. #theta in #gamma-p frame',
                      40, -0.5, 3.5, 35, -1., 6.),
            LogAxis(0, 0, 1), 'colz')
#        self.phi'] = THWrapper(
#            ROOT.TH1D('phi', '#phi in #gamma-p frame', 35, -0.4, 6.6),
#            LogAxis(0, 0, 0), '')
        # Feynman x
        self.xf = THWrapper(
            ROOT.TH1D('xf', 'Feynman x (final-state only)', 110, -1.1, 1.1),
            LogAxis(0, 0, 0), '')

    def output(self, pad, filename):
        '''
        Print all histograms to an output PostScript file.
        '''
        initial = ROOT.gPad
        pad.cd()
        hists = [self.KS, self.child, self.em, self.etay, self.gamma, self.orig,
            self.ptpz, self.pypx, self.thetaphi, self.vxy, self.xf, self.zp]
        for j in hists:
            j.draw(pad)
            pad.Print(filename)
        ROOT.gPad = initial

    def fill(self, tracks):
        '''
        Fill all histograms for a list of tracks
        '''
        # OK, this is a bit ugly looking, but localising all the track
        # and histogram methods outside the for loop saves quite a lot
        # of time when looping over millions of particles.
        # Here are all the Particle accessor methods we use:
        Particle = ROOT.erhic.ParticleMC
        status = Particle.GetStatus
        parent = Particle.GetParentIndex
        child1 = Particle.GetChild1Index
        childN = Particle.GetChildNIndex
        px = Particle.GetPx
        py = Particle.GetPy
        pz = Particle.GetPz
        pt = Particle.GetPt
        vert = Particle.GetVertex
        phi = Particle.GetPhi
        theta = Particle.GetTheta
        rap = Particle.GetRapidity
        eta = Particle.GetEta
        thetag = Particle.GetThetaVsGamma
        ptg = Particle.GetPtVsGamma
        e = Particle.GetE
        p = Particle.GetP
        m = Particle.GetM
        z = Particle.GetZ
        xf = Particle.GetXFeynman
        # Here are all the histogram Fill operations we need
        hks = self.KS.h.Fill
        horig = self.orig.h.Fill
        hchild = self.child.h.Fill
        hpypx = self.pypx.h.Fill
        hptpz = self.ptpz.h.Fill
        hvxy = self.vxy.h.Fill
        hthetaphi = self.thetaphi.h.Fill
        hetay = self.etay.h.Fill
        hgamma = self.gamma.h.Fill
        hem = self.em.h.Fill
        hzp = self.zp.h.Fill
        hxf = self.xf.h.Fill
        # We skip certain species for some of the histograms
        beams = [2212, 2112, 11]
        # Now we loop over the tracks
        for track in tracks:
            # Only fill final-state particles here as intermediate
            # state properties (e.g. negative mass for virtual photon)
            # can make it harder to spot issues.
            if track.GetStatus() != 1:
                continue
            hks(status(track))
            horig(parent(track))
            hchild(child1(track), childN(track))
            hpypx(px(track), py(track))
            hptpz(pz(track), pt(track))
            vertex = vert(track)
            hvxy(vertex.x(), vertex.y())
            hthetaphi(phi(track), theta(track))
            hetay(rap(track), eta(track))
            hgamma(thetag(track), ptg(track))
            # Skip particles that are not in the final state, or may
            # be scattered beams, as some quantities don't make sense then.
            if track.GetPdgCode().Code() not in beams:
                hem(m(track), e(track))
                hzp(p(track), z(track))
                hxf(xf(track))


def print_progress(n, max):
    print 'Event', str(n).rjust(len(str(max))), 'of', max

def execute(tree, outname, maxevents):
    '''
    '''
    eicpy.root_logon()
    outfile = ROOT.TFile(outname.replace('.ps', '.root'), 'recreate')
    h = EventHists()
    t = TrackHists()
    nevents = min(maxevents, tree.GetEntries())
    for i in xrange(0, nevents):
        tree.GetEntry(i)
        h.fill(tree.event)
        tracks = [track for track in tree.event.GetTracks()]
        t.fill(tracks)
        if (i + 1) % (nevents / 10) == 0:
            print_progress(i + 1, nevents)
    window = ROOT.TCanvas()
    window.Print(''.join([outname, '[']))
    h.output(window, outname)
    t.output(window, outname)
    window.Print(''.join([outname, ']']))
    outfile.Write()

if __name__ == '__main__':
    # Process the command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', help='Input ROOT file')
    parser.add_argument('outfile', help='Output PostScript file')
    parser.add_argument('--nevents', type=int, default=1000000000, help='Maximum number of events')
    args = parser.parse_args()
    file = ROOT.TFile(args.infile, 'read')
    execute(file.EICTree, args.outfile, args.nevents)
