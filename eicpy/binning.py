import array
import math

import ROOT

class Bins(object):
    '''
    Base class for generating bin boundaries in some range.
    '''
    def __init__(self, n, lower, upper):
        '''
        Initialise with the number of ranges for which
        to generate bin edges, in the range [lower, upper].
        '''
        self.n = n
        self.lower = lower
        self.upper = upper


class BinsLog10(Bins):
    '''
    Generator class for bins with equal intervals in log10.
    '''
    def __init__(self, n, lower, upper):
        super(BinsLog10, self).__init__(n, lower, upper)
        self.width = math.log10(self.upper / self.lower) / self.n
        # Redefine lower and upper as log10 equivalents
        self.lower = math.log10(self.lower)
        self.upper = math.log10(self.upper)

    def edges(self):
        '''
        Generator function for (n+1) bin edges.
        '''
        for i in range(0, self.n + 1):
            yield math.pow(10., self.lower + i * self.width)

    def lower_edges(self):
        '''
        Generator function for n bin lower edges.
        '''
        for i in range(0, self.n):
            yield math.pow(10., self.lower + i * self.width)

    def upper_edges(self):
        '''
        Generator function for n bin upper edges.
        '''
        for i in range(1, self.n + 1):
            yield math.pow(10., self.lower + i * self.width)


def bin_log10(obj):
    '''
    Rebins an axis such that bins span equal intervals in log_10(x).
    '''
    if isinstance(obj, ROOT.TAxis):
        binner = BinsLog10(obj.GetNbins(), obj.GetXmin(), obj.GetXmax())
        # Need array('d') to pass to TAxis.Set() in place of Double_t*
        edgearray = array.array('d', [x for x in binner.edges()])
        obj.Set(obj.GetNbins(), edgearray)
    elif isinstance(obj, ROOT.TH1):
        axes = [obj.GetXaxis(), obj.GetYaxis(), obj.GetZaxis()]
        for i in range(0, obj.GetDimension()):
            bin_log10(axes[i])
