""" Control Time Calculation

Summary
-------
This module is used to calculate the control time for a source with
input observer light curve given the limiting detection magnitudes for
a series of observations.  The control time is the time duration over
which a new occurance of a transient will be detected.

Given a limiting magnitude the phases at which the input light curve
can be detected are determined.  When associated with a date of
observation, the phases translate to starting dates of transients that
can be detected.  The union of these dates for all observations is
used to calculate the control time.

The phases at which the input light curve can be detected is
determined through finding roots.  This is done through the three
functions _rootsearch, _bisect, and _roots.  The use of a standard
package may produce more robust and faster root finding.  A failure
mode is the algorithm not finding all the roots.  Performance could be
enhanced for certain functional forms of the light curves that allow
analytic determination of roots.

The class _Interval is used to contain a single phase and time
interval.  The class _Intervals is used to contain sets of _Intervals
that are non-overlapping.

Each filter/limiting magnitude pair has its own corresponding
discovery phases.  In a large survey, the limiting magnitudes of the
observations can be very similar.  In order to save compute time, the
code uses previously-computed discovery phases from another limiting
magnitude within a user-controlled precision.  The class
ControlTimeMetric is designed for the calculation of control times for
a single set of light curves.  It chaches discovery phases.  The
method ControlTimeMetric.calcControlTime calculates or retrieves
discovery phases for a set of input observations and returns the
control time.

The implementation is set up to be efficient by having one ControlTimeMetric object per input
light curve.  That one object calculates the control time for
different surveys.

"""

import numpy
import matplotlib.pyplot as plt

def _rootsearch(f,a,b, band, dx):
    x1 = a; f1 = f(a, band)
    x2 = a + dx; f2 = f(x2, band)
    while f1*f2 > 0.0:
        if x1 >= b:
            return None,None
        x1 = x2; f1 = f2
        x2 = x1 + dx; f2 = f(x2, band)
    return x1,x2

def _bisect(f,x1,x2,band,switch=0,epsilon=1.0e-9):
    f1 = f(x1,band)
    if f1 == 0.0:
        return x1
    f2 = f(x2,band)
    if f2 == 0.0:
        return x2
    if f1*f2 > 0.0:
        print('Root is not bracketed')
        return None
    n = numpy.ceil(numpy.log(abs(x2 - x1)/epsilon)/numpy.log(2.0))
    for i in range(int(n)):
        x3 = 0.5*(x1 + x2); f3 = f(x3,band)
        if (switch == 1) and (abs(f3) >abs(f1)) and (abs(f3) > abs(f2)):
            return None
        if f3 == 0.0:
            return x3
        if f2*f3 < 0.0:
            x1 = x3
            f1 = f3
        else:
            x2 =x3
            f2 = f3
    return (x1 + x2)/2.0

def _roots(f, a, b, band, eps=1e-6):
    ans=[]
    while 1:
        x1,x2 = _rootsearch(f,a,b,band, eps)
        if x1 != None:
            a = x2
            root = _bisect(f,x1,x2,band, 1)
            if root != None:
                ans.append(root)
        else:
            break
    return ans

class _Interval:
    def __init__(self,a,b):
        self.start=a
        self.end=b

    def sum(self,x):
        return _Interval(self.start+x,self.end+x)

    def length(self):
        return self.end-self.start

    def __str__(self):
        return '[{}, {}]\n'.format(self.start,self.end)
        
class _Intervals:
    def __init__(self,interval_list):
        self.intervals=[]
        if len(interval_list) > 0:
            self.union(interval_list)

    def __str__(self):
        out = '['
        first = True
        for i in self.intervals:
            if first:
                comma=''
                first = False
            else:
                comma=','
            out = out+'{} [{}, {}]'.format(comma, i.start,i.end)
        out= out + ']\n'
        return out

    def union(self,interval_list):
        allIntervals=[]
        if len(self.intervals) > 0:
            allIntervals.extend(self.intervals)
        if len(interval_list) > 0:
            allIntervals.extend(interval_list)
        if len(allIntervals) > 0:
            allIntervals.sort(key = lambda self: self.start)
            y=[ allIntervals[0] ]
            for x in allIntervals[1:]:
                if y[-1].end < x.start:
                    y.append(x)
                else:
                    y[-1].end = max(y[-1].end,x.end)
            self.intervals=y

    def sum(self, x):
        inter_=[]
        for i in self.intervals:
            inter_.append(i.sum(x))
        return _Intervals(inter_)
        
    def length(self):
        ans=0
        for i in self.intervals:
            ans=ans+i.length()
        return ans

class ControlTimeMetric:

    """ blah blah
    """
    def __init__(self, in_f,trange, magPrecision=0.01):

        self.magPrecision=magPrecision
        self.ranges=dict()
        
        # the function and range
        self.in_f = in_f
        self.trange=trange

    def calcControlTime(self, limmags, dates, bands):

        self._calcRanges(limmags, bands)
        ans=0
        intervals=_Intervals([])
        for l,b,d in zip(limmags,bands,dates):
            range= self.ranges[(numpy.round(l/self.magPrecision),b)]
            intervals.union(range.sum(-d).intervals)
        return intervals.length()
    
    def _calcRanges(self, limmags, bands):
        for l,b in zip(limmags,bands):
            if (numpy.round(l/self.magPrecision),b) not in self.ranges:
                self.ranges[(numpy.round(l/self.magPrecision),b)]=self.getDetectionRange(l,b)
        

    def cutePlot(self, limmags, bands):       
        x_plt=numpy.arange(self.trange[0],self.trange[1],0.01)
        plt.plot(x_plt,self.in_f(x_plt,'r'))
        for l,b, i in zip(limmags,bands,numpy.linspace(0,1, len(limmags))):
            intervals = self.ranges[(numpy.round(l/self.magPrecision),b)]
            for interval in intervals.intervals:
                plt.plot(numpy.array([interval.start,interval.end]),numpy.array([l,l]),color=plt.cm.RdYlBu(i))
        plt.show()

        
    def getDetectionRange(self, limmag, band):
        f = lambda x, y : self.in_f(x,y)-limmag

        root=_roots(f, self.trange[0], self.trange[1],band, eps=1e-2)
        if len(root) == 0:
            #check to see if the function is brighter than lim mag
            if f(numpy.sum(self.trange)/2,band) <= 0:
                #if so the entire range is detected
                return _Intervals([_Interval(self.trange[0], self.trange[1])])
            else:
                return _Intervals([])

        root = numpy.append(root,self.trange)
        root = numpy.unique(root)

        #check to see if the first interval is brighter than limiting mag
        if f(numpy.sum(root[0:2])/2,band) <= 0:
            #if so the first range should be returned
            firstindex=0
        else:
            firstindex=1
        lastindex  = firstindex+ 2*((len(root)-firstindex)/2)

        ans = []
        for ind in xrange(firstindex,lastindex,2):
            ans.append(_Interval(root[ind],root[ind+1]))
        
        return _Intervals(ans)
