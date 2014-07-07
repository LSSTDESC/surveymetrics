#!/usr/bin/env python
import numpy
import matplotlib.pyplot as plt

def rootsearch(f,a,b, band, dx):
    x1 = a; f1 = f(a, band)
    x2 = a + dx; f2 = f(x2, band)
    while f1*f2 > 0.0:
        if x1 >= b:
            return None,None
        x1 = x2; f1 = f2
        x2 = x1 + dx; f2 = f(x2, band)
    return x1,x2

def bisect(f,x1,x2,band,switch=0,epsilon=1.0e-9):
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

def roots(f, a, b, band, eps=1e-6):
    #print ('The roots on the interval [%f, %f] are:' % (a,b))
    ans=[]
    while 1:
        x1,x2 = rootsearch(f,a,b,band, eps)
        if x1 != None:
            a = x2
            root = bisect(f,x1,x2,band, 1)
            if root != None:
                ans.append(root)
                #print (round(root,-int(math.log(eps, 10))))
        else:
            #print ('\nDone')
            break
    return ans

class Interval:
    def __init__(self,a,b):
        self.start=a
        self.end=b

    def sum(self,x):
        return Interval(self.start+x,self.end+x)

    def length(self):
        return self.end-self.start

    def __str__(self):
        return '[{}, {}]\n'.format(self.start,self.end)
        
class Intervals:
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
        return Intervals(inter_)
        
    def length(self):
        ans=0
        for i in self.intervals:
            ans=ans+i.length()
        return ans

class ControlTimeMetric:

    def __init__(self,limmags, dates, bands, in_f,trange):

        self.magPrecision=0.01
        
        # the function
        self.in_f = in_f
        self.trange=trange
    
        self.limmags = limmags
        self.dates = dates
        self.bands =bands

        self.ranges=dict()

        self.calcRanges()
        
        #self.cutePlot()

    def calcControlTime(self):
        ans=0
        intervals=Intervals([])
        for l,b,d in zip(self.limmags,self.bands,self.dates):
            range= self.ranges[(numpy.round(l/self.magPrecision),b)]
            intervals.union(range.sum(-d).intervals)
        return intervals.length()
    
    def calcRanges(self):
        for l,b in zip(self.limmags,self.bands):
            if (numpy.round(l/self.magPrecision),b) not in self.ranges:
                self.ranges[(numpy.round(l/self.magPrecision),b)]=self.getDetectionRange(l,b)
        

    def cutePlot(self):       
        x_plt=numpy.arange(self.trange[0],self.trange[1],0.01)
        plt.plot(x_plt,self.in_f(x_plt,'r'))
        for l,b, i in zip(self.limmags,self.bands,numpy.linspace(0,1, len(self.limmags))):
            intervals = self.ranges[(numpy.round(l/self.magPrecision),b)]
            for interval in intervals.intervals:
                plt.plot(numpy.array([interval.start,interval.end]),numpy.array([l,l]),color=plt.cm.RdYlBu(i))
        plt.show()

        
    def getDetectionRange(self, limmag, band):
        f = lambda x, y : self.in_f(x,y)-limmag

        root=roots(f, self.trange[0], self.trange[1],band, eps=1e-2)
        if len(root) == 0:
            #check to see if the function is brighter than lim mag
            if f(numpy.sum(self.trange)/2,band) <= 0:
                #if so the entire range is detected
                return Intervals([Interval(self.trange[0], self.trange[1])])
            else:
                return Intervals([])

        root = numpy.append(root,self.trange)
        root = numpy.unique(root)

        #check to see if the first interval is brighter than limiting mag
        #print root[0], root[1], numpy.sum(root[0:2])/2, self.in_f(numpy.sum(root[0:2])/2,'r'), limmag,f(numpy.sum(root[0:1])/2,band)
        if f(numpy.sum(root[0:2])/2,band) <= 0:
            #if so the first range should be returned
            firstindex=0
        else:
            firstindex=1
        lastindex  = firstindex+ 2*((len(root)-firstindex)/2)

        ans = []
        for ind in xrange(firstindex,lastindex,2):
            ans.append(Interval(root[ind],root[ind+1]))
        
        return Intervals(ans)
