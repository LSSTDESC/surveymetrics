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

#f=lambda x, y:x*numpy.cos(x-4) +1
#print roots(f, -3, 3, 'b')

#wefew
class CTMetric:

    def __init__(self):

        # the function
        self.in_f = lambda x, y :x*numpy.cos(x-4)

        self.trange=numpy.array([0,10]) # range of applicability of the function
    
        self.limmags = numpy.array([-3,-1,1,3])
        self.dates = numpy.array([1,3,7,8])
        self.bands = ['r','r','r','r']

        x_plt=numpy.arange(self.trange[0],self.trange[1],0.01)
        for l,b in zip(self.limmags,self.bands):
            ranges = self.getNoDetectionRange(l,b)
            plt.plot(x_plt,self.in_f(x_plt,'r'))
            for range in ranges:
                plt.plot(range,numpy.array([l,l]))
            plt.show()

        
    def getNoDetectionRange(self, limmag, band):
        f = lambda x, y : self.in_f(x,y)-limmag

        #print self.in_f(2,'r')
        #print f(2,'r')

        root=roots(f, self.trange[0], self.trange[1],band, eps=1e-2)
        if len(root) == 0:
            #check to see if the function is fainter than lim mag
            if f(numpy.sum(trange)/2,band) > 0:
                #if so the entire range is not detected
                return trange
            else:
                return None

        root = numpy.append(root,self.trange)
        root = numpy.unique(root)

        #check to see if the first interval is fainter than limiting mag
        if f(numpy.sum(root[0:1])/2,band) > 0:
            #if so the first range should be returned
            firstindex=0
        else:
            firstindex=1

        ans = []

        for ind in xrange(firstindex,len(root),2):
            ans.append(numpy.array([root[ind],root[ind+1]]))
            #print ans[-1]
        return ans

shit = CTMetric()
