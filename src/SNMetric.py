import numpy
import metric

from lsst.sims.maf.metrics.baseMetric import BaseMetric

class SNMetric(BaseMetric):

    def __init__(self, metricName='SNMetric', nightcol='night', filtercol='filter',
                 m5col='fivesigma_modified', units='', 
                 T0=30., 
                 uniqueBlocks=False, **kwargs):
        """
            """
        
        cols=[nightcol,m5col,filtercol]
        self.nightcol = nightcol
        self.m5col = m5col
        self.filtercol = filtercol
        super(SNMetric, self).__init__(cols, metricName, units=units, **kwargs)
        self.filterNames = numpy.array(['u','g','r','i','z','y'])
        self.filterTargetMag = numpy.array([22,22.5,22.5,22.5,22.5,22]) # temporary need separate for Wide and DD
        self.T0=T0
        self.Omega=1.
        self.zmax=.5
        self.snr=20
        self.tau=4.
        self.dt0=28.
        self.a=kwargs['solidangle']
        print self.a
        self.metric_calc=metric.Metric.OneFieldBandMetric(self.Omega,self.zmax,self.snr,self.T0,self.tau,self.dt0,self.a)
        
    def run(self, dataSlice):
        ans=0
        #Figure out the seasons for the pixel
        seasons= SNMetric.splitBySeason(dataSlice,self.T0)
        for s in seasons:
            ans=ans+self.seasonMetric(s, dataSlice)
        return {'result':ans}
                               

    @staticmethod
    def splitBySeason(dataSlice,T0):
        dates= numpy.unique(dataSlice['night'])
        gaps = dates - numpy.roll(dates,1)
        wheretosplit = numpy.where(gaps > T0)[0]
        wheretosplit= numpy.append(wheretosplit,[0,len(dates)-1])
        wheretosplit=numpy.unique(wheretosplit)
        out=[]
        for i in xrange(len(wheretosplit)-1):
            out.append(numpy.array(dates[wheretosplit[i]:wheretosplit[i+1]]))
        return out

    def seasonMetric(self, s, dataSlice):
        ans=0
        for filter, targetMag in zip(self.filterNames,self.filterTargetMag):
            night_in=[]
            snr_in=[]
            for night in s:
                m5col_ = dataSlice['fivesigma_modified'][numpy.logical_and(dataSlice['night']==night,dataSlice['filter'] == filter)]
                if len(m5col_) > 0:
                    night_in.append(night)
                    snrs_=5*10**((-targetMag+m5col_)/2.5)
                    snr_=numpy.sqrt(numpy.sum(numpy.power(snrs_,2)))
                    snr_in.append(snr_)
            if len(night_in) == 0:
                return 0
            night_in.append(night_in[-1]+self.tau)
            nignt_in=numpy.array(night_in, dtype=numpy.float)
            snr_in=numpy.array(snr_in)
            ans= ans +self.metric_calc.metric(snr_in,night_in)
        return ans
