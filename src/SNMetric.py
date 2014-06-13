import numpy
import metric

from lsst.sims.maf.metrics.baseMetric import BaseMetric

class SNMetric(BaseMetric):

    def __init__(self, metricName='SNMetric',  filterNames=numpy.array(['u','g','r','i','z','y']),
                 filterTargetMag= numpy.array([22,22.5,22.5,22.5,22.5,22]),
                 T0=30., Omega =1., zmax = .5 ,snr=50., tau=4., dt0=28., a=1.,
                 uniqueBlocks=False, **kwargs):
        """
            """
        
        cols=['night','filter','fivesigma_modified']
        super(SNMetric, self).__init__(cols, metricName, **kwargs)

        #survey args
        self.filterNames = filterNames
        self.filterTargetMag = filterTargetMag 

        #metric args
        # for now these are all scalars though in principle could be filter-dependent
        self.T0=T0
        self.Omega=Omega
        self.zmax=zmax
        self.snr=snr
        self.tau=tau
        self.dt0=dt0
        self.a=a
        
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
                    
            #if a filter is not measured in a season the metric for the season is zero
            if len(night_in) == 0:
                return 0

            #calculate the metric
            night_in.append(night_in[-1]+self.tau)
            nignt_in=numpy.array(night_in, dtype=numpy.float)
            snr_in=numpy.array(snr_in)
            ans= ans +self.metric_calc.metric(snr_in,night_in)
        return ans
