import numpy
import metric
import sncosmo
import lsst.sims.photUtils
from astropy import units as u
from astropy import coordinates as coord
import ephem

from lsst.sims.maf.metrics.baseMetric import BaseMetric

class SNMetric(BaseMetric):

    def __init__(self, metricName='SNMetric',  filterNames=numpy.array(['u','g','r','i','z','y']),
                 filterTargetMag= numpy.array([22,22.5,22.5,22.5,22.5,22]),filterWave = numpy.array([375.,476.,621.,754.,870.,980.]),
                 T0=30., Omega =1., zmax = .5 ,snr=50., tau=4., dt0=28., a=1.,
                 uniqueBlocks=False, **kwargs):
        """
            """
        
        cols=['night','filter','fivesigma_modified','fieldRA','fieldDec']
        super(SNMetric, self).__init__(cols, metricName, **kwargs)

        self.longestTolerableGap=15.
        
        self.metricDtype = 'float'
        #survey args
        self.filterinfo=dict()
        for name, mag, wave in zip(filterNames, filterTargetMag, filterWave):
            self.filterinfo[name]=dict()
            self.filterinfo[name]['targetMag']=mag
            self.filterinfo[name]['wave']=wave
        
        #self.filterNames = filterNames
        #self.filterTargetMag = filterTargetMag 
        #self.filterWave=filterWave
        
        #metric args
        # for now these are all scalars though in principle could be filter-dependent
        self.T0=T0
        self.Omega=Omega
        self.zmax=zmax
        self.snr=snr
        self.tau=tau
        self.dt0=dt0
        self.a=a

        self.ccm=sncosmo.CCM89Dust()
        self.map = lsst.sims.photUtils.EBV.EbvMap()
        self.metric_calc=metric.Metric.OneFieldBandMetric(self.Omega,self.zmax,self.snr,self.T0,self.tau,self.dt0,self.a)

    def run(self, dataSlice):
        ans=0.
        #Figure out the seasons for the pixel
        seasons= self.splitBySeason(dataSlice)
        for s in seasons:
            ans=ans+self.seasonMetric(s, dataSlice)
        #return {'result':ans}
        return ans
                           
    def splitBySeason(self,dataSlice):
        dates= numpy.unique(dataSlice['night'])
        gaps = dates - numpy.roll(dates,1)
        wheretosplit = numpy.where(gaps > self.longestTolerableGap)[0]
        wheretosplit= numpy.append(wheretosplit,[0,len(dates)-1])
        wheretosplit=numpy.unique(wheretosplit)
        out=[]
        for i in xrange(len(wheretosplit)-1):
            out.append(numpy.array(dates[wheretosplit[i]:wheretosplit[i+1]]))
        return out

    def seasonMetric(self, s, dataSlice):
        ans=0
        for filter in self.filterinfo.keys():
            night_in=[]
            snr_in=[]
            targetMag=self.getTargetMag(dataSlice, filter)           
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

    def getTargetMag(self, dataSlice, filter):
    # fainten the extragalactic magnitude due to Galactic dust
        dec=numpy.maximum(dataSlice['fieldDec'][0],-numpy.pi/2)
        co = coord.ICRS(ra=dataSlice['fieldRA'][0], dec=dec, unit=(u.rad, u.rad))
        ebv=sncosmo.get_ebv_from_map(co, cache=True)
        #np=ephem.Equatorial(dataSlice['fieldRA'][0],dataSlice['fieldDec'][0])
        #g=ephem.Galactic(np)
        #ebv=map.generateEbv(g.lon,g.lat)
        self.ccm.set(ebv=ebv)
        av=-2.5*numpy.log10(self.ccm.propagate(numpy.array([self.filterinfo[filter]['wave']*10]),numpy.array([1])))
        
#        print co, ccm.parameters, av
#        return self.filterinfo[filter]['targetMag']+av
        return self.filterinfo[filter]['targetMag']+av
