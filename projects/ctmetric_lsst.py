import numpy
import surveymetrics.ctmetric
import sncosmo
import lsst.sims.photUtils
from astropy import units as u
from astropy import coordinates as coord
import ephem

from lsst.sims.maf.metrics.baseMetric import BaseMetric

class ctmetric_lsst(BaseMetric):

    def __init__(self, function, range, metricName='ctmetric_lsst', magPrecision=0.01,
                 uniqueBlocks=False, **kwargs):
        """
            """
        
        cols=['mjd','filter','fivesigma_modified','fieldRA','fieldDec']
        super(snmetric_lsst, self).__init__(cols, metricName, **kwargs)
       
        self.metricDtype = 'float'
        
        #survey args
        self.filterinfo=dict()
        filterNames=numpy.array(['u','g','r','i','z','y'])
        filterWave = numpy.array([375.,476.,621.,754.,870.,980.])
        for name, mag, wave in zip(filterNames, filterTargetMag, filterWave):
            self.filterinfo[name]=dict()
            self.filterinfo[name]['wave']=wave
        
        #metric args
        self.function=function
        self.range = range
        self.magPrecision=magPrecision

        self.ccm=sncosmo.CCM89Dust()
        self.map = lsst.sims.photUtils.EBV.EbvMap()
        self.metric_calc=surveymetrics.ctmetric.ControlTimeMetric(self.function,self.range,magPrecision=self.magPrecision)

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
        ebv=sncosmo.get_ebv_from_map(co, mapdir='../data',cache=True)
        #np=ephem.Equatorial(dataSlice['fieldRA'][0],dataSlice['fieldDec'][0])
        #g=ephem.Galactic(np)
        #ebv=map.generateEbv(g.lon,g.lat)
        self.ccm.set(ebv=ebv)
        av=-2.5*numpy.log10(self.ccm.propagate(numpy.array([self.filterinfo[filter]['wave']*10]),numpy.array([1])))
        
#        print co, ccm.parameters, av
#        return self.filterinfo[filter]['targetMag']+av
        return self.filterinfo[filter]['targetMag']+av
