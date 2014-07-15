import numpy
import surveymetrics.ctmetric
import lsst_info
import sncosmo
import lsst.sims.photUtils
from astropy import units as u
from astropy import coordinates as coord
import ephem

from lsst.sims.maf.metrics.baseMetric import BaseMetric

class ctmetric_lsst(BaseMetric):

    def __init__(self, function, range, metricName='ctmetric_lsst', magPrecision=0.01,
                 limitingsnr=5.,
                 uniqueBlocks=False, **kwargs):
        """
            """

        self.metriccache=dict()
        self.ebvprecision=100.
        self.limitingsnr=limitingsnr
        
        cols=['mjd','filter','fivesigma_modified','fieldRA','fieldDec']
        super(snmetric_lsst, self).__init__(cols, metricName, **kwargs)
       
        self.metricDtype = 'float'
        
        #survey args
        self.filterinfo=lsst.LSST.filterinfo
        
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
                           

    def run(self, dataSlice):
       #get ebv of this slice
        dec=numpy.maximum(dataSlice['fieldDec'][0],-numpy.pi/2)
        co = coord.ICRS(ra=dataSlice['fieldRA'][0], dec=dec, unit=(u.rad, u.rad))
        ebv=sncosmo.get_ebv_from_map(co, mapdir='../data',cache=True)

        #check to see of metric is already available for this ebv
        label=round(ebv/ebvprecision)
        if label in self.metriccache:
            #if so use it
            metric = self.metriccache[label]
        else:
            #if not create it
            ##### To fix with magnitude later
            fn = lambda x, y : self.function(x,y)
            metric = surveymetric.ctmetric.ControlTimeMetric(fn ,self.range, self.magPrecision)
            self.metriccache[label]=metric

        #convert limiting magnitude to threshold assuming sky limit
        m5col_ = dataSlice['fivesigma_modified']
        mtarget = m5col_-2.5*numpy.log10(self.limitingsnr/5.)
        mjds = dataSlice['mjd']
        bands= dataSlice['bands']

        #calculate
        return metric.calcControlTime(mtarget, dates, bands)
