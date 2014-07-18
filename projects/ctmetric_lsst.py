import numpy
import surveymetrics.ctmetric
import lsst_info
import sncosmo
import lsst.sims.photUtils
from astropy import units as u
from astropy import coordinates as coord
import ephem
import importlib
import string

from lsst.sims.maf.metrics.baseMetric import BaseMetric

class ctmetric_lsst(BaseMetric):
    
    def __init__(self,  metricName='ctmetric_lsst', magPrecision=0.01,
                 limitingsnr=5.,
                 uniqueBlocks=False, **kwargs):
        """
            """
        cols=['expMJD','filter','fivesigma_modified'] #,'fieldRA','fieldDec']

        super(ctmetric_lsst, self).__init__(cols, metricName, **kwargs)

        self.metriccache=dict()
        self.ebvprecision=10.
        self.limitingsnr=limitingsnr
        
       
        self.metricDtype = 'float'
        
        #survey args
        self.filterinfo=lsst_info.LSST.filterinfo
        
        #metric args
        self.magPrecision=magPrecision

        self.ccm=sncosmo.CCM89Dust()
        self.map = lsst.sims.photUtils.EBV.EbvMap()

        dotlocation = string.rfind(kwargs['lc.class'],'.')
        module_name = kwargs['lc.class'][0:dotlocation]
        class_name = kwargs['lc.class'][dotlocation+1:]
        somemodule = importlib.import_module(module_name)
        
        self.lightcurve =getattr(somemodule, class_name)(**kwargs)
                           

    def run(self, dataSlice, slicePoint=None):
       #get ebv of this slice
#        dec=numpy.maximum(dataSlice['fieldDec'][0],-numpy.pi/2)
#        co = coord.ICRS(ra=dataSlice['fieldRA'][0], dec=dec, unit=(u.rad, u.rad))

        co = coord.ICRS(ra=slicePoint['ra'], dec=slicePoint['dec'], unit=(u.rad, u.rad))
        
#        ebv=sncosmo.get_ebv_from_map(co, mapdir='../data',cache=True)
        ebv=0.
        
        #check to see of metric is already available for this ebv
        label=round(ebv/self.ebvprecision)
        if label in self.metriccache:
            #if so use it
            metric = self.metriccache[label]
        else:
            #if not create it
            fn = lambda x,y: self.lightcurve.mag(x,y)
            
            ##### To fix with magnitude later
            metric = surveymetrics.ctmetric.ControlTimeMetric(fn,self.lightcurve.trange,magPrecision=self.magPrecision)
            self.metriccache[label]=metric

        #convert limiting magnitude to threshold assuming sky limit
        m5col_ = dataSlice['fivesigma_modified']
        mtarget = m5col_-2.5*numpy.log10(self.limitingsnr/5.)
        mjds = dataSlice['expMJD']
        bands= dataSlice['filter']

        #calculate
        return metric.calcControlTime(mtarget, mjds, bands)
