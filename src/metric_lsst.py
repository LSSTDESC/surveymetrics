#!/usr/bin/env python
import numpy
from lsst.sims.operations.maf.metrics.simpleMetrics import SimpleScalarMetric
from metric import Metric

desfields=['SN-E1','SN-E2','SN-X1','SN-X2','SN-X3','SN-C1','SN-C2','SN-C3','SN-S1','SN-S2']
dessite='TOLOLO4M'
desetime=dict()
desetime['shallow']=dict()
desetime['shallow']['g']=175.
desetime['shallow']['r']=150.
desetime['shallow']['i']=200.
desetime['shallow']['z']=2*200.
desetime['shallow']['zmax']=0.6
desetime['deep']=dict()
desetime['deep']['g']=3*200.
desetime['deep']['r']=3*400.
desetime['deep']['i']=5*360.
desetime['deep']['z']=11*330.
desetime['deep']['zmax']=0.8
longitude = 1.23595297054

class SNMetricForLSST(SimpleScalarMetric):

    def __init__(self,input):
        dum = 0
        
    def run(self, dataSlice):
        iqr = np.percentile(dataSlice[self.colname],75)-np.percentile(dataSlice[self.colname],25)
        rms = iqr/1.349 #approximation
        return rms

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser( description = "Remove targets close to database" )
    parser.add_argument( "filename"       , help = "input filename with target coadd_objects_id"   )
    args = parser.parse_args()
    pdict=vars(args)

    dum = SNMetricForLSST('blah')

#./metric_des  ../data/schedule15/obssim-19/obs.txt
