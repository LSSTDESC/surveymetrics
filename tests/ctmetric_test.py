#!/usr/bin/env python
import numpy
import unittest
import surveymetrics.ctmetric
import matplotlib.pyplot as plt
import sncosmo
from astropy.cosmology import WMAP9 as cosmo


class TestMetric(unittest.TestCase):

    def test_sn(self):
        z=0.1
        model = sncosmo.Model(source='salt2-extended')
        model.set(z=z, t0=55000.)
        band = sncosmo.get_bandpass('desg')
        model.set_source_peakabsmag(-19.3,band,'ab')

        def fn(x,y):
            return model.bandmag(y,'ab', x)

        metric=surveymetrics.ctmetric.ControlTimeMetric(fn,numpy.array([55000.,55010]))

        print metric.calcControlTime(numpy.array([20., 20.]),numpy.array([0.,10]),['desr','desr'])
        
    def test_0(self):

        def fn(x,c):
            return -numpy.sin(x)
        metric=surveymetrics.ctmetric.ControlTimeMetric(fn,numpy.array([-10.5,13.5]))
    

        def fn(x,c):
            return numpy.sin(x)
            
        metric=surveymetrics.ctmetric.ControlTimeMetric(fn,numpy.array([-10.5,13.5]))

    def test_1(self):

        def fn(x,c):
            return 100.* ((x<0.) | (x>1.))
            
        metric=surveymetrics.ctmetric.ControlTimeMetric(fn,numpy.array([-.5,1.5]))
        self.assertTrue(numpy.abs(1.5-metric.calcControlTime(numpy.array([1., 1.]),numpy.array([0.,.5]),['r','r'])) < 1e-8)
             
        self.assertTrue(numpy.abs(1.-metric.calcControlTime(numpy.array([1.]),numpy.array([0.]),['r'])) < 1e-8)

        self.assertTrue(numpy.abs(2.-metric.calcControlTime(numpy.array([1., 1.]),numpy.array([0.,10]),['r','r'])) < 1e-8)

        metric=surveymetrics.ctmetric.ControlTimeMetric(fn,numpy.array([-.5,1.5]))
        self.assertTrue(numpy.abs(1.-metric.calcControlTime(numpy.array([1., -1.]),numpy.array([0.,10]),['r','r'])) < 1e-8)


    def test_2(self):

        def fn(x,c):
            return 100.* ((x < 0) | ((x>1) & (x < 2)) | (x>3))
            
        metric=surveymetrics.ctmetric.ControlTimeMetric(fn,numpy.array([-.5,3.5]))
        self.assertTrue(numpy.abs(3-metric.calcControlTime(numpy.array([1., 1.]),numpy.array([0.,.5]),['r','r'])) < 1e-8)
             
        self.assertTrue(numpy.abs(2.-metric.calcControlTime(numpy.array([1.]),numpy.array([0.]),['r'])) < 1e-8)

        self.assertTrue(numpy.abs(4.-metric.calcControlTime(numpy.array([1., 1.]),numpy.array([0.,10]),['r','r'])) < 1e-8)

        self.assertTrue(numpy.abs(2.-metric.calcControlTime(numpy.array([1., -1.]),numpy.array([0.,10]),['r','r'])) < 1e-8)

        
if __name__ == '__main__':
    unittest.main()
