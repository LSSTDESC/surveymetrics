#!/usr/bin/env python
import numpy
import unittest
import surveymetrics.ctmetric
import matplotlib.pyplot as plt

class TestMetric(unittest.TestCase):

        
    def test_0(self):

        def fn(x,c):
            return -numpy.sin(x)
        metric=surveymetrics.ctmetric.ControlTimeMetric(numpy.array([.5, .5]),numpy.array([0.,.5]),['r','r'],fn,numpy.array([-10.5,13.5]))
    

        def fn(x,c):
            return numpy.sin(x)
            
        metric=surveymetrics.ctmetric.ControlTimeMetric(numpy.array([.5, .5]),numpy.array([0.,.5]),['r','r'],fn,numpy.array([-10.5,13.5]))

    def test_1(self):

        def fn(x,c):
            return 100.* ((x<0.) | (x>1.))
            
        metric=surveymetrics.ctmetric.ControlTimeMetric(numpy.array([1., 1.]),numpy.array([0.,.5]),['r','r'],fn,numpy.array([-.5,1.5]))
        self.assertTrue(numpy.abs(1.5-metric.calcControlTime()) < 1e-8)
             
        metric=surveymetrics.ctmetric.ControlTimeMetric(numpy.array([1.]),numpy.array([0.]),['r'],fn,numpy.array([-.5,1.5]))
        self.assertTrue(numpy.abs(1.-metric.calcControlTime()) < 1e-8)

        metric=surveymetrics.ctmetric.ControlTimeMetric(numpy.array([1., 1.]),numpy.array([0.,10]),['r','r'],fn,numpy.array([-.5,1.5]))
        self.assertTrue(numpy.abs(2.-metric.calcControlTime()) < 1e-8)

        metric=surveymetrics.ctmetric.ControlTimeMetric(numpy.array([1., -1.]),numpy.array([0.,10]),['r','r'],fn,numpy.array([-.5,1.5]))
        self.assertTrue(numpy.abs(1.-metric.calcControlTime()) < 1e-8)


    def test_2(self):

        def fn(x,c):
            return 100.* ((x < 0) | ((x>1) & (x < 2)) | (x>3))
            
        metric=surveymetrics.ctmetric.ControlTimeMetric(numpy.array([1., 1.]),numpy.array([0.,.5]),['r','r'],fn,numpy.array([-.5,3.5]))
        self.assertTrue(numpy.abs(3-metric.calcControlTime()) < 1e-8)
             
        metric=surveymetrics.ctmetric.ControlTimeMetric(numpy.array([1.]),numpy.array([0.]),['r'],fn,numpy.array([-.5,3.5]))
        self.assertTrue(numpy.abs(2.-metric.calcControlTime()) < 1e-8)

        metric=surveymetrics.ctmetric.ControlTimeMetric(numpy.array([1., 1.]),numpy.array([0.,10]),['r','r'],fn,numpy.array([-.5,3.5]))
        self.assertTrue(numpy.abs(4.-metric.calcControlTime()) < 1e-8)

        metric=surveymetrics.ctmetric.ControlTimeMetric(numpy.array([1., -1.]),numpy.array([0.,10]),['r','r'],fn,numpy.array([-.5,3.5]))
        self.assertTrue(numpy.abs(2.-metric.calcControlTime()) < 1e-8)

        
if __name__ == '__main__':
    unittest.main()
