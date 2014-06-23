#!/usr/bin/env python
import numpy
import unittest
import metric
#import metric_des
import matplotlib.pyplot as plt

class TestMetric(unittest.TestCase):

    def setUp(self):
        self.desfields=['SN-E1','SN-E2','SN-X1','SN-X2','SN-X3','SN-C1','SN-C2','SN-C3','SN-S1','SN-S2']
        self.longitude = 1.23595297054
        
    # test reading an astac file with known contents
    def test_read_astac(self):
        file='../data/schedule15/obssim-19/obs.txt'
        obs=metric.Metric.read_astac(file)
        self.assertEquals(len(obs), 22998)
        self.assertEquals(obs[0]['start_mjd'], '56888.19921573579')
        self.assertEquals(obs[22997]['hexname'], '1091-247')
        

    def test_prune_by_fieldname(self):
        file='../data/schedule15/obssim-19/obs.txt'

        obs=metric.Metric.read_astac(file)
        metric.Metric.prune_by_fieldname(obs,self.desfields)
        self.assertEquals(len(obs), 3206)
        obs = metric.Metric.get_by_band(obs,'i')
        self.assertEquals(len(obs),626)
        obs = metric.Metric.get_by_fieldname(obs,'SN-X3')
        self.assertEquals(len(obs),162)
        groups=metric.Metric.group_by_night(obs,self.longitude)
        
    # a simple case with an almost analytic answer
    def test_onefieldband(self):
        omega=1.
        zmax=0
        t0=0
        tau=4
        dt0=4
        snr_0=20.
        a=1.
        
        ndata=20.
        epoch1=numpy.arange(0,tau*(ndata+1),tau)

        t =epoch1[-1]-epoch1[0]
        
        ston1=numpy.zeros(ndata)+snr_0
        
        onemetric= metric.Metric.OneFieldBandMetric(omega,zmax,snr_0, t0, tau, dt0, a)

        #target SNR
        metric_1 = onemetric.metric(ston1,epoch1)
        ans=omega*ndata*tau*(1-numpy.exp(-1))#/  (tau/dt0*(1-numpy.exp(-dt0/tau)))
        self.assertTrue(numpy.abs(ans-metric_1) < tau)

        #lower SNR
        metric_1 = onemetric.metric(ston1/2,epoch1)
        self.assertTrue(numpy.abs(ans/4-metric_1) < tau)

        #higher SNR
        metric_1 = onemetric.metric(ston1*2,epoch1)
        self.assertTrue(numpy.abs(ans-metric_1) < tau)

        #higher SNR randomized SNR
        ndata2=10000.
        epoch2=numpy.arange(0,tau*(ndata2+1),tau)
        t2 =epoch2[-1]-epoch2[0]
        ston2= numpy.random.normal(loc=snr_0/2,scale=1,size=ndata2)
        metric_1 = onemetric.metric(ston2,epoch2)
        ans=omega*t2*(tau/t2*ndata2/4.*(1-numpy.exp(-1))-a/snr_0**2 * (tau/dt0*(1-numpy.exp(-dt0/tau))))
        self.assertTrue(numpy.abs(ans-metric_1) < 0.02*ans)


        
                
    def test_onefieldband_plot(self):
        tau=4.
        dt0=4.
        snr_0=20.
        ndata=20
        epoch1=numpy.arange(0,tau*(ndata+1),tau)
        ston1=numpy.zeros(ndata)+snr_0
        
        onemetric= metric.Metric.OneFieldBandMetric(1.,1.,snr_0, 30, tau, dt0,1)
        metric_1 = onemetric.metric(ston1,epoch1)

        metrics_2=[]
        for i in xrange(200):
            epoch2 = numpy.arange(1,tau*ndata+-1)
            numpy.random.shuffle(epoch2)
            epoch2=epoch2[0:ndata-1]
            epoch2=numpy.append(epoch2,[0,tau*ndata])
            epoch2=numpy.sort(epoch2)
            metrics_2.append(onemetric.metric(ston1,epoch2))
        metrics_2=numpy.array(metrics_2)

        metrics_3=[]
        for i in xrange(200):
            ston2 = numpy.random.normal(loc=20,scale=1,size=ndata)
            metrics_3.append(onemetric.metric(ston2,epoch1))
        metrics_3=numpy.array(metrics_3)
        plt.hist([metrics_2,metrics_3],label=['Random Phase','Random S/N'])
        plt.axvline(metric_1,linewidth=4,label='Uniform',color='black')
        plt.xlabel('Metric')
        plt.legend(loc=2)
        plt.savefig('fom.eps')
        self.assertTrue(True)

        plt.clf()
        for i in xrange(len(epoch2)-1):
            x=numpy.arange(epoch2[i],epoch2[i+1],0.02)
            plt.fill_between(x,ston2[i]*numpy.exp(-(x-epoch2[i])/4.),numpy.zeros(len(x)))
        plt.xlabel('date')
        plt.ylabel('SNR$^2$')
        plt.savefig('saw.eps')

if __name__ == '__main__':
    unittest.main()
