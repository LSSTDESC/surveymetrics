#!/usr/bin/env python
import numpy
#from pyslalib import slalib

class Metric:

    #the metric in the note
    class OneFieldBandMetric:

        #the metric parameters
        def __init__(self, omega, zmax,snr_0, T0, tau, dt0, a):
            self.omega=omega
            self.zmax=zmax
            self.T0=T0
            self.snr_0=snr_0
            self.dt0=dt0
            self.tau=tau
            self.a=a

        #the method that calculates the metric
        #ston and epoch must be sorted according to epoch
        def metric(self, ston, epoch):
            T=epoch[-1]-epoch[0]
            deltat=numpy.roll(epoch,-1)-epoch
            deltat=deltat[0:len(deltat)-1]
            ans=self.omega*(T-(1+self.zmax)*self.T0)
            if (ans <= 0):
		return 0
            ans = ans * (self.tau/T*numpy.sum(numpy.minimum(ston**2/self.snr_0**2,1)*(1-numpy.exp(-deltat/self.tau)))-self.a*numpy.var(numpy.minimum(ston,self.snr_0))/self.snr_0**2)
            return ans / (self.tau/self.dt0*(1-numpy.exp(-self.dt0/self.tau))) # normalize to ideal survey

    #def __init__(self):
        #initialize dictionary of OneFieldBandMetrics

    @staticmethod
    def prune_by_fieldname(obs,names):
        for i in xrange(len(obs)-1,-1,-1):
            if obs[i]['hexname'] not in names:
                obs.pop(i)

    @staticmethod
    def get_by_fieldname(obs,field):
        out=[]
        for i in xrange(len(obs)):
            if obs[i]['hexname'] == field:
                out.append(obs[i])
        return out

    @staticmethod
    def get_by_band(obs,band):
        out=[]
        for i in xrange(len(obs)):
            if obs[i]['filter'] == band:
                out.append(obs[i])
        return out

    @staticmethod
    def group_by_night(obs,long):
#        long = slalib.sla_obs(0, site)[2]
        localnoon=long/2/numpy.pi # when noon at site time at GMT
        #a single night is divided into days split by noon
        day=numpy.array([float(o['start_mjd']) for o in obs])
        day = day - localnoon -0.5
        day = numpy.floor(day)
        o_and_d=zip(obs,day)
        
        uniqueday = sorted(set(day))
        
        ans=dict()
        for uday in uniqueday:
            dum=[]
            for o,d in o_and_d:
                if d==uday:
                    dum.append(o)
            ans[uday]= dum
        return ans

    @staticmethod
    def read_astac(file):
        with open(file) as inf:
            head=inf.readline().split()
            observations = []
            for line in inf:
                ob=dict()
                parts = line.split()
                for h,a in zip(head,parts):
                    ob[h]=a
                observations.append(ob)
        return observations
