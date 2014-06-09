#!/usr/bin/env python
import numpy
#from pyslalib import slalib
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

def desanal(file):
    omega=3.
    T0=30.
    tau=5.
    dt0=4.
    a=1.
    
    bands = ['g','r','i','z']

    ans=0.
    
    # get data
    obs = Metric.read_astac(file)

    # pick out only those in SN fields
    Metric.prune_by_fieldname(obs,desfields)

    # do one field at a time
    for f in desfields:
        obs_f=Metric.get_by_fieldname(obs, f)
        if f[-1]=='3':
            zmax=desetime['deep']['zmax']
        else:
            zmax=desetime['shallow']['zmax']

        # do one band at a time
        for b in bands:
            if f[-1]=='3':
                snr0=numpy.sqrt(desetime['deep'][b])
            else:
                snr0=numpy.sqrt(desetime['shallow'][b])
            
            metric = Metric.OneFieldBandMetric(omega,zmax,snr0,T0, tau, dt0,a)
          
            obs_0=Metric.get_by_band(obs_f,b)
 
            # observations into nights
            obs_0 = Metric.group_by_night(obs_0,longitude)

            snrs=[]
            dates=[]
            # group each nights observations into a single effective observation
            for date in obs_0.keys():
                os = obs_0[date]
                efftime=0
                for o in os:
                    if (o['bad'] == 'False') and (float(o['exptime'])>10):
                        efftime=efftime+float(o['t_eff'])*float(o['exptime'])

                if efftime !=0:
                    snrs.append(numpy.sqrt(efftime))
                    dates.append(date)
                                
            snrs=numpy.array(snrs)
            dates=numpy.array(dates)

            dsort=dates.argsort()
            snrs=snrs[dsort]
            dates=dates[dsort]

            dates = numpy.append(dates,dates[-1]+tau)
            
            #calculate metric for observation in bands
            ans=ans+ metric.metric(snrs,dates)
    return ans

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser( description = "Remove targets close to database" )
    parser.add_argument( "filename"       , help = "input filename with target coadd_objects_id"   )
    args = parser.parse_args()
    pdict=vars(args)

    print desanal(pdict['filename'])

#./metric_des  ../data/schedule15/obssim-19/obs.txt
