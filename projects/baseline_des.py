#!/usr/bin/env python
import numpy
#from pyslalib import slalib
import surveymetrics.snmetric

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

header = "start_time  start_mjd   start_clock sun_zd  moon_zd moon_brightness tactician   duration    hexname tiling  filter  exptime slew_time   slew_angle  ra  decl    airmass moon_angle  skymag  seeing  clouds  bad t_eff   pred_t_eff  pred_seeing\n"
dum1 = "2014-08-19T04:49:42Z    "
firstmjd = 56888.20118332838
dum2 = "   1408423782.239572   162.62244648490318  113.44810579407248  0.04270449897197033 SupernovaDeadmanTactician   204.0   "
#field exptime band
dum3= "   175.0   29  0.0 9.5 -43.998 1.2152378703662599  87.19906119354471   21.926952290815812  0.9426964355542654  0   False   "
# teff
dum4 = " 0.6108206352904961  0.8913507052003073"

deltamjd = 5.
nepoch = 25
bands =['g','r','i','z']

f = open('baseline.txt', 'w')
f.write(header)

for i in xrange(nepoch):
    for j in desfields:
        if j[-1]=='3':
            key = 'deep'
        else:
            key = 'shallow'
        for k in bands:
            f.write(dum1+str(firstmjd+i*deltamjd)+dum2+j+" "+str(desetime[key][k])+" "+k+dum3+str(desetime[key][k])+dum4+"\n")
f.close()

