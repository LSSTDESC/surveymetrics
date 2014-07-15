import numpy

class LSST:

    filterinfo=dict()
    _filterNames=numpy.array(['u','g','r','i','z','y'])
    _filterWave = numpy.array([375.,476.,621.,754.,870.,980.])
    for name, wave in zip(_filterNames, _filterWave):
        filterinfo[name]=dict()
        filterinfo[name]['wave']=wave
