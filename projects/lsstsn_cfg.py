# Here is an example of a very very simple MAF configuration driver script
# to run:
# runDriver.py most_simple_cfg.py

# This script uses the LSST pex_config.  This is executed as a python script, but only things that 
#  start with 'root.' are passed on to the driver script.

# Import MAF helper functions 
from lsst.sims.maf.driver.mafConfig import makeBinnerConfig, makeMetricConfig, makeDict

# Set the output directory
root.outputDir = './snmetric_lsst_out'
# Set the database to use (the example db included in the git repo)
#root.dbAddress = {'dbAddress':'sqlite:///opsimblitz2_1060_sqlite.db'}
root.dbAddress = {'dbAddress':'sqlite:////Users/akim/Downloads/opsimblitz2_1060_sqlite.db'}
# Name of this run (filename base)
root.opsimName = 'snmetric_lsst'

nside=16

import healpy
solidangle = healpy.nside2pixarea(nside, degrees=True)
import numpy
solidangle = numpy.asscalar(solidangle)

# Configure a metric to run. Compute the mean on the final delivered seeing.  Once the mean seeing has been computed everywhere on the sky, compute the RMS as a summary statistic.

metric = makeMetricConfig('snmetric_lsst',summaryStats={'RmsMetric':{}})

# Configure a binner.  Use the Healpixbinner to compute the metric at points in the sky.  Set the constraint as an empty string so all data is returned.
binner = makeBinnerConfig('HealpixBinner',kwargs={"nside":nside},  metricDict=makeDict(metric),
                          constraints=[''])

root.binners = makeDict(binner)

