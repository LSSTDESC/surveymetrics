#!/usr/bin/env python
from distutils.core import setup

setup(name='surveymetrics',
      version='0.0',
      description='Metrics to Assess Astronomical Surveys',
      author='Alex Kim',
      author_email='agkim@lbl.gov',
      packages=['surveymetrics', 'projects'],
      scripts=['projects/snmetric_des']
     )
