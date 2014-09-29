# script to install RSAT dependencies
# questions to fcaprist@stanford.edu

import sys
import os


import commands


compiler = 'gfortran' # or 'gfortran or intelem'

status, mainPath = commands.getstatusoutput('pwd')


print '\n INSTALLING CYTHON STATS ROUTINE-KDE AND QUADTREE IMPLEMENTATION'
STATSPATH = mainPath+'/STATS/kernelQuadCpp/'

os.chdir(STATSPATH)

os.system('python setup.py install --install-lib=../')

print 'INSTALLING STATS F2PY ROUTINE'

os.chdir(mainPath+'/STATS')
os.system('f2py --fcompiler='+compiler+' --quiet -m meancov -c meanCov2D.f90')
os.system('f2py --fcompiler='+compiler+' --quiet -m lonlatChecks -c lonlatChecks.f90')

print '\n INSTALLING SAFETY METRIC F2PY ROUTINES'

SAFETYPATH = mainPath+'/SafetyMetrics'
os.chdir(SAFETYPATH)
os.system('f2py --fcompiler='+compiler+' --quiet -m safetyMetrics -c safetyMetricFortran.f90')

print 'INSTALLING TRAJECTORY AND DEBRIS PROPAGATION F2PY ROUTINES'
TRAJPATH = mainPath+'/TrajectoryPropagation/SourceCode'
os.chdir(TRAJPATH)
if compiler=='intelem':
    os.system('ifort -c -fPIC ode.f90')
else:
    os.system('gfortran -c -fPIC ode.f90')

os.system('f2py --fcompiler='+compiler+' --quiet -m orbitProp -c moduleConstants.f90 modulePlanetEarth.f90 denEstimator.f90 linearInterpolation.f90 speedOfSound.f90 ode.o trajProp.f90')

os.system('f2py --fcompiler='+compiler+' --quiet -m debrisPropagation -c moduleConstants.f90 modulePlanetEarth.f90 denEstimator.f90 linearInterpolation.f90 speedOfSound.f90 ode.o debrisProp.f90')
os.system('f2py --fcompiler='+compiler+' --quiet -m debrisPropagationAllTime -c moduleConstants.f90 modulePlanetEarth.f90 denEstimator.f90 linearInterpolation.f90 speedOfSound.f90 ode.o debrisPropAllTime.f90')
os.system('f2py --fcompiler='+compiler+' --quiet -m malFunc -c moduleConstants.f90 modulePlanetEarth.f90 denEstimator.f90 linearInterpolation.f90 speedOfSound.f90 ode.o malfuncTurn.f90')

os.system('mv orbitProp.so ../orbitProp.so')
os.system('mv debrisPropagation.so ../debrisPropagation.so')
os.system('mv debrisPropagationAllTime.so ../debrisPropagationAllTime.so')
os.system('mv malFunc.so ../malFunc.so')

