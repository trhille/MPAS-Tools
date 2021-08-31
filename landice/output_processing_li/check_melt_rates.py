#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Aug 31, 2021
This checks mean melt rates under floating ice.
@author: trevorhillebrand
"""
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser

print("** Gathering information.  (Invoke with --help for more details. All arguments are optional)")
parser = OptionParser(description=__doc__)
parser.add_option("-f", dest="filenames", help="filename for plotting", metavar="FILENAME")
options, args = parser.parse_args()

filenames = options.filenames.split(',') #create list of filenames to loop over

floatValue = 4 # value in cellMask denoting floating ice
dynamicValue = 2 # value in cellMask denoting dynamic ice.
fig, ax = plt.subplots(1,1)
for filename in filenames:
    melt = Dataset(filename, 'r')
    daysSinceStart = melt.variables["daysSinceStart"][:]
    cellMask = melt.variables["cellMask"][:]
    areaCell = melt.variables["areaCell"][:]
    dynamicMask = (cellMask & dynamicValue) // dynamicValue
    floatMask = (cellMask & floatValue) // floatValue
    basalMassBal = melt.variables["basalMassBal"][:]
    dynamicFloatMask = floatMask & dynamicMask
    meanFloatMelt = 0. * daysSinceStart
    floatMelt = (basalMassBal * dynamicFloatMask) < 0.
    for t in np.arange(0, len(daysSinceStart)):
        meanFloatMelt[t] = np.sum(areaCell[floatMelt[t,:]] *
                             basalMassBal[t, floatMelt[t,:]]) / np.sum(areaCell[floatMelt[t,:]])

    ax.plot(daysSinceStart / 365.,  meanFloatMelt * 3.154e7 / 910., label=filename)
    melt.close()

ax.set_xlabel('Year')
ax.set_ylabel('Mean melt rate under\n dynamic floating ice (m/yr)')
ax.grid('on')
ax.legend(fontsize='small')
ax.set_ylim(top=0)
plt.show()
