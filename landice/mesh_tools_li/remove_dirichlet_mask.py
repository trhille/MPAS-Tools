#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 16:05:22 2020
This script resets dirichletVelocityMask to zero everywhere.
@author: trevorhillebrand
"""

import netCDF4
from optparse import OptionParser


print("== Gathering information.  (Invoke with --help for more details. All arguments are optional)\n")
parser = OptionParser()
parser.description = __doc__
parser.add_option("-f", "--file", dest="inputFile", help="name of file to be modified.", default="landice_grid.nc", metavar="FILENAME")
for option in parser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
options, args = parser.parse_args()


print("  Input file: {}".format(options.inputFile))

f = netCDF4.Dataset(options.inputFile, 'r+')

# Remove dirichlet boundary conditions
dirichletVelocityMask = f.variables['dirichletVelocityMask'][:] * 0

f.variables['dirichletVelocityMask'][:] = dirichletVelocityMask

f.close()

