#!/usr/bin/env python3

from datetime import datetime
from netCDF4 import Dataset
import pandas as pd
from scipy.interpolate import LinearNDInterpolator


ghf = pd.read_csv("Antarctic_GHF.xyz", sep='\s+', header=None,
                  names=["x", "y", "Heat Flux mW/m^2"])

ghf_x = ghf.x.values
ghf_y = ghf.y.values
ghf_val = ghf["Heat Flux mW/m^2"].values / 1000.0

crane = Dataset("Crane.nc", "r+")
x_crane = crane.variables["xCell"][:]
y_crane = crane.variables["yCell"][:]

interp = LinearNDInterpolator(list(zip(ghf_x, ghf_y)), ghf_val)
crane.variables["basalHeatFlux"][0, :] = interp(x_crane, y_crane)

new_history_str = datetime.now().strftime("%a %b %d %H:%M:%S %Y") + \
    " Interpolated geothermal flux from Martos et al. (2017) " + \
    "using scipy.interpolate.LinearNDInterpolator"
newhist = '\n'.join([new_history_str, getattr(crane, 'history')])

setattr(crane, 'history', newhist)

crane.close()
