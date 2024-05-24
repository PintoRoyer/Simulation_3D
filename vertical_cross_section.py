# -*- coding: utf-8 -*-
"""
Created on Thu May 16 13:45:56 2024

@author: manon
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.colors import LinearSegmentedColormap

import xarray as xr
# import wrf


files = xr.open_dataset("../Donnees/CORSE.1.SEG01.013.vars.nc")

# RRT = files.data_vars["PABST"][0, 0: 10, 700:900, 930] * 1e-2

x_start = 1322
x_end   = 1322
y_start = 500 #‗1300
y_end   = 1800 #1450
z_start = 0
z_end   = 85


clouds = (
    files.data_vars["RRT"][0, z_start: z_end, x_start, y_start: y_end] +
    files.data_vars["RIT"][0, z_start: z_end, x_start, y_start: y_end] +
    files.data_vars["RGT"][0, z_start: z_end, x_start, y_start: y_end] +
    files.data_vars["RCT"][0, z_start: z_end, x_start, y_start: y_end] +
    files.data_vars["RST"][0, z_start: z_end, x_start, y_start: y_end]
) * 1e3

PABST = files.data_vars["PABST"][0, z_start: z_end, x_start, y_start: y_end] * 1e-2

wind_w = np.array(files.data_vars["WT"][0, z_start: z_end, x_start, y_start: y_end] * 3.6)
wind_v = np.array(files.data_vars["VT"][0, z_start: z_end, x_start, y_start: y_end] * 3.6)
wind_norm_vw = np.sqrt(wind_w**2 + wind_v**2)

ni = files.coords["ni"][x_start: x_end]
nj = files.coords["nj"][y_start: y_end]
ALT = files.coords['level'][z_start: z_end] * 1e-3
lat = files.coords['latitude'][y_start: y_end, x_start].values


cmap = LinearSegmentedColormap.from_list("cmap2", [
    (0, (1, 1, 1, 0)),
    (0.1, (0, 0, 1, 0.5)),
    (0.3, (0.3, 0, 0.7, 0.8)),
    (0.6, (0.5, 0, 0.5, 0.8)),
    (1, (1, 0, 0, 1))
])


cmap1 = LinearSegmentedColormap.from_list("cmap2", [
    (0, (1, 1, 1, 0)),
    (0.2, (1, 1, 1, 0.1)),
    (0.25, (1, 1, 1, 0.4)),
    (0.55, "blue"),
    (1, "red")
])


plt.figure(figsize=(18,5))

#cf = plt.contourf(nj, ALT, PABST, levels=np.linspace(PABST.min(), PABST.max(), 100), cmap="turbo")
c = plt.contour(nj, ALT, clouds, levels=(0.1, ), colors=("white", ))
cf = plt.contourf(nj, ALT, wind_norm_vw, levels=np.linspace(wind_norm_vw.min(), wind_norm_vw.max(), 100), cmap="viridis")

index_per_quiver = 30
plt.quiver(
    nj[::index_per_quiver],
    ALT[::5],
    (wind_v/wind_norm_vw)[::5, ::index_per_quiver],
    (wind_w/wind_norm_vw)[::5, ::index_per_quiver]
)

#plt.colorbar(cf, label="Somme des rapports de mélange\ndes états condensés de l'eau (g/kg)") 
plt.colorbar(cf, label="Vitesse verticale du vent (km/h)")
plt.xticks(nj[::60], np.round(lat[::60], 2))

plt.grid(axis="y")

# RRT_section = wrf.vertcross( field3d=RRT.values, vert=ALT.values, start_point=wrf.CoordPair(674, 884), end_point=wrf.CoordPair(876, 884))

plt.savefig("cross_section_test.png")
#plt.savefig("cross_section_wind_w.png")