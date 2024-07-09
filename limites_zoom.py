# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 13:07:24 2024

@author: manon
"""

import json
#from collections.abc import Iterable

# import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

from plots import Map
from readers_3D import MesoNH, get_mesonh, get_time_index, index_to_lonlat
from matplotlib.animation import ArtistAnimation

DX250_ZOOM = (
    ((600, 860), (497, 1397), 5, 0, 0),
    ((950, 1250), (966, 1400), 6, 30, 1),
    ((1000, 1370), (1200, 1500), 7, 0, 2),
    ((1150, 1450), (1260, 1530), 7, 15, 3),
    ((1440, 1790), (1530, 1730), 8, 15, 4),
    ((1470, 1940), (1650, 1930), 8, 45, 5),
)

def sum_clouds(rct, rit, rgt, rst):
    """Add up different thickness of the condensed states of water."""
    return rct + rit + rgt + rst

def norm_wind(ut, vt):
    return np.sqrt(ut**2 + vt**2)

resol_dx = 250
index_echeance = 0
zoom = DX250_ZOOM

mesonh = get_mesonh(resol_dx)
i_lim, j_lim, hour, minute, file_index = zoom[index_echeance]

mesonh.get_data(file_index)

#var = (mesonh.get_var("PABST", level=range(0, 92))[:, j_lim[0]:j_lim[1], i_lim[0]:i_lim[1]])/100
var = (mesonh.get_var("RCT", "RIT", "RGT", "RST", func=sum_clouds, level=range(0, 92))[:, j_lim[0]:j_lim[1], i_lim[0]:i_lim[1]])*100
#var = (mesonh.get_var("RVT", level=range(0, 92))[:, j_lim[0]:j_lim[1], i_lim[0]:i_lim[1]])*100
#var = (mesonh.get_var("UT", "VT", func=norm_wind, level=range(0, 92))[:, j_lim[0]:j_lim[1], i_lim[0]:i_lim[1]])*3.6
#var = (mesonh.get_var("WT", level=range(0, 92))[:, j_lim[0]:j_lim[1], i_lim[0]:i_lim[1]])
#var = (mesonh.get_var("THT", level=range(0, 92))[:, j_lim[0]:j_lim[1], i_lim[0]:i_lim[1]])-273.15
#var = (mesonh.get_var("TKET", level=range(0, 92))[:, j_lim[0]:j_lim[1], i_lim[0]:i_lim[1]])



print(var.min(), var.max())


