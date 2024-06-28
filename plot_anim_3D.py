"""
Zooms
=====

Description
-----------
Plot animations of clouds, pressure and wind at different levels at one time step for the different resolution of the simulation.
"""

import json
#from collections.abc import Iterable

# import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

from plots import Map
from readers_3D import MesoNH, get_mesonh, get_time_index, index_to_lonlat
from matplotlib.animation import ArtistAnimation

# Matplotlib configuration to have LaTeX style
# plt.rcParams.update({"text.usetex": True, "font.family": "serif", "font.size": 15})

NO_ZOOM = (
    ((0, -1), (0, -1), 5, 0, 0),
    ((0, -1), (0, -1), 6, 30, 1),
    ((0, -1), (0, -1), 7, 0, 2),
    ((0, -1), (0, -1), 7, 15, 3),
    ((0, -1), (0, -1), 8, 15, 4),
    ((0, -1), (0, -1), 8, 45, 5),
)


#PRESSURE
def plot_pressure(mesonh, my_map, resol_dx, level, *, anim = False):

    var = mesonh.get_var("PABST", level=level)/100
    contourf = my_map.plot_contourf(
            var, cmap="turbo",extend="both", levels=np.linspace(lim_min, lim_max, 100)
        )
    
    level_value = int(mesonh.level[level])
    title = f"Simulation Méso-NH du 2022-08-18 à {time} (DX = {resol_dx} m)\nPression au niveau {level} (~{level_value} m) "
    title = my_map.set_title(title)
    
    
    if anim == False:
        cbar = plt.colorbar(contourf, label="Pression (hPa)")
        cbar.set_ticks(np.round(np.linspace(lim_min, lim_max, 8)))
        plt.savefig(f"pressure_lvl{level}_{time}_dx{resol_dx}.png")

    return contourf, title, time, lim_min, lim_max

#CLOUDS
def sum_clouds(rct, rit, rgt, rst):
    """Add up different thickness of the condensed states of water."""
    return rct + rit + rgt + rst

def plot_clouds(mesonh, my_map, resol_dx, level, *, anim = False):
    var = mesonh.get_var("RCT", "RIT", "RGT", "RST", func=sum_clouds, level=level)
    contourf = my_map.plot_contourf(
            var, cmap=LinearSegmentedColormap.from_list("cmap2", ["black", "white", "blue", "red"]), extend="both", levels=np.linspace(var.min(), var.max(), 100)
        )
    
    level_value = int(mesonh.level[level])
    title = f"Simulation Méso-NH du 2022-08-18 à {time} (DX = {resol_dx} m)\nPression au niveau {level} (~{level_value} m) "
    title = my_map.set_title(title)
    
    
    if anim == False:
        cbar = plt.colorbar(contourf, label="Pression (hPa)")
        cbar.set_ticks(np.round(np.linspace(var.min(), var.max(), 8)))
        plt.savefig(f"clouds_lvl{level}_{time}_dx{resol_dx}.png")

    return contourf, title, time, lim_min, lim_max


def plot_anim(mesonh, my_map, resol_dx, func):
    frame = []
    
    for level in range(0, len(mesonh.level)):
        #print(mesonh.get_limits("PABST"))
        contourf, title, time, lim_min, lim_max = func(mesonh, my_map, resol_dx, level, anim = True)
        
        frame.append([contourf, title])
        
    cbar = plt.colorbar(contourf, label="Pression (hPa)")
    cbar.set_ticks(np.round(np.linspace(lim_min, lim_max, 8)))
    
    animation = ArtistAnimation(fig, frame, interval=250)
    animation.save(f"pressure_anim_{time}_dx{resol_dx}.gif")
    
    
resol_dx = 1000
level = 80
index_echeance = 2


mesonh = get_mesonh(resol_dx)
i_lim, j_lim, hour, minute, file_index = NO_ZOOM[index_echeance]

mesonh.get_data(file_index)
time = f"{str(hour).zfill(2)}h{str(minute).zfill(2)}"

# Creating Map instance
my_map = Map(mesonh.longitude, mesonh.latitude)

# Information on zoom
lon = [0, 0]  # bornes min max lon
lat = [0, 0]  # bornes min max lat
lon[0], lat[0] = index_to_lonlat(mesonh, i_lim[0], j_lim[0])
lon[1], lat[1] = index_to_lonlat(mesonh, i_lim[1], j_lim[1])

# # Limits for colorbars
with open(f"limits_3D/lim3D_{resol_dx}m.json", "r", encoding="utf-8") as file:
    lim = json.loads(file.read())
lim_min = lim["pressure"][0]
lim_max = lim["pressure"][1]

fig, axes, _ = my_map.init_axes(fig_kw={"figsize": (8, 5), "layout": "compressed"})
axes.set_extent([lon[0], lon[1], lat[0], lat[1]])


#plot_pressure(mesonh, my_map, resol_dx, level, )
#plot_anim(mesonh, my_map, resol_dx, plot_pressure)

#plot_clouds(mesonh, my_map, resol_dx, level, )
plot_anim(mesonh, my_map, resol_dx, plot_clouds)
