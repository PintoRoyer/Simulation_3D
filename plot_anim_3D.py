"""
Zooms
=====

Description
-----------
Plot animations of clouds, pressure and wind at different levels at one time step for the different resolution of the simulation.
"""

import json
from collections.abc import Iterable

# import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
# from matplotlib.colors import LinearSegmentedColormap

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
def plot_pressure(resol_dx, level, args: Iterable):
    """
    Plot and save at format .png pressure maps at a certain level, timestep and resolution

    Parameters
    ----------
    resol_dx : int
        The resolution of the simulation.
    level : int
        The level of altitude at which the variables are observed 
    args : iterable
        The arguments to be given to ``plot_pressure`` it should be like:

            args = (
                ((i_min, i_max), (j_min, j_max), hour, minute),
            )
    """
    
    mesonh = get_mesonh(resol_dx)
    

    for i_lim, j_lim, hour, minute, file_index in args:
        mesonh.get_data(file_index)
        time = f"{str(hour).zfill(2)}h{str(minute).zfill(2)}"

        plt.close("all")

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
    
        fig, axes, _ = my_map.init_axes(fig_kw={"figsize": (8, 5), "layout": "compressed"})
        axes.set_extent([lon[0], lon[1], lat[0], lat[1]])
      
    
        var = mesonh.get_var("PABST", level=level)/100
        contourf = my_map.plot_contourf(
                var, cmap="turbo",extend="both", levels=np.linspace(lim["pressure"][0], lim["pressure"][1], 100)
            )
        
        cbar = plt.colorbar(contourf, label="Pression (hPa)")
        cbar.set_ticks(np.round(np.linspace(lim["pressure"][0], lim["pressure"][1], 8)))
        
        print(var.min())
        
        level_value = int(mesonh.level[level])
        plt.title(f"Simulation Méso-NH du 2022-08-18 à {time}\nPression au niveau {level} (~{level_value} m) ")
        plt.savefig(f"pressue_lvl{level}_{time}_dx{resol_dx}.png")


def plot_pressure_anim(resol_dx, args: Iterable):
    """
    Plot and save at format .gif pressure maps at a certain level, timestep and resolution

    Parameters
    ----------
    resol_dx : int
        The resolution of the simulation.
    args : iterable
        The arguments to be given to ``plot_pressure`` it should be like:

            args = (
                ((i_min, i_max), (j_min, j_max), hour, minute),
            )
    """
    mesonh = get_mesonh(resol_dx)

    for i_lim, j_lim, hour, minute, file_index in args:
        mesonh.get_data(file_index)
        time = f"{str(hour).zfill(2)}h{str(minute).zfill(2)}"
        
        #print(mesonh.get_limits("PABST"))
        
        plt.close("all")
    
        # Creating Map instance
        my_map = Map(mesonh.longitude, mesonh.latitude)
    
        # Information on zoom
        lon = [0, 0]  # bornes min max lon
        lat = [0, 0]  # bornes min max lat
        lon[0], lat[0] = index_to_lonlat(mesonh, i_lim[0], j_lim[0])
        lon[1], lat[1] = index_to_lonlat(mesonh, i_lim[1], j_lim[1])
    
        # Limits for colorbars
        with open(f"limits_3D/lim3D_{resol_dx}m.json", "r", encoding="utf-8") as file:
            lim = json.loads(file.read())
        
        fig, axes, _ = my_map.init_axes(fig_kw={"figsize": (8, 5), "layout": "compressed"})
        axes.set_extent([lon[0], lon[1], lat[0], lat[1]])
    
        frame = []    
    
        for level in range(0, len(mesonh.level)):
            var = mesonh.get_var("PABST", level=level)/100
            contourf = my_map.plot_contourf(
                var, cmap="turbo", extend="both", levels=np.linspace(lim["pressure"][0], lim["pressure"][1], 100)
            )
            level_value = int(mesonh.level[level])
            title = f"Simulation Méso-NH du 2022-08-18 à {time}\nPression au niveau {level} (~{level_value} m) "
            title = my_map.set_title(title)
            
            frame.append([contourf, title])
            
        cbar = plt.colorbar(contourf, label="Pression (hPa)")
        cbar.set_ticks(np.round(np.linspace(lim["pressure"][0], lim["pressure"][1], 8)))
        
        animation = ArtistAnimation(fig, frame, interval=250)
        
        animation.save(f"pressure_anim_{time}_dx{resol_dx}.gif")




#plot_pressure(1000, 91, NO_ZOOM[0: 1])
plot_pressure_anim(1000, NO_ZOOM[0:1])

