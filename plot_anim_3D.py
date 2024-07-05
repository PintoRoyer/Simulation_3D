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
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

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


DX250_ZOOM = (
    ((600, 860), (497, 1397), 5, 0, 0),
    ((950, 1250), (966, 1400), 6, 30, 1),
    ((1000, 1370), (1200, 1500), 7, 0, 2),
    ((1150, 1450), (1260, 1530), 7, 15, 3),
    ((1440, 1790), (1530, 1730), 8, 15, 4),
    ((1470, 1940), (1650, 1930), 8, 45, 5),
)

DX1000_ZOOM = (
    ((148, 215), (123, 348), 5, 0, 0),
    ((230, 312), (242, 362), 6, 30, 1),
    ((247, 342), (298, 381), 7, 0, 2),
    ((275, 363), (314, 389), 7, 15, 3),
    ((320, 446), (377, 453), 8, 15, 4),
    ((383, 483), (414, 472), 8, 45, 5),
)

#PRESSURE
def plot_pressure(mesonh, my_map, zoom_name, resol_dx, level, *, anim = False):
    # # Limits for colorbars
    with open(f"limits_3D/lim3D_{zoom_name}_{resol_dx}m.json", "r", encoding="utf-8") as file:
        lim = json.loads(file.read())
    lim_min = lim["pressure"][0]
    lim_max = lim["pressure"][1]
    
    
    
    var = mesonh.get_var("PABST", level=level)/100
    var_lim = var[j_lim[0]:j_lim[1], i_lim[0]:i_lim[1]]
    contourf = my_map.plot_contourf(
            var, cmap="turbo", levels=np.linspace(var_lim.min(), var_lim.max(), 100)
        )
    
    level_value = int(mesonh.level[level])
    title = f"Simulation Méso-NH du 2022-08-18 à {time} (DX = {resol_dx} m)\nPression au niveau {level} (~{level_value} m) "
    title = my_map.set_title(title)
    
    variable = "pressure"
    label = "Pression (hPa)"
    
    if anim == False:
        cbar = plt.colorbar(contourf, label=label)
        cbar.set_ticks(np.linspace(var_lim.min(), var_lim.max(), 8))
        plt.savefig(f"pressure_lvl{level}_{time}_dx{resol_dx}.png")

    return (contourf, title), time, lim_min, lim_max, variable, label


#CLOUDS
def sum_clouds(rct, rit, rgt, rst):
    """Add up different thickness of the condensed states of water."""
    return rct + rit + rgt + rst


def plot_clouds(mesonh, my_map, zoom_name, resol_dx, level, *, anim = False):
    # # Limits for colorbars
    with open(f"limits_3D/lim3D_{resol_dx}m.json", "r", encoding="utf-8") as file:
        lim = json.loads(file.read())
    lim_min = 0
    lim_max = 19.55
    
    var = mesonh.get_var("RCT", "RIT", "RGT", "RST", func=sum_clouds, level=level)*1000
    contourf = my_map.plot_contourf(
            var, cmap=LinearSegmentedColormap.from_list("cmap2", ["black", "white", "blue", "red"]), levels=np.linspace(lim_min, lim_max, 100)
        )
    
    level_value = int(mesonh.level[level])
    title = f"Simulation Méso-NH du 2022-08-18 à {time} (DX = {resol_dx} m)\nSomme des rapports de mélange des états condensés\nde l'eau au niveau {level} (~{level_value} m) "
    title = my_map.set_title(title)
    
    variable = "clouds"
    label = "Somme des rapports de mélange (g/kg)"
    
    if anim == False:
        cbar = plt.colorbar(contourf, label=label)
        cbar.set_ticks(np.linspace(lim_min, lim_max, 8))
        plt.savefig(f"clouds_lvl{level}_{time}_dx{resol_dx}.png")

    return (contourf, title), time, lim_min, lim_max, variable, label


#WATER VAPOR
def plot_vapor(mesonh, my_map, zoom_name, resol_dx, level, *, anim = False):
    # # Limits for colorbars
    with open(f"limits_3D/lim3D_{resol_dx}m.json", "r", encoding="utf-8") as file:
        lim = json.loads(file.read())
    lim_min = lim["vapor"][0]
    lim_max = lim["vapor"][1]
    
    var = mesonh.get_var("RVT", level=level)*1000
    contourf = my_map.plot_contourf(
            var, cmap=LinearSegmentedColormap.from_list("cmap2", ["black", "white", "blue", "red"]), levels=np.linspace(lim_min, lim_max, 100)
        )
    
    level_value = int(mesonh.level[level])
    title = f"Simulation Méso-NH du 2022-08-18 à {time} (DX = {resol_dx} m)\nRapport de mélange de vapeur d'eau au niveau {level} (~{level_value} m) "
    title = my_map.set_title(title)
    
    variable = "vapor"
    label = "Somme des rapports de mélange (g/kg)"
    
    if anim == False:
        cbar = plt.colorbar(contourf, label=label)
        cbar.set_ticks(np.linspace(lim_min, lim_max, 8))
        plt.savefig(f"vapor_lvl{level}_{time}_dx{resol_dx}.png")

    return (contourf, title), time, lim_min, lim_max, variable, label


#HORIZONTAL WIND
def norm_wind(ut, vt):
    return np.sqrt(ut**2 + vt**2)

def norm_quiver(ut, vt):
    norm = norm_wind(ut, vt)
    return ut / norm, vt / norm

def plot_horiz_wind(mesonh, my_map, zoom_name, resol_dx, level, *, anim = False):
    # # Limits for colorbars
    
    with open(f"limits_3D/lim3D_{resol_dx}m.json", "r", encoding="utf-8") as file:
        lim = json.loads(file.read())
    lim_min = lim["horiz_wind"][0]
    lim_max = lim["horiz_wind"][1]
    
    var = mesonh.get_var("UT", "VT", func=norm_wind, level=level)*3.6
    contourf = my_map.plot_contourf(
            var, cmap=LinearSegmentedColormap.from_list("cmap2",
                [
                    "white",
                    (240 / 255, 248 / 255, 255 / 255),
                    "darkcyan",
                    "yellow",
                    "orange",
                    "red",
                    "purple",
                    "black",
                ]), levels=np.linspace(lim_min, lim_max, 100)
        )
    if resol_dx == 1000 :    
        mesh=20
    elif resol_dx == 500 :
        mesh=40
    else : 
        mesh=80
    kwargs = {
        "x_mesh": mesh,
        "y_mesh": mesh,
        "width": 0.015,
        "scale": 7,
        "scale_units": "xy",
        "units": "xy",
    }
    UT_norm, VT_norm = mesonh.get_var("UT", "VT", func=norm_quiver, level=level)
    quiver = my_map.plot_quiver(UT_norm, VT_norm, **kwargs)
    
    level_value = int(mesonh.level[level])
    title = f"Simulation Méso-NH du 2022-08-18 à {time} (DX = {resol_dx} m)\nNorme et direction du vent horizontal au niveau {level} (~{level_value} m) "
    title = my_map.set_title(title)
    
    variable = "horiz_wind"
    label = "Norme du vent horizontal (km/h)"
    
    if anim == False:
        cbar = plt.colorbar(contourf, label=label)
        cbar.set_ticks(np.linspace(lim_min, lim_max, 8))
        plt.savefig(f"{variable}_lvl{level}_{time}_dx{resol_dx}.png")

    return (contourf, title, quiver), time, lim_min, lim_max, variable, label


#VERTICAL WIND
def plot_vert_wind(mesonh, my_map, zoom_name, resol_dx, level, *, anim = False):
    # # # Limits for colorbars
    # with open(f"limits_3D/lim3D_{resol_dx}m.json", "r", encoding="utf-8") as file:
    #     lim = json.loads(file.read())
    # lim_min = lim["vert_wind"][0]
    # lim_max = lim["vert_wind"][1]
    lim_min = -15
    lim_max = 30
    
    var = mesonh.get_var("WT", level=level)
    contourf = my_map.plot_contourf(
            var, cmap="seismic", norm=TwoSlopeNorm(0), extend='both', levels=np.linspace(lim_min, lim_max, 100)
        )
    
    level_value = int(mesonh.level[level])
    title = f"Simulation Méso-NH du 2022-08-18 à {time} (DX = {resol_dx} m)\nNorme du vent vertical au niveau {level} (~{level_value} m) "
    title = my_map.set_title(title)
    
    variable = "vert_wind"
    label = "Norme du vent vertical (m/s)"
    
    if anim == False:
        cbar = plt.colorbar(contourf, label=label)
        cbar.set_ticks(np.linspace(lim_min, lim_max, 8))
        plt.savefig(f"{variable}_lvl{level}_{time}_dx{resol_dx}.png")

    return (contourf, title), time, lim_min, lim_max, variable, label


def plot_anim(mesonh, my_map, zoom_name, resol_dx, func):
    # print(mesonh.get_limits("PABST")) #/100
    # print(mesonh.get_limits("RCT", "RIT", "RGT", "RST", func=sum_clouds)) #*1000
    # print(mesonh.get_limits("RVT")) #*1000
    # print(mesonh.get_limits("UT", "VT", func=norm_wind)) #*3.6
    # print(mesonh.get_limits("WT"))
    
    frame = []
    for level in range(22, 84): #len(mesonh.level)
        frame_content, time, lim_min, lim_max, variable, label = func(mesonh, my_map, zoom_name, resol_dx, level, anim = True)
        
        frame.append(frame_content)
        
    cbar = plt.colorbar(frame_content[0], label=label)
    cbar.set_ticks(np.linspace(116, 911, 8))
    
    animation = ArtistAnimation(fig, frame, interval=250)
    animation.save(f"{variable}_anim_{time}_dx{resol_dx}.gif")
    
    
resol_dx = 1000
level = 86
index_echeance = 0
zoom = DX1000_ZOOM
zoom_name = "DX1000_ZOOM"

mesonh = get_mesonh(resol_dx)
i_lim, j_lim, hour, minute, file_index = zoom[index_echeance]

mesonh.get_data(file_index)
time = f"{str(hour).zfill(2)}h{str(minute).zfill(2)}"

# Creating Map instance
my_map = Map(mesonh.longitude, mesonh.latitude)

# Information on zoom
lon = [0, 0]  # bornes min max lon
lat = [0, 0]  # bornes min max lat
lon[0], lat[0] = index_to_lonlat(mesonh, i_lim[0], j_lim[0])
lon[1], lat[1] = index_to_lonlat(mesonh, i_lim[1], j_lim[1])

fig, axes, _ = my_map.init_axes(fig_kw={"figsize": (8, 5), "layout": "compressed"}, feature_kw={"color" : "gray"})
axes.set_extent([lon[0], lon[1], lat[0], lat[1]])


#plot_pressure(mesonh, my_map, zoom_name, resol_dx, level, )
plot_anim(mesonh, my_map, zoom_name, resol_dx, plot_pressure)

#plot_clouds(mesonh, my_map, zoom_name, resol_dx, level, )
#plot_anim(mesonh, my_map, zoom_name, resol_dx, plot_clouds)

#plot_vapor(mesonh, my_map, zoom_name, resol_dx, level, )
#plot_anim(mesonh, my_map, zoom_name, resol_dx, plot_vapor)

#plot_horiz_wind(mesonh, my_map, zoom_name, resol_dx, level, )
#plot_anim(mesonh, my_map, zoom_name, resol_dx, plot_horiz_wind)

#plot_vert_wind(mesonh, my_map, zoom_name,  resol_dx, level, )
#plot_anim(mesonh, my_map, zoom_name, resol_dx, plot_vert_wind)