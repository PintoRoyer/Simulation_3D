# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 11:39:36 2024

@author: manon
"""

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from matplotlib.animation import ArtistAnimation

data = xr.open_mfdataset("../Données/DX250/CORSE.1.SEG01.OUT.*.nc")

fig = plt.figure()
axes = plt.axes(projection=ccrs.PlateCarree())
axes.coastlines(linewidth=1, color="black", alpha=0.4)
axes.add_feature(cfeature.BORDERS, color="black", linewidth=1, alpha=0.4)

glines = axes.gridlines(draw_labels=True, linewidth=0.5)
glines.top_labels = glines.right_labels = False

levels = np.linspace(data.ACPRR.min(), data.ACPRR.max(), 25)

frame = []

for time in range(data.ACPRR.shape[0]):
    n = int(time * 50 / data.ACPRR.shape[0])
    print("[" + n * "#" + (50 - n) * " " + "]", end="\r")
    contourf = axes.contourf(
        data.longitude[0], data.latitude[:, 0], data.ACPRR[time], cmap="Blues", levels=levels
    )
    date, hour = str(data.time[time].values).split("T")
    title = axes.text(
        0.5,
        1.05,
        f"Simulation Méso-NH le {date} à {hour[: 5]} TU\n(dx = 250 m)",
        ha="center",
        transform=axes.transAxes,
    )
    frame.append([contourf, title])


plt.colorbar(contourf, label="Taux de précipitations accumulées (mm)", format="%.0e", fraction=0.03)

animation = ArtistAnimation(fig, frame, interval=250)

animation.save("essai.gif")

plt.show()

 #plt.savefig(f"pressure_{time}_{zoom}_{resol_dx}m.png")

 # # Clouds
 # axes = my_map.init_axes(
 #     fig_kw={"figsize": (8, 5), "layout": "compressed"},
 #     feature_kw={"linewidth": 1, "alpha": 0.5, "color": "white"},
 # )[1]
 # axes.set_extent([lon[0], lon[1], lat[0], lat[1]])

 # var = mesonh.get_var("THCW", "THRW", "THIC", "THSN", "THGR", func=sum_clouds)
 # contourf = my_map.plot_contourf(
 #     var,
 #     cmap=LinearSegmentedColormap.from_list("cmap2", ["black", "white", "blue", "red"]),
 #     levels=np.linspace(lim["clouds"][0], lim["clouds"][1], 100),
 # )
 # cbar = plt.colorbar(contourf, label="Épaisseur nuageuse (mm)")
 # cbar.set_ticks(np.linspace(lim["clouds"][0], lim["clouds"][1], 8))
 # plt.savefig(f"clouds_{time}_{zoom}_{resol_dx}m.png")

 # # Wind speed
 # axes = my_map.init_axes(
 #     fig_kw={"figsize": (8, 5), "layout": "compressed"}, feature_kw={"color": "black"}
 # )[1]
 # axes.set_extent([lon[0], lon[1], lat[0], lat[1]])
 # var = mesonh.get_var("WIND10", func=lambda x: x * 3.6)
 # contourf = my_map.plot_contourf(
 #     var,
 #     cmap=LinearSegmentedColormap.from_list(
 #         "cmap2",
 #         [
 #             "white",
 #             (240 / 255, 248 / 255, 255 / 255),
 #             "darkcyan",
 #             "yellow",
 #             "orange",
 #             "red",
 #             "purple",
 #             "black",
 #         ],
 #     ),
 #     levels=np.linspace(lim["wind"][0], lim["wind"][1], 100),
 # )
 # cbar = plt.colorbar(contourf, label="Vitesse du vent horizontal (km/h)")
 # cbar.set_ticks(np.linspace(lim["wind"][0], lim["wind"][1], 8))

 # # Wind direction
 # if i_lim == (0, -1):
 #     kwargs = {"scale": 40}

 # else:
 #     mesh = 25
 #     if resol_dx == 500:
 #         mesh = 12
 #     elif resol_dx == 1000:
 #         mesh = 6
 #     kwargs = {
 #         "x_mesh": mesh,
 #         "y_mesh": mesh,
 #         "width": 0.005,
 #         "scale": 20,
 #         "scale_units": "xy",
 #         "units": "xy",
 #     }

 # wind_u, wind_v = mesonh.get_var("UM10", "VM10", "WIND10", func=norm_wind)
 # my_map.plot_quiver(wind_u, wind_v, **kwargs)
 # plt.savefig(f"wind_{time}_{zoom}_{resol_dx}m.png")

def sum_clouds(thcw, thrw, thic, thsn, thgr):
    """Add up different thickness of the condensed states of water."""
    return thcw + thrw + thic + thsn + thgr


def norm_wind(um10, vm10, wind10):
    """Normalize the wind components."""
    return um10 / wind10, vm10 / wind10