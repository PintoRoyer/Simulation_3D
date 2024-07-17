from datetime import datetime, timedelta

import numpy as np
from netCDF4 import Dataset

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
import cartopy.crs as ccrs
import cartopy.feature as cfeature


plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 13,
    # "figure.facecolor":  (0, 0, 0, 0)
})


def get_level_index_from_pressure(pressure, press_target: float = 850.):
    nb_levels = pressure.shape[0]
    p_mean = np.zeros(nb_levels)
    for lvl in range(nb_levels):
        p_mean[lvl] = np.mean(pressure[lvl])

    delta = np.abs(press_target - p_mean)
    index = np.array(np.where(delta == delta.min())).T
    return index[0][0]


def calc_vorticity_divergence(filename, i_min=0, i_max=-1, j_min=0, j_max=-1, barbs_mesh=100, pressure=850):
    # Open Meso-NH file
    data = Dataset(filename)
    time = datetime.strptime("2022-08-18 00:00:00", "%Y-%m-%d %H:%M:%S") + timedelta(
        seconds=int(data.variables["time"][0])
    )
    level = get_level_index_from_pressure(data.variables["PABST"][0, :, :, :] * 1e-2, pressure)
    print(level)

    vorticity = np.zeros((2050, 2050))
    divergence = np.zeros((2050, 2050))
    u_wind = data.variables["UT"][0, level, :, :]
    v_wind = data.variables["VT"][0, level, :, :]

    for j in range(j_min, j_max + 1):
        for i in range(i_min, i_max + 1):
            vorticity[j, i] = (
                (v_wind[j, i + 1] - v_wind[j, i]) / 250
                - (u_wind[j + 1, i] - u_wind[j, i]) / 250
            )
            divergence[j, i] = (
                (u_wind[j, i + 1] - u_wind[j, i]) / 250
                + (v_wind[j + 1, i] - v_wind[j, i]) / 250
            )

    fig = plt.figure(figsize=(16, 10), layout="compressed")
    fig.suptitle(f"Simulation Méso-NH du {time} (DX = 250 m) au niveau {level} (~{data.variables['level'][level]:.1f} m, ~{pressure} hPa)")

    # Plotting vorticity
    axes = fig.add_subplot(1, 2, 1, projection=ccrs.PlateCarree())
    axes.add_feature(cfeature.COASTLINE, color="gray")
    axes.add_feature(cfeature.BORDERS, color="gray")
    glines = axes.gridlines(draw_labels=True, linewidth=0.4, alpha=0.4)
    glines.top_labels = glines.right_labels = False

    cf = axes.contourf(
        data.variables["longitude"][j_min: j_max, i_min: i_max],
        data.variables["latitude"][j_min: j_max, i_min: i_max],
        vorticity[j_min: j_max, i_min: i_max],
        cmap="seismic",
        norm=TwoSlopeNorm(0),
        levels=np.linspace(-0.1, 0.2, 100),
        extend="both"
    )
    plt.colorbar(cf, label="Vorticité potentielle du vent horizontal (1/s)", format="%.3f", location="bottom")
    axes.barbs(
        data.variables["longitude"][j_min: j_max: barbs_mesh, i_min: i_max: barbs_mesh],
        data.variables["latitude"][j_min: j_max: barbs_mesh, i_min: i_max: barbs_mesh],
        u_wind[j_min: j_max: barbs_mesh, i_min: i_max:barbs_mesh],
        v_wind[j_min: j_max: barbs_mesh, i_min: i_max: barbs_mesh],
        color="black",
        length=5,
        alpha=0.5
    )
    axes.set_title("Vorticité potentielle du vent horizontal")

    # Plotting divergence
    axes = fig.add_subplot(1, 2, 2, projection=ccrs.PlateCarree())
    axes.add_feature(cfeature.COASTLINE, color="gray")
    axes.add_feature(cfeature.BORDERS, color="gray")
    glines = axes.gridlines(draw_labels=True, linewidth=0.4, alpha=0.4)
    glines.top_labels = glines.right_labels = False

    cf = axes.contourf(
        data.variables["longitude"][j_min: j_max, i_min: i_max],
        data.variables["latitude"][j_min: j_max, i_min: i_max],
        divergence[j_min: j_max, i_min: i_max],
        cmap="seismic",
        norm=TwoSlopeNorm(0),
        levels=np.linspace(-0.05, 0.05, 100),
        extend="both"
    )
    plt.colorbar(cf, label="Divervence du vent horizontal (1/s)", format="%.3f", location="bottom")
    axes.barbs(
        data.variables["longitude"][j_min: j_max: barbs_mesh, i_min: i_max: barbs_mesh],
        data.variables["latitude"][j_min: j_max: barbs_mesh, i_min: i_max: barbs_mesh],
        u_wind[j_min: j_max: barbs_mesh, i_min: i_max:barbs_mesh],
        v_wind[j_min: j_max: barbs_mesh, i_min: i_max: barbs_mesh],
        color="black",
        length=5,
        alpha=0.5
    )
    axes.set_title("Divervence du vent horizontal")

    plt.savefig(f"DX250_{str(time.hour).zfill(2)}{str(time.minute).zfill(2)}Z_rot_{pressure}hPa.png")


# i_min, i_max, j_min, j_max, barbs_mesh, file_id
DX250_ZOOM = (
    (600 , 860 , 497 , 1397, 30, "004"),  # 0500Z
    (950 , 1250, 966 , 1400, 20, "010"),  # 0630Z
    (1000, 1370, 1200, 1500, 20, "012"),  # 0700Z
    (1150, 1450, 1260, 1530, 20, "013"),  # 0715Z
    (1440, 1790, 1530, 1730, 20, "017"),  # 0815Z
    (1470, 1940, 1650, 1930, 30, "019"),  # 0845Z
)

for press in (850, 500, 300):
    calc_vorticity_divergence(f"CORSE.1.SEG01.{DX250_ZOOM[3][-1]}.vars.nc", *DX250_ZOOM[3][:-1], pressure=press)
