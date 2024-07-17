from datetime import datetime, timedelta

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import cartopy.crs as ccrs
import cartopy.feature as cfeature


plt.rcParams.update({"text.usetex": True, "font.family": "serif", "font.size": 13})

DX250_ZOOM = (
    ((600, 860), (497, 1397), 5, 0),
    ((950, 1250), (966, 1400), 6, 30),
    ((1000, 1370), (1200, 1500), 7, 0),
    ((1150, 1450), (1260, 1530), 7, 15),
    ((1440, 1790), (1530, 1730), 8, 15),
    ((1470, 1940), (1650, 1930), 8, 45),
)


i_min, i_max, j_min, j_max = 1150, 1450, 1260, 1530

data = Dataset("CORSE.1.SEG01.013.vars.nc")
clouds = (
    data.variables["RIT"][0, 60, :, :]
    + data.variables["RST"][0, 60, :, :]
    + data.variables["RGT"][0, 60, :, :]
) * 1e3
lon = data.variables["longitude"][:]
lat = data.variables["latitude"][:]


fig = plt.figure(figsize=(8, 5), layout="compressed")
time = datetime.strptime("2022-08-18 00:00:00", "%Y-%m-%d %H:%M:%S") + timedelta(
    seconds=int(data.variables["time"][0])
)
# fig.suptitle(f"Simulation Méso-NH du {time} (DX = 250 m)\nSomme des rapports de mélanges des états condensés à {data.variables['level'][60]:.1f} m")

axes = fig.add_subplot(projection=ccrs.PlateCarree())


axes.add_feature(cfeature.COASTLINE, color="gray")
axes.add_feature(cfeature.BORDERS, color="gray")
glines = axes.gridlines(draw_labels=True, linewidth=0.4, alpha=0.5)
glines.top_labels = glines.right_labels = False


cf = axes.contourf(
    lon,
    lat,
    clouds,
    cmap=LinearSegmentedColormap.from_list("", ["black", "white", "blue", "red"]),
    levels=np.linspace(clouds.min(), clouds.max(), 100)
)
cb = fig.colorbar(
    cf,
    label="Somme des rapports de mélanges des états condensés (g/kg)"
)
cb.set_ticks(np.linspace(clouds.min(), clouds.max(), 8))

lon_min = lon[j_min, i_min]
lon_max = lon[j_max, i_max]
lat_min = lat[j_min, i_min]
lat_max = lat[j_max, i_max]

axes.plot(
    [lon_min, lon_min, lon_max, lon_max, lon_min],
    [lat_min, lat_max, lat_max, lat_min, lat_min],
    color="red",
    linewidth=2
)

plt.show()