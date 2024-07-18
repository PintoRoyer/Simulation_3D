"""Plot the equivalent potential temperature."""
from datetime import datetime, timedelta

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
import cartopy.crs as ccrs
import cartopy.feature as cfeature

plt.rcParams.update({"text.usetex": True, "font.family": "serif", "font.size": 13})


# Openning data
mnh_data = Dataset("CORSE.1.SEG01.013.vars.nc")
data = Dataset("DX250_0715Z_EPT_LVL28.nc")
ept = data.variables["equivalent_potential_temperature"][0, 0, :, :]

# Creating fig
fig = plt.figure(figsize=(8, 5), layout="compressed")
axes = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

# Set title
time = datetime.strptime("2022-08-18 00:00:00", "%Y-%m-%d %H:%M:%S") + timedelta(
    seconds=int(data.variables["time"][0])
)
# fig.suptitle(
#     f"Simulation Méso-NH du {time} (DX = 250 m) à 850 hPa\nTempérature potentielle équivalente\n"
#     f"Ascendances et subsidences\nContour nuageux"
# )
# fig.suptitle(
#     f"Simulation Méso-NH du {time} (DX = 250 m) à 850 hPa\nVent vertical et horizontal"
# )

# Configuring axes
axes.add_feature(cfeature.COASTLINE, color="gray")
axes.add_feature(cfeature.BORDERS, color="gray")
glines = axes.gridlines(draw_labels=True, linewidth=0.4, alpha=0.4)
glines.top_labels = glines.right_labels = False

# Plot equivalent potential temperature
contourf = axes.contourf(
    data.variables["longitude"][:, :],
    data.variables["latitude"][:, :],
    ept,
    cmap="turbo",
    levels=np.linspace(ept.min(), ept.max(), 100)
)
cbar = plt.colorbar(contourf, label="Température potentielle équivalente (K)")
cbar.set_ticks(np.linspace(ept.min(), ept.max(), 10))

# Vertical winds
# contourf = axes.contourf(
#     mnh_data.variables["longitude"][:, :],
#     mnh_data.variables["latitude"][:, :],
#     mnh_data.variables["WT"][0, 28, :, :],
#     cmap=LinearSegmentedColormap.from_list("", (
#         (0.0, (0, 0, 1, 1)),
#         (0.5, (1, 1, 1, 0)),
#         (1.0, (1, 0, 0, 1))
#     )),
#     levels=np.linspace(-2, 10, 100),
#     extend="both",
#     norm=TwoSlopeNorm(0),
#     # alpha=1
# )
# cbar = plt.colorbar(contourf, label="Vitesse du vent vertical (m/s)")
# cbar.set_ticks(np.linspace(-2, 10, 10))

# Horizontal wind
# axes.barbs(
#     mnh_data.variables["longitude"][::80, ::80],
#     mnh_data.variables["latitude"][::80, ::80],
#     mnh_data.variables["UT"][0, 28, ::80, ::80],
#     mnh_data.variables["VT"][0, 28, ::80, ::80],
#     color="black",
#     length=5,
#     alpha=0.5
# )

# Clouds and precipitations
# axes.contour(
#     mnh_data.variables["longitude"][:, :],
#     mnh_data.variables["latitude"][:, :],
#     (
#         mnh_data.variables["RCT"][0, 28, :, :]
#         + mnh_data.variables["RST"][0, 28, :, :]
#         + mnh_data.variables["RIT"][0, 28, :, :]
#         + mnh_data.variables["RGT"][0, 28, :, :]
#     ) * 1e3,
#     levels=[0.1],
#     colors=["black"],
#     linewidths=1,
#     linestyles="dashed",
#     alpha=1
# )
# axes.contour(
#     mnh_data.variables["longitude"][:, :],
#     mnh_data.variables["latitude"][:, :],
#     mnh_data.variables["RRT"][0, 28, :, :] * 1e3,
#     levels=[0.1],
#     colors=["green"],
#     linewidths=2,
#     linestyles="dashed",
#     alpha=1
# )

# Display fig
plt.show()
