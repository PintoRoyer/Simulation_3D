import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset


data = Dataset("CORSE.1.SEG01.013.vars.nc")
ni = data.variables["ni"][:]
nj = data.variables["nj"][:]
clouds = (
    data.variables["RCT"][0, :, :, :]
    + data.variables["RIT"][0, :, :, :]
    + data.variables["RST"][0, :, :, :]
    + data.variables["RGT"][0, :, :, :]
) * 1e3


fig = plt.figure()
axes = fig.add_subplot(1, 1, 1, projection="3d")
for i in range(0, clouds.shape[0] // 2):
    lvl = np.where(clouds[i], i, np.nan)
    axes.scatter(ni, nj, lvl, color="gray")

plt.show()
