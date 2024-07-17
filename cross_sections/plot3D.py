import numpy as np
from netCDF4 import Dataset
import pygmt


data = Dataset("CORSE.1.SEG01.004.vars.nc")
LON = data.variables["longitude"][:, :]
LAT = data.variables["latitude"][:, :]
LVL = data.variables["level"][:]

clouds = (
    + data.variables["RCT"][0, :, :, :]
    + data.variables["RIT"][0, :, :, :]
    + data.variables["RST"][0, :, :, :]
    + data.variables["RGT"][0, :, :, :]
) * 1e3

clouds = np.where(clouds > 0.1, 1, 0)
for i in range(clouds.shape[0]):
    clouds[i] = LVL[i] * clouds[i]
  

grid = pygmt.datasets.load_earth_relief(
    resolution="10m",
    region=[LON.min(), LON.max(), LAT.min(), LAT.max()]
)

fig = pygmt.Figure()
fig.grdview(
    grid=grid,
    perspective=[130, 30],
    frame=["xa", "yaf", "WSnE"],
    projection="M15c",
    zsize="1.5c",
    surftype="s",
    cmap="geo",
)



fig.show()
    
