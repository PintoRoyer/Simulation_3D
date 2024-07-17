import matplotlib.pyplot as plt
from netCDF4 import Dataset

LVL_INDEX = 20
i_min, i_max = 1150, 1450
j_min, j_max = 1260, 1530

data = Dataset("CORSE.1.SEG01.013.vars.nc")
press = data.variables["PABST"][0, LVL_INDEX, i_min: i_max, j_min: j_max] * 1e-2
print(data.variables["level"][LVL_INDEX])


fig = plt.figure()
axes = fig.add_subplot(1, 1, 1, projection="3d")
axes.plot_wireframe(
    data.variables["longitude"][i_min: i_max, j_min: j_max],
    data.variables["latitude"][i_min: i_max, j_min: j_max],
    press
)
plt.show()
