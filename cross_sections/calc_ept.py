from datetime import datetime, timedelta

import numpy as np
from netCDF4 import Dataset

import metpy.calc as mpcalc
from metpy.units import units
from collections import namedtuple


Variable = namedtuple("Variable", ["values", "unit"])


def calc(func, *variables):
    result = []
    for j in range(variables[0].values.shape[0]):
        line = []
        for i in range(variables[0].values.shape[1]):
            line.append(func(*[var.values[j, i] * var.unit for var in variables]).magnitude)
        result.append(line)
    return np.array(result)


def level_at_850hpa(pressure):
    nb_levels = pressure.shape[0]
    p_mean = np.zeros(nb_levels)
    for lvl in range(nb_levels):
        p_mean[lvl] = np.mean(pressure[lvl])

    delta = np.abs(850 - p_mean)
    index = np.array(np.where(delta == delta.min())).T
    return index[0][0]



def calc_ept(filename):
    # Open Meso-NH file
    data = Dataset(f"/mesonh/panf/ADASTRA/CORSE/DX250/{filename}")
    time = datetime.strptime("2022-08-18 00:00:00", "%Y-%m-%d %H:%M:%S") + timedelta(
        seconds=int(data.variables["time"][0])
    )

    # Defining the level and zoom
    lvl_index = level_at_850hpa(data.variables["PABST"][0, :, :, :] * 1e-2)

    # Extract variables
    potential_temp = data.variables["THT"][0, lvl_index, :, :]
    pressure = data.variables["PABST"][0, lvl_index, :, :]
    mixing_ratio = data.variables["RVT"][0, lvl_index, :, :]

    # Compute the equivalent potential temperature
    temperature = calc(
        mpcalc.temperature_from_potential_temperature,
        Variable(pressure, units("Pa")),
        Variable(potential_temp, units("K"))
    )

    relative_humidity = calc(
        mpcalc.relative_humidity_from_mixing_ratio,
        Variable(pressure, units("Pa")),
        Variable(temperature, units("K")),
        Variable(mixing_ratio, units("kg/kg"))
    )

    dewpoint = calc(
        mpcalc.dewpoint_from_relative_humidity,
        Variable(temperature, units("K")),
        Variable(relative_humidity, units("%"))
    )

    equivalent_potential_temperature = calc(
        mpcalc.equivalent_potential_temperature,
        Variable(pressure, units("Pa")),
        Variable(temperature, units("K")),
        Variable(dewpoint, units("°C"))
    )

    with open(
        f"DX250_{str(time.hour).zfill(2) + str(time.minute).zfill(2)}Z_EPT_LVL28.npy", "wb"
    ) as file:
        np.save(file, equivalent_potential_temperature)


# for filename in (
#     "CORSE.1.SEG01.004.vars.nc",
#     "CORSE.1.SEG01.010.vars.nc",
#     "CORSE.1.SEG01.012.vars.nc",
#     "CORSE.1.SEG01.013.vars.nc",
#     "CORSE.1.SEG01.017.vars.nc",
#     "CORSE.1.SEG01.019.vars.nc"
# ):
#     calc_ept(filename)

data = Dataset("CORSE.1.SEG01.013.vars.nc")
idx = level_at_850hpa(data.variables["PABST"][0, :, :, :] * 1e-2, 300)

print(idx)
print(data.variables["level"][idx])
print(np.mean(data.variables["PABST"][0, idx] * 1e-2))