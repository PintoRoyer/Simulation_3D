"""Compute and save the equivalent potential temperature from Meso-NH file."""

from collections import namedtuple
from collections.abc import Callable
from datetime import datetime, timedelta

import metpy.calc as mpcalc
import numpy as np
from metpy.units import units
from netCDF4 import Dataset

Variable = namedtuple("Variable", ["values", "unit"])


def calc(func: Callable, *variables):
    """
    Apply a metpy function over given 2D variables.

    Parameters
    ----------
    func : Callable
        The function to apply.
    *variables
        The Variables instances to be given to ``func``.

    Returns
    -------
    np.array
        The array of the result.
    """
    result = []
    for j in range(variables[0].values.shape[0]):
        line = []
        for i in range(variables[0].values.shape[1]):
            line.append(func(*[var.values[j, i] * var.unit for var in variables]).magnitude)
        result.append(line)
    return np.array(result)


def level_at_850hpa(pressure: np.array):
    """
    Find the index for the level closest to 850hPa. As the pressure variates, ``level_at_850hpa``
    averages it for each level.

    Parameters
    ----------
    pressure : np.array
        The 3D array of pressure.

    Returns
    -------
    int
        The index of the level closest to 850hPa.
    """
    nb_levels = pressure.shape[0]
    p_mean = np.zeros(nb_levels)
    for lvl in range(nb_levels):
        p_mean[lvl] = np.mean(pressure[lvl])

    delta = np.abs(850 - p_mean)
    index = np.array(np.where(delta == delta.min())).T
    return index[0][0]


def calc_ept(filename: str):
    """
    Compute equivalent potential temperature from a given file. The result will be stored into a
    *.npy file, please see NumPy documentation for more information about *.npy file handling.

    .. warning::
        This function is designed to work with the data on NUWA, you should change path to use it
        outside.

    Parameters
    ----------
    filename : str
        The name of the Meso-NH file to open.
    """
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

    # Compute the temperature
    temperature = calc(
        mpcalc.temperature_from_potential_temperature,
        Variable(pressure, units("Pa")),
        Variable(potential_temp, units("K")),
    )

    # Compute the relative humidity
    relative_humidity = calc(
        mpcalc.relative_humidity_from_mixing_ratio,
        Variable(pressure, units("Pa")),
        Variable(temperature, units("K")),
        Variable(mixing_ratio, units("kg/kg")),
    )

    # Compute the dewpoint
    dewpoint = calc(
        mpcalc.dewpoint_from_relative_humidity,
        Variable(temperature, units("K")),
        Variable(relative_humidity, units("%")),
    )

    # Compute the equivalent potential temperature
    equivalent_potential_temperature = calc(
        mpcalc.equivalent_potential_temperature,
        Variable(pressure, units("Pa")),
        Variable(temperature, units("K")),
        Variable(dewpoint, units("°C")),
    )

    # Save the data
    with open(
        f"DX250_{str(time.hour).zfill(2) + str(time.minute).zfill(2)}Z_EPT_LVL28.npy", "wb"
    ) as file:
        np.save(file, equivalent_potential_temperature)


if __name__ == "__main__":
    # Be careful: this script takes several hours and a lot of RAM
    for mesonh_file in (
        "CORSE.1.SEG01.004.vars.nc",
        "CORSE.1.SEG01.010.vars.nc",
        "CORSE.1.SEG01.012.vars.nc",
        "CORSE.1.SEG01.013.vars.nc",
        "CORSE.1.SEG01.017.vars.nc",
        "CORSE.1.SEG01.019.vars.nc",
    ):
        calc_ept(mesonh_file)
