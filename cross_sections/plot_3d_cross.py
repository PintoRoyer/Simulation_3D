"""Some functions to plot cross sections with different variables from Meso-NH file."""

from collections import namedtuple
from collections.abc import Callable
from datetime import datetime, timedelta

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import metpy.calc as mpcalc
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from metpy.units import units
from netCDF4 import Dataset

from misc_functions import oblique_proj

plt.rcParams.update(
    {"text.usetex": True, "font.family": "serif", "font.size": 13, "figure.facecolor": (0, 0, 0, 0)}
)


DX250_ZOOM = (
    ((600, 860), (497, 1397), 5, 0),
    ((950, 1250), (966, 1400), 6, 30),
    ((1000, 1370), (1200, 1500), 7, 0),
    ((1150, 1450), (1260, 1530), 7, 15),
    ((1440, 1790), (1530, 1730), 8, 15),
    ((1470, 1940), (1650, 1930), 8, 45),
)

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


def add_contourf(data: Dataset, var: np.array, axes: plt.Axes, label: str, **contourf_kw):
    """
    Add a contourf to the given axes.

    Parameters
    ----------
    data : Dataset
        The Meso-NH Dataset.
    var : np.array
        The variable to plot on the axes as a contourf.
    axes : plt.Axes
        The axes instance to plot ``var`` on.
    label : str
        The label of the colorbar.
    **contourf_kw
        Some keyword arguments to be passed to ``plt.contourf``.
    """
    if "levels" not in contourf_kw:
        contourf_kw["levels"] = np.linspace(var.min(), var.max(), 100)

    contourf = axes.contourf(
        data.variables["longitude"][:, :], data.variables["latitude"][:, :], var, **contourf_kw
    )
    cbar = plt.colorbar(contourf, label=label, location="bottom")
    cbar.set_ticks(np.linspace(contourf_kw["levels"].min(), contourf_kw["levels"].max(), 10))


def add_cross_contourf(
    data: Dataset,
    var: np.array,
    axes: plt.Axes,
    label: str,
    i_min: int,
    j_min: int,
    i_max: int,
    j_max: int,
    **contourf_kw,
):
    """
    Add a cross section contourf on the given axes.

    Parameters
    ----------
    data : Dataset
        The Meso-NH Dataset.
    var : np.array
        The variable to plot on the axes as a contourf.
    axes : plt.Axes
        The axes instance to plot ``var`` on.
    label : str
        The label of the colorbar.
    i_min : int
        The minimum index on x-axis.
    i_max : int
        The maximum index on x-axis.
    j_min : int
        The minimum index on y-axis.
    j_max : int
        The minimum index on y-axis.
    **contourf_kw
        Some keyword arguments to be passed to ``plt.contourf``.
    """
    level = data.variables["level"][:]
    var_cross = oblique_proj(
        var, data.variables["ni"][:], data.variables["nj"][:], level, i_min, j_min, i_max, j_max
    )[1]

    if "levels" not in contourf_kw:
        contourf_kw["levels"] = np.linspace(var_cross.min(), var_cross.max(), 100)
    elif isinstance(contourf_kw["levels"], int):
        contourf_kw["levels"] = np.linspace(var_cross.min(), var_cross.max(), contourf_kw["levels"])

    contourf = axes.contourf(
        np.array(range(var_cross.shape[1])),
        np.array(range(var_cross.shape[0])),
        var_cross,
        **contourf_kw,
    )
    cbar = plt.colorbar(contourf, ax=axes, label=label, location="bottom")
    cbar.set_ticks(np.linspace(contourf_kw["levels"].min(), contourf_kw["levels"].max(), 6))

    axes.set_xticks([0, var_cross.shape[1] - 1], ["A", "B"])
    ticks_level = range(0, var_cross.shape[0], var_cross.shape[0] // 10)
    axes.set_yticks(ticks_level, [f"{int(i)} m" for i in level[ticks_level]])
    axes.set_ylabel("Altitude (m)")


def rh_cross(data: Dataset, axes: plt.Axes, i_min: int, j_min: int, i_max: int, j_max: int):
    """
    Relative humidity cross section

    Parameters
    ----------
    data : Dataset
        The Meso-NH dataset.
    axes : plt.Axes
        The axes to draw on.
    i_min : int
        The minimum index on x-axis.
    i_max : int
        The maximum index on x-axis.
    j_min : int
        The minimum index on y-axis.
    j_max : int
        The minimum index on y-axis.
    """
    pressure = oblique_proj(
        data.variables["PABST"][0, :, :, :],
        data.variables["ni"][:],
        data.variables["nj"][:],
        data.variables["level"][:],
        i_min,
        j_min,
        i_max,
        j_max,
    )[1]
    potential_temperature = oblique_proj(
        data.variables["THT"][0, :, :, :],
        data.variables["ni"][:],
        data.variables["nj"][:],
        data.variables["level"][:],
        i_min,
        j_min,
        i_max,
        j_max,
    )[1]
    mixing_ratio = oblique_proj(
        data.variables["RVT"][0, :, :, :],
        data.variables["ni"][:],
        data.variables["nj"][:],
        data.variables["level"][:],
        i_min,
        j_min,
        i_max,
        j_max,
    )[1]

    temperature = calc(
        mpcalc.temperature_from_potential_temperature,
        Variable(pressure, units("Pa")),
        Variable(potential_temperature, units("K")),
    )
    print("Temperature: ok")

    relative_humidity = calc(
        mpcalc.relative_humidity_from_mixing_ratio,
        Variable(pressure, units("Pa")),
        Variable(temperature, units("K")),
        Variable(mixing_ratio, units("kg/kg")),
    )
    print("Relative humidity: ok")

    # Relative humidity contourf
    contourf = axes.contourf(
        np.array(range(relative_humidity.shape[1])),
        np.array(range(relative_humidity.shape[0])),
        relative_humidity,
        cmap=LinearSegmentedColormap.from_list("", ["black", "white", "blue", "red"]),
        levels=np.linspace(0, 100, 100),
    )
    cbar = plt.colorbar(contourf, ax=axes, label=r"Humidité relative (\%)", location="bottom")
    cbar.set_ticks(np.linspace(0, 100, 6))

    # Clouds contour
    var = (
        data.variables["RCT"][0, :, :, :]
        + data.variables["RIT"][0, :, :, :]
        + data.variables["RST"][0, :, :, :]
        + data.variables["RGT"][0, :, :, :]
    ) * 1e3
    var_cross = oblique_proj(
        var,
        data.variables["ni"][:],
        data.variables["nj"][:],
        data.variables["level"][:],
        i_min,
        j_min,
        i_max,
        j_max,
    )[1]
    axes.contour(
        np.array(range(var_cross.shape[1])),
        np.array(range(var_cross.shape[0])),
        var_cross,
        levels=[0.1],
        colors=["black"],
    )

    # Wind barbs
    # wt_cross = oblique_proj(
    #     data.variables["WT"][0, :, :, :],
    #     data.variables["ni"][:],
    #     data.variables["nj"][:],
    #     data.variables["level"][:],
    #     i_min,
    #     j_min,
    #     i_max,
    #     j_max,
    # )[1]

    # ut = data.variables["UT"][0, :, :, :]
    # ut = ut - np.mean(Dataset("CORSE.1.SEG01.019.vars.nc").variables["UT"][0, 28, 0: 1500, 0: 1250])
    # ut_cross = oblique_proj(
    #     ut,
    #     data.variables["ni"][:],
    #     data.variables["nj"][:],
    #     data.variables["level"][:],
    #     i_min,
    #     j_min,
    #     i_max,
    #     j_max,
    # )[1]
    # axes.barbs(
    #     np.array(range(0, ut_cross.shape[1], 40)),
    #     np.array(range(0, ut_cross.shape[0], 5)),
    #     ut_cross[::5, ::40],
    #     wt_cross[::5, ::40],
    #     color="gray",
    #     length=5,
    # )

    axes.set_xticks([0, relative_humidity.shape[1] - 1], ["A", "B"])
    ticks_level = range(0, relative_humidity.shape[0], relative_humidity.shape[0] // 10)
    axes.set_yticks(ticks_level, [f"{int(i)} m" for i in data.variables["level"][ticks_level]])
    axes.set_ylabel("Altitude (m)")


def rvt_map(data: Dataset, axes: plt.Axes):
    """
    Plot vapour mixing ratio map.

    Parameters
    ----------
    data : Dataset
        The Meso-NH dataset.
    axes : plt.Axes
        The axes to draw on.
    """

    var = data.variables["RVT"][0, :, :, :] * 1e3
    add_contourf(
        data,
        var[42],
        axes,
        "Rapport de mélange en vapeur d'eau à 700 hPa (g/kg)",
        cmap=LinearSegmentedColormap.from_list("", ["black", "white", "blue", "red"]),
    )


def rvt_cross(data: Dataset, axes: plt.Axes, i_min: int, j_min: int, i_max: int, j_max: int):
    """
    Plot vapour mixing ratio cross section.

    Parameters
    ----------
    data : Dataset
        The Meso-NH dataset.
    axes : plt.Axes
        The axes to draw on.
    i_min : int
        The minimum index on x-axis.
    i_max : int
        The maximum index on x-axis.
    j_min : int
        The minimum index on y-axis.
    j_max : int
        The minimum index on y-axis.
    """
    var = data.variables["RVT"][0, :, :, :] * 1e3
    add_cross_contourf(
        data,
        var,
        axes,
        "Rapport de mélange en vapeur d'eau (g/kg)",
        i_min,
        j_min,
        i_max,
        j_max,
        cmap=LinearSegmentedColormap.from_list("", ["black", "white", "blue", "red"]),
    )

    var = data.variables["PABST"][0, :, :, :] * 1e-2
    var_cross = oblique_proj(
        var,
        data.variables["ni"][:],
        data.variables["nj"][:],
        data.variables["level"][:],
        i_min,
        j_min,
        i_max,
        j_max,
    )[1]
    c = axes.contour(
        np.array(range(var_cross.shape[1])),
        np.array(range(var_cross.shape[0])),
        var_cross,
        levels=[700],
        colors=["gray"],
    )
    plt.clabel(c, levels=[700])


def pressure_map(data: Dataset, axes: plt.Axes):
    """Plot
    pressure map.

    Parameters
    ----------
    data : Dataset
        The Meso-NH dataset.
    axes : plt.Axes
        The axes to draw on.
    """
    var = data.variables["PABST"][0, :, :, :] * 1e-2
    add_contourf(
        data,
        var[0],
        axes,
        "Pression à la surface (hPa)",
        cmap="turbo",
    )


def pressure_cross(data: Dataset, axes: plt.Axes, i_min: int, j_min: int, i_max: int, j_max: int):
    """
    Plot pressure cross section.

    Parameters
    ----------
    data : Dataset
        The Meso-NH dataset.
    axes : plt.Axes
        The axes to draw on.
    i_min : int
        The minimum index on x-axis.
    i_max : int
        The maximum index on x-axis.
    j_min : int
        The minimum index on y-axis.
    j_max : int
        The minimum index on y-axis.
    """
    var = data.variables["PABST"][0, :, :, :] * 1e-2
    add_cross_contourf(
        data,
        var,
        axes,
        "Pression (hPa)",
        i_min,
        j_min,
        i_max,
        j_max,
        cmap="turbo",
        extend="both",
    )


def tket_map(data: Dataset, axes: plt.Axes):
    """P
    lot turbulent kinetic energy map.

    Parameters
    ----------
    data : Dataset
        The Meso-NH dataset.
    axes : plt.Axes
        The axes to draw on.
    """
    var = data.variables["TKET"][0, :, 1260:1530, 1150:1450]
    add_contourf(
        data,
        var[77],
        axes,
        "Énergie synétique turbulente (m²/s²)",
        cmap="inferno",
        extend="max",
        levels=np.linspace(0, 120, 100),
    )


def clouds_map(data: Dataset, axes: plt.Axes):
    """Plo
    t clouds map.

    Parameters
    ----------
    data : Dataset
        The Meso-NH dataset.
    axes : plt.Axes
        The axes to draw on.
    """
    var = (
        data.variables["RIT"][0, :, :, :]
        + data.variables["RST"][0, :, :, :]
        + data.variables["RGT"][0, :, :, :]
    )[:, 1260:1530, 1150:1450] * 1e3
    add_contourf(
        data,
        var[77],
        axes,
        ("Somme des rapports de mélanges des états condensés"),
        cmap=LinearSegmentedColormap.from_list("", ["black", "white", "blue", "red"]),
    )


def clouds_cross(data: Dataset, axes: plt.Axes, i_min: int, j_min: int, i_max: int, j_max: int):
    """
    Plot clouds cross section.

    Parameters
    ----------
    data : Dataset
        The Meso-NH dataset.
    axes : plt.Axes
        The axes to draw on.
    i_min : int
        The minimum index on x-axis.
    i_max : int
        The maximum index on x-axis.
    j_min : int
        The minimum index on y-axis.
    j_max : int
        The minimum index on y-axis.
    """
    var = (
        data.variables["RCT"][0, :, :, :]
        + data.variables["RIT"][0, :, :, :]
        + data.variables["RST"][0, :, :, :]
        + data.variables["RGT"][0, :, :, :]
    ) * 1e3
    add_cross_contourf(
        data,
        var,
        axes,
        "Somme des rapports de mélanges des états condensés (g/kg)",
        i_min,
        j_min,
        i_max,
        j_max,
        cmap=LinearSegmentedColormap.from_list("", ["black", "white", "blue", "red"]),
    )


def wind_w_map(data: Dataset, axes: plt.Axes):
    """Plo
    t vertical wind speed map.

    Parameters
    ----------
    data : Dataset
        The Meso-NH dataset.
    axes : plt.Axes
        The axes to draw on.
    """
    var = data.variables["WT"][0, :, :, :]
    add_contourf(
        data,
        np.mean(var, axis=0),
        axes,
        "Vitesse du vent vertical\nmoyennée sur la verticale (m/s)",
        cmap="seismic",
        norm=TwoSlopeNorm(0),
    )


def wind_w_cross(data: Dataset, axes: plt.Axes, i_min: int, j_min: int, i_max: int, j_max: int):
    """
    Plot vertical wind speed cross section.

    Parameters
    ----------
    data : Dataset
        The Meso-NH dataset.
    axes : plt.Axes
        The axes to draw on.
    i_min : int
        The minimum index on x-axis.
    i_max : int
        The maximum index on x-axis.
    j_min : int
        The minimum index on y-axis.
    j_max : int
        The minimum index on y-axis.
    """
    wt = data.variables["WT"][0, :, :, :]
    wt_cross = oblique_proj(
        wt,
        data.variables["ni"][:],
        data.variables["nj"][:],
        data.variables["level"][:],
        i_min,
        j_min,
        i_max,
        j_max,
    )[1]

    ut = data.variables["UT"][0, :, :, :]
    ut = ut - np.mean(Dataset("CORSE.1.SEG01.019.vars.nc").variables["UT"][0, 28, 0:1500, 0:1250])

    ut_cross = oblique_proj(
        ut,
        data.variables["ni"][:],
        data.variables["nj"][:],
        data.variables["level"][:],
        i_min,
        j_min,
        i_max,
        j_max,
    )[1]

    add_cross_contourf(
        data,
        wt,
        axes,
        "Vitesse du vent vertical (m/s)",
        i_min,
        j_min,
        i_max,
        j_max,
        cmap="seismic",
        norm=TwoSlopeNorm(0),
    )

    mesh = ut_cross.shape[1] // 20
    axes.quiver(
        np.array(range(0, ut_cross.shape[1], mesh)),
        np.array(range(0, ut_cross.shape[0], 5)),
        ut_cross[::5, ::mesh],
        wt_cross[::5, ::mesh],
    )

    var = (
        data.variables["RCT"][0, :, :, :]
        + data.variables["RIT"][0, :, :, :]
        + data.variables["RST"][0, :, :, :]
        + data.variables["RGT"][0, :, :, :]
    ) * 1e3
    var_cross = oblique_proj(
        var,
        data.variables["ni"][:],
        data.variables["nj"][:],
        data.variables["level"][:],
        i_min,
        j_min,
        i_max,
        j_max,
    )[1]
    axes.contour(
        np.array(range(var_cross.shape[1])),
        np.array(range(var_cross.shape[0])),
        var_cross,
        levels=[0.1],
        colors=["black"],
    )


def rain_cross(data: Dataset, axes: plt.Axes, i_min: int, j_min: int, i_max: int, j_max: int):
    """
    Plot rain cross section.

    Parameters
    ----------
    data : Dataset
        The Meso-NH dataset.
    axes : plt.Axes
        The axes to draw on.
    i_min : int
        The minimum index on x-axis.
    i_max : int
        The maximum index on x-axis.
    j_min : int
        The minimum index on y-axis.
    j_max : int
        The minimum index on y-axis.
    """
    # Rain
    var = data.variables["RRT"][0, :, :, :] * 1e3
    add_cross_contourf(
        data,
        var,
        axes,
        "Rapport de mélange en eau précipitante (g/kg)",
        i_min,
        j_min,
        i_max,
        j_max,
        cmap="Blues",
    )

    # Clouds contour
    var = (
        data.variables["RCT"][0, :, :, :]
        + data.variables["RIT"][0, :, :, :]
        + data.variables["RST"][0, :, :, :]
        + data.variables["RGT"][0, :, :, :]
    ) * 1e3
    var_cross = oblique_proj(
        var,
        data.variables["ni"][:],
        data.variables["nj"][:],
        data.variables["level"][:],
        i_min,
        j_min,
        i_max,
        j_max,
    )[1]
    axes.contour(
        np.array(range(var_cross.shape[1])),
        np.array(range(var_cross.shape[0])),
        var_cross,
        levels=[0.1],
        colors=["black"],
    )


def snow_cross(data: Dataset, axes: plt.Axes, i_min: int, j_min: int, i_max: int, j_max: int):
    """
    Plot snow cross section.

    Parameters
    ----------
    data : Dataset
        The Meso-NH dataset.
    axes : plt.Axes
        The axes to draw on.
    i_min : int
        The minimum index on x-axis.
    i_max : int
        The maximum index on x-axis.
    j_min : int
        The minimum index on y-axis.
    j_max : int
        The minimum index on y-axis.
    """
    var = data.variables["RST"][0, :, :, :] * 1e3
    add_cross_contourf(
        data,
        var,
        axes,
        "Rapport de mélange en neige (g/kg)",
        i_min,
        j_min,
        i_max,
        j_max,
        cmap="Blues",
    )

    # Clouds contour
    var = (
        data.variables["RCT"][0, :, :, :]
        + data.variables["RIT"][0, :, :, :]
        + data.variables["RST"][0, :, :, :]
        + data.variables["RGT"][0, :, :, :]
    ) * 1e3
    var_cross = oblique_proj(
        var,
        data.variables["ni"][:],
        data.variables["nj"][:],
        data.variables["level"][:],
        i_min,
        j_min,
        i_max,
        j_max,
    )[1]
    axes.contour(
        np.array(range(var_cross.shape[1])),
        np.array(range(var_cross.shape[0])),
        var_cross,
        levels=[0.1],
        colors=["black"],
    )


def ice_cross(data: Dataset, axes: plt.Axes, i_min: int, j_min: int, i_max: int, j_max: int):
    """
    Plot ice cross section.

    Parameters
    ----------
    data : Dataset
        The Meso-NH dataset.
    axes : plt.Axes
        The axes to draw on.
    i_min : int
        The minimum index on x-axis.
    i_max : int
        The maximum index on x-axis.
    j_min : int
        The minimum index on y-axis.
    j_max : int
        The minimum index on y-axis.
    """
    var = data.variables["RIT"][0, :, :, :] * 1e3
    add_cross_contourf(
        data,
        var,
        axes,
        "Rapport de mélange en glace (g/kg)",
        i_min,
        j_min,
        i_max,
        j_max,
        cmap="Blues",
    )

    # Clouds contour
    var = (
        data.variables["RCT"][0, :, :, :]
        + data.variables["RIT"][0, :, :, :]
        + data.variables["RST"][0, :, :, :]
        + data.variables["RGT"][0, :, :, :]
    ) * 1e3
    var_cross = oblique_proj(
        var,
        data.variables["ni"][:],
        data.variables["nj"][:],
        data.variables["level"][:],
        i_min,
        j_min,
        i_max,
        j_max,
    )[1]
    axes.contour(
        np.array(range(var_cross.shape[1])),
        np.array(range(var_cross.shape[0])),
        var_cross,
        levels=[0.1],
        colors=["black"],
    )


def graupel_cross(data: Dataset, axes: plt.Axes, i_min: int, j_min: int, i_max: int, j_max: int):
    """
    Plot graupel cross section.

    Parameters
    ----------
    data : Dataset
        The Meso-NH dataset.
    axes : plt.Axes
        The axes to draw on.
    i_min : int
        The minimum index on x-axis.
    i_max : int
        The maximum index on x-axis.
    j_min : int
        The minimum index on y-axis.
    j_max : int
        The minimum index on y-axis.
    """
    var = data.variables["RGT"][0, :, :, :] * 1e3
    add_cross_contourf(
        data,
        var,
        axes,
        "Rapport de mélange en graupel (g/kg)",
        i_min,
        j_min,
        i_max,
        j_max,
        cmap="Blues",
    )

    # Clouds contour
    var = (
        data.variables["RCT"][0, :, :, :]
        + data.variables["RIT"][0, :, :, :]
        + data.variables["RST"][0, :, :, :]
        + data.variables["RGT"][0, :, :, :]
    ) * 1e3
    var_cross = oblique_proj(
        var,
        data.variables["ni"][:],
        data.variables["nj"][:],
        data.variables["level"][:],
        i_min,
        j_min,
        i_max,
        j_max,
    )[1]
    axes.contour(
        np.array(range(var_cross.shape[1])),
        np.array(range(var_cross.shape[0])),
        var_cross,
        levels=[0.1],
        colors=["black"],
    )


def cloudwater_cross(data: Dataset, axes: plt.Axes, i_min: int, j_min: int, i_max: int, j_max: int):
    """
    Plot cloud water cross section.

    Parameters
    ----------
    data : Dataset
        The Meso-NH dataset.
    axes : plt.Axes
        The axes to draw on.
    i_min : int
        The minimum index on x-axis.
    i_max : int
        The maximum index on x-axis.
    j_min : int
        The minimum index on y-axis.
    j_max : int
        The minimum index on y-axis.
    """
    var = data.variables["RCT"][0, :, :, :] * 1e3
    add_cross_contourf(
        data,
        var,
        axes,
        "Rapport de mélange en eau nuageuse (g/kg)",
        i_min,
        j_min,
        i_max,
        j_max,
        cmap="Blues",
    )

    # Clouds contour
    var = (
        data.variables["RCT"][0, :, :, :]
        + data.variables["RIT"][0, :, :, :]
        + data.variables["RST"][0, :, :, :]
        + data.variables["RGT"][0, :, :, :]
    ) * 1e3
    var_cross = oblique_proj(
        var,
        data.variables["ni"][:],
        data.variables["nj"][:],
        data.variables["level"][:],
        i_min,
        j_min,
        i_max,
        j_max,
    )[1]
    axes.contour(
        np.array(range(var_cross.shape[1])),
        np.array(range(var_cross.shape[0])),
        var_cross,
        levels=[0.1],
        colors=["black"],
    )


def tht_map(data: Dataset, axes: plt.Axes):
    """
    Plot potential temperature map.

    Parameters
    ----------
    data : Dataset
        The Meso-NH dataset.
    axes : plt.Axes
        The axes to draw on.
    """
    var = data.variables["THT"][0, :, :, :]
    add_contourf(
        data,
        np.mean(var, axis=0),
        axes,
        "Tempérture potentielle\nmoyennée sur la vertical (K)",
        cmap="turbo",
    )


def tht_cross(data: Dataset, axes: plt.Axes, i_min: int, j_min: int, i_max: int, j_max: int):
    """
    Plot potential temperature cross section

    Parameters
    ----------
    data : Dataset
        The Meso-NH dataset.
    axes : plt.Axes
        The axes to draw on.
    i_min : int
        The minimum index on x-axis.
    i_max : int
        The maximum index on x-axis.
    j_min : int
        The minimum index on y-axis.
    j_max : int
        The minimum index on y-axis.
    """
    var = data.variables["THT"][0, :, :, :]
    add_cross_contourf(
        data,
        var,
        axes,
        "Température potentielle (K)",
        i_min,
        j_min,
        i_max,
        j_max,
        cmap="turbo",
        levels=np.linspace(288, 333, 100),
        extend="both",
    )

    # Clouds contour
    var = (
        data.variables["RCT"][0, :, :, :]
        + data.variables["RIT"][0, :, :, :]
        + data.variables["RST"][0, :, :, :]
        + data.variables["RGT"][0, :, :, :]
    ) * 1e3
    var_cross = oblique_proj(
        var,
        data.variables["ni"][:],
        data.variables["nj"][:],
        data.variables["level"][:],
        i_min,
        j_min,
        i_max,
        j_max,
    )[1]
    axes.contour(
        np.array(range(var_cross.shape[1])),
        np.array(range(var_cross.shape[0])),
        var_cross,
        levels=[0.1],
        colors=["black"],
    )


def plot_figure_hydro(filename: str, resol_dx: int, i_min: int, j_min: int, i_max: int, j_max: int):
    """
    Show a cross-section of the different hydrometeor types and their distribution within the cloud.

    Parameters
    ----------
    filename : str
        The Meso-NH file to open.
    resol_dx : int
        The horizontal resolution (only for display purpose).
    i_min : int
        The minimum index on x-axis.
    i_max : int
        The maximum index on x-axis.
    j_min : int
        The minimum index on y-axis.
    j_max : int
        The minimum index on y-axis.
    """
    # Open data and get time
    data = Dataset(filename)
    time = datetime.strptime("2022-08-18 00:00:00", "%Y-%m-%d %H:%M:%S") + timedelta(
        seconds=int(data.variables["time"][0])
    )

    # Figure creation
    fig = plt.figure(figsize=(24, 10), layout="compressed")
    fig.suptitle(f"Meso-NH simulation on {time} (DX = {resol_dx} m)\nVertical cross section")
    gridspec = plt.GridSpec(2, 3, figure=fig)

    # Draw the cross section on the other fig
    axes = fig.add_subplot(gridspec[0, 0])
    rain_cross(data, axes, i_min, j_min, i_max, j_max)

    axes = fig.add_subplot(gridspec[0, 1])
    snow_cross(data, axes, i_min, j_min, i_max, j_max)

    axes = fig.add_subplot(gridspec[1, 0])
    ice_cross(data, axes, i_min, j_min, i_max, j_max)

    axes = fig.add_subplot(gridspec[1, 1])
    graupel_cross(data, axes, i_min, j_min, i_max, j_max)

    axes = fig.add_subplot(gridspec[0, 2])
    cloudwater_cross(data, axes, i_min, j_min, i_max, j_max)


def plot_figure(
    filename: str,
    resol_dx: int,
    func_map: Callable,
    func_cross: Callable,
    i_min: int,
    j_min: int,
    i_max: int,
    j_max: int,
):
    """
    Plot a cross-section and the associated map.

    Parameters
    ----------
    filename : str
        The Meso-NH file to open.
    resol_dx : int
        The horizontal resolution (only for display purpose).
    func_map : Callable
        The function to be used to plot the map.
    func_cross : Callable
        The function to be used to plot the cross section.
    i_min : int
        The minimum index on x-axis.
    i_max : int
        The maximum index on x-axis.
    j_min : int
        The minimum index on y-axis.
    j_max : int
        The minimum index on y-axis.
    """
    # Open data and get time
    data = Dataset(filename)
    time = datetime.strptime("2022-08-18 00:00:00", "%Y-%m-%d %H:%M:%S") + timedelta(
        seconds=int(data.variables["time"][0])
    )

    # Figure creation
    fig = plt.figure(layout="compressed", figsize=(16, 10))
    # fig.suptitle(f"Meso-NH simulation on {time} (DX = {resol_dx} m)\nVertical cross section")

    # Draw cross section on a map
    axes = fig.add_subplot(1, 2, 1, projection=ccrs.PlateCarree())
    axes.add_feature(cfeature.COASTLINE, color="gray")
    axes.add_feature(cfeature.BORDERS, color="gray")
    glines = axes.gridlines(draw_labels=True, linewidth=0.4, alpha=0.5)
    glines.top_labels = glines.right_labels = False
    func_map(data, axes)

    # Extract lon/lat cross section
    lon_min = data.variables["longitude"][j_min, i_min]
    lat_min = data.variables["latitude"][j_min, i_min]
    lon_max = data.variables["longitude"][j_max, i_max]
    lat_max = data.variables["latitude"][j_max, i_max]

    # Display cross section
    axes.plot([lon_min, lon_max], [lat_min, lat_max], color="red", linewidth=2)
    plt.text(lon_min, lat_min, "A", color="red", fontsize=15)
    plt.text(lon_max, lat_max, "B", color="red", fontsize=15)

    # Draw the cross section on the other fig
    axes = fig.add_subplot(1, 2, 2)
    func_cross(data, axes, i_min, j_min, i_max, j_max)


def quick_view(filename: str, i_min: int, j_min: int, i_max: int, j_max: int):
    """
    Display a map with the cross section on it to quickly see where the cross section will pass.

    Parameters
    ----------
    filename : str
        The Meso-NH file to open.
    i_min : int
        The minimum index on x-axis.
    i_max : int
        The maximum index on x-axis.
    j_min : int
        The minimum index on y-axis.
    j_max : int
        The minimum index on y-axis.
    """

    plt.figure()
    data = Dataset(filename)

    # var = np.sum((
    #     data.variables["RIT"][0, :, :, :]
    #     + data.variables["RST"][0, :, :, :]
    #     + data.variables["RGT"][0, :, :, :]
    # ) * 1e3, axis=0)

    # var = data.variables["RVT"][0, 42, :, :] * 1e3

    var = np.sqrt(data.variables["UT"][0, 28, :, :] ** 2 + data.variables["VT"][0, 28, :, :] ** 2)

    plt.plot([0, 0, 1250, 1250, 0], [0, 1500, 1500, 0, 0], color="red")

    contourf = plt.contourf(
        var,
        levels=np.linspace(var.min(), var.max(), 100),
        cmap=LinearSegmentedColormap.from_list("", ["black", "white", "blue", "red"]),
    )
    plt.colorbar(contourf)

    plt.quiver(
        np.arange(0, 2050, 50),
        np.arange(0, 2050, 50),
        data.variables["UT"][0, 28, ::50, ::50],
        data.variables["VT"][0, 28, ::50, ::50],
    )

    plt.plot([i_min, i_max], [j_min, j_max], color="red", linewidth=2)
    plt.text(i_min, j_min, "A", color="red", fontsize=15)
    plt.text(i_max, j_max, "B", color="red", fontsize=15)


if __name__ == "__main__":
    i_min, j_min = 350, 614
    i_max, j_max = 1100, 614
    filename = "CORSE.1.SEG01.004.vars.nc"

    # quick_view(filename, i_min, j_min, i_max, j_max)

    # plot_figure_hydro(filename, 250, clouds_map, i_min, j_min, i_max, j_max)
    # plot_figure(filename, 250, wind_w_map, wind_w_cross, i_min, j_min, i_max, j_max)
    # plot_figure(filename, 250, tht_map, tht_cross, i_min, j_min, i_max, j_max)
    plot_figure(filename, 250, rvt_map, rh_cross, i_min, j_min, i_max, j_max)
    # plot_figure(filename, 250, rvt_map, rvt_cross, i_min, j_min, i_max, j_max)
    # plot_figure(filename, 250, tket_map, clouds_cross, i_min, j_min, i_max, j_max)

    plt.show()
