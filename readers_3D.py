"""
Readers
=======

Description
-----------
This module provides a number of reader classes for managing different types of time and variables.
These classes are all structured in the same way, so that they can be used in the same way
regardless of the files or variables managed.

It also gives some basic functions to manipulate index and conversion from lon/lat to index and
index to lat/lon.

Create a reader class
---------------------
To create a reader class, you should have the following attributes:

* files     : the list of files

* data      : the current open file

* longitude : the longitude of the measures

* latitude  : the latitude of the measures

and the following methods:

* get_data(file_index: int) to open a file

* get_var(*varnames, func: Callable = lambda x: x) to get a variable

* get_limits(*varnames, func: Callable = lambda x: x) to get the min and max of a variable
"""

from collections.abc import Callable, Iterable

import numpy as np
from netCDF4 import Dataset


class MesoNH:
    """
    MesoNH is a reader class that reads Meso-NH files.

    Attributes
    ----------
    files : list
        The list of the Meso-NH files.
    data : Dataset
        The current file open.
    longitude : np.array
        The longitudes.
    latitude : np.array
        The latitudes.
    """

    def __init__(self, files):
        """Constructor method."""
        self.files = files
        self.data = None

        data = Dataset(self.files[0])
        self.longitude = data.variables["longitude"][:, :]
        self.latitude = data.variables["latitude"][:, :]
        self.level = data.variables["level"][:]

    def get_data(self, file_index: int):
        """
        Open the file corresponding to the given file_index. The
        Dataset object is stored into MesoNH.data.

        Parameters
        ----------
        file_index : int
            The index of the the file to open.
        """
        self.data = Dataset(self.files[file_index])

    def get_var(self, *varnames, level, func: Callable = lambda x: x):
        """
        Returns the NumPy array corresponding to result given by func applied on the given
        variables.

        Parameters
        ----------
        varnames : str
            Variable names to give to func.
        level : int, keyword-only
            Level of altitude at which the variables are observed 
        func : Callable, keyword-only, optionnal
            The function to apply to the given variables.

        Returns
        -------
        out : np.array
            The result given by func.
        """
        args = []
        for varname in varnames:
            args.append(self.data.variables[varname][0, level, :, :])

        return func(*args)

    def get_stats(
        self,
        index_i: int,
        index_j: int,
        *varnames,
        func: Callable = lambda x: x,
        t_range: iter = None,
        size: int = 1,
    ):
        """
        Compute spatio-temporal limits average and standard deviation over a group of pixel centered
        on a given position and over a given time range.

        .. note:
            The standard deviation is calculated for each time step and it returns the average
            standard deviation.

        Parameters
        ----------
        index_i : int
            The index on the x-axis.
        index_j : int
            The index on the y-axis
        varnames : str
            Variable names to give to func.
        func : Callable, keyword-only, optionnal
            The function to apply to the given variables.
        t_range : iter, keyword-only, optionnal
            The temporal range over wich the average is to be calculated. By default all the
            available time interval will be taken.
            You can give a range or the list of the index you want.
        size : int, keyword-only, optionnal
            The size of the spatial average in index. By default it's set on 1, so the average
            will be calculated of the given case and over the cases around in the four
            directions:

                      size ├───┤
                ───┌───┬───┬───┐ ┬
                j+1│   │   │   │ │ size
                ───├───┼───┼───┤ ┴
                 j │   │i;j│   │
                ───├───┼───┼───┤
                j-1│   │   │   │
                ───└───┴───┴───┘
                   │i-1│ i │i+1│

        Returns
        -------
        out : tuple
            A tuple that contains the limits, the average and the standard deviation.
        """
        if not t_range:
            t_range = range(0, len(self.files))

        mean_per_timestep = []
        std_per_timestep = []
        var_min, var_max = np.inf, -np.inf
        for i in t_range:
            data = Dataset(self.files[i])

            args = []
            for varname in varnames:
                args.append(data.variables[varname][0])
            array = func(*args)
            array = array[index_j - size : index_j + size, index_i - size : index_i + size]

            var_min = min(var_min, array.min())
            var_max = max(var_max, array.max())

            mean_per_timestep.append(np.mean(array))
            std_per_timestep.append(np.std(array))

        return (var_min, var_max), np.mean(mean_per_timestep), np.mean(std_per_timestep)

    def get_limits(self, *varnames, func: Callable = lambda x: x):
        """
        Search min and max of a given variable.

        Parameters
        ----------
        varnames : str
            The names of the variables.
        func : Callable, keyword-only, optionnal
            The function to apply to the given variables.

        Returns
        -------
        out : tuple
            A tuple containing two elements: (var_min, var_max).
        """
        var_min = np.inf
        var_max = -np.inf
        for _, file in enumerate(self.files):
            data = Dataset(file)

            args = []
            for varname in varnames:
                args.append(data.variables[varname][0])
            result = func(*args)

            current_min = result.min()
            current_max = result.max()

            var_min = min(var_min, current_min)
            var_max = max(var_max, current_max)

        return var_min, var_max


def get_mesonh(resolution_dx: int, path: str = "../Donnees/"):
    """
    Returns a Meso-NH reader for the given resolution and path.

    Parameters
    ----------
    resolution_dx : int
        The wanted resolution.
    path : str, optionnal
        The path of the netCDF files. By default it's set on ../Donnees/.

    Returns
    -------
    out : MesoNH
        The Meso-NH reader.
    """
    files = [
        f"{path}DX{resolution_dx}/CORSE.1.SEG01.{str(t).zfill(3)}.vars.nc" for t in (4, 10, 12, 13, 17, 19)
    ]

    return MesoNH(files)


def get_index(array: np.array, target: float):
    """
    Search and return the index of the value closest to target in the given array. This function
    can handle n-dimensionnal arrays.

    Parameters
    ----------
    array : np.array
        The array in which to search.
    target : float
        The value to search in array.

    Returns
    -------
    out : np.array
        The index of the value closest to target in array. If seraval indexes matche, it
        only returns the first one.
    """
    delta = np.abs(target - array)
    index = np.array(np.where(delta == delta.min()))
    return index[:, 0]


def get_index_from_vect(x_array: np.array, y_array: np.array, value: Iterable[float, float]):
    """
    Search for the given vector value in the x and y-array.

    Parameters
    ----------
    x_array : np.array
        The components on the x-axis.
    y_array : np.array
        The components on the y-axis.
    value : Iterable[float, float]
        The vector to search for.

    Returns
    -------
    out : np.array
        The index on x- and y-axis.
    """
    norms = np.sqrt((x_array - value[0]) ** 2 + (y_array - value[1]) ** 2)
    index = np.array(np.where(norms == norms.min()))
    return index[:, 0]


def index_to_lonlat(reader, i: int, j: int):
    """
    Get the longitudes and latitudes from given limits indexes.

    Parameters
    ----------
    reader
        An instance of reader.
    i : innt
        The index on x-axis.
    j : int
        The index on y-axis.

    Returns
    -------
    out : tuple
        A tuple that contains two tuples: (longitude_min, longitude_max) and
        (latitude_min, latitude_max).
    """
    lon = reader.longitude[j, i]
    lat = reader.latitude[j, i]
    return (lon, lat)


def lonlat_to_index(reader, lon: float, lat: float):
    """
    Get the indexes from given limit longitudes and latitudes.

    .. warning::
        Due to projection, the returned indexes may not match perfectly match the given lon/lat.

    Parameters
    ----------
    reader
        An instance of reader.
    lon : tuple
        The longitude to search.
    lat : tuple
        The latitude to search.

    Returns
    -------
    out : tuple
        A tuple that contains two elements: (i, j).
    """
    j, i = get_index_from_vect(reader.longitude, reader.latitude, (lon, lat))
    return i, j


def get_time_index(hour: int, minute: int):
    """
    Compute the index of the Meso-NH file from hour and minute.

    Parameters
    ----------
    hour : int
        The hours.
    minute : int
        The minutes.

    Returns
    -------
    out : int
        The index of the file corresponding to the given timestamp.
    """
    return int((hour - 4) * 4 + (minute/60) * 4 - 4)
