import numpy as np
from netCDF4 import Dataset



def write_ept(hour, ncid):
    data = Dataset(f"DX250_{hour}Z_EPT_LVL28.nc", "w", format="NETCDF4")
    
    data.createDimension("time", 1)
    data.createDimension("level", 1)
    data.createDimension("ni", 2050)
    data.createDimension("nj", 2050)

    time = data.createVariable("time", "d", ("time",))
    level = data.createVariable("level", "d", ("level",))
    lon = data.createVariable("longitude", "f8", ("nj", "ni"))
    lat = data.createVariable("latitude", "f8", ("nj", "ni"))
    ept = data.createVariable("equivalent_potential_temperature", "f8", ("time", "level", "nj", "ni"))
    ept.unit = "K"


    ref = Dataset(f"CORSE.1.SEG01.{ncid}.vars.nc")
    with open(f"DX250_{hour}Z_EPT_LVL28.npy", "rb") as file:
        ept_data = np.load(file)
    time[:] = ref.variables["time"][:]
    level[:] = ref.variables["level"][28]
    lon[:, :] = ref.variables["longitude"][:, :]
    lat[:, :] = ref.variables["latitude"][:, :]
    ept[:, :] = np.array([ept_data])
    
    data.close()


args = (
    ("0500", "004"),
    ("0630", "010"),
    ("0700", "012"),
    ("0715", "013"),
    ("0815", "017"),
    ("0845", "019")
)

for hour, ncid in args:
    write_ept(hour, ncid)

