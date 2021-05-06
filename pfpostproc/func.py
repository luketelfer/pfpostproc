# packages
from parflowio.pyParflowio import PFData
import rasterio as rio
import numpy as np
import xarray as xr
from pfpostproc.attrs import *
import pandas as pd
from datetime import datetime,timedelta
from memory_profiler import profile

# open pfb, return array
def pfb_arr(fpath):
    pfb = PFData(fpath)
    pfb.loadHeader();
    pfb.loadData();
    arr = pfb.copyDataArray()
    arr = arr.squeeze()
    pfb.close()
    return arr

# open raster, return array
def raster_arr(fpath):
    raster = rio.open(fpath)
    arr = raster.read(1).astype(float)
    arr = np.flip(arr,axis=0)
    return arr
    