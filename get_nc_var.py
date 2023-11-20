###
# Install first dask ,then xarray>=2023.0.6
###
FLOAT_LEN = 4  # bytes
THRESHOLD = 2000000000  # about 2GB
VARS = ["EXTCOF55","T2","ALBEDO"]
URL = 'https://stromboli.le.isac.cnr.it/opendap/fabio_dottorato/case_dust_june_2020/simulazioni/D_ATL_01/wrfout_d01_2020-06-04_000000.nc'
URL2 = 'https://stromboli.le.isac.cnr.it/opendap/fabio_dottorato/dati_satellitari/viirs/without_post_processing/AERDB_D3_VIIRS_SNPP.A2020156.002.2023083210306.nc'

import numpy as np
import xarray as xr
import netCDF4
from pydap.client import open_url
from pydap.cas.urs import setup_session
from functools import reduce
from math import ceil
import glob, os

class Piece:
   
    var = None
    n_piece = 0
    dims = None
    shape = None  # Tuple 
    concat_dim = ''  # 
    
    def __init__(self, offset=None, end=None):
        self.offset = []  # per ogni dimensione, è l'offset da cui far partire la richiesta di dati
        self.end = []  # per ogni dimensione, il punto finale della richiesta
        self.d = []
        for i in range(len(Piece.dims)):
            if offset == None:
                self.offset.append(0)
            else:
                self.offset.append(offset[i])
            if end == None:
                self.end.append(1)
            else:
                self.end.append(end[i])
    def size(self):
        _tot = 1.
        for i in range(len(Piece.dims)):
            _tot = _tot * ( self.end[i] - self.offset[i])
        _tot = _tot * FLOAT_LEN
        return _tot
    
    def ask_file(self):
        print(f"Ask opendap server ... -- offset:{self.offset} -- end:{self.end}\n")
        # Basandosi su offset e ends, chiede al server la variabile in questione
        ds = xr.open_dataset(URL, engine='pydap')  
        ds = ds[Piece.var]
        dimension_slices = {}
        for index, dim in enumerate(Piece.dims):
            # dimensions_slices.append(f"{dim}":(self.offset[index], self.end[index]))
            dimension_slices[dim]=(self.offset[index],self.end[index])
        selection_dict = {dim: slice(start, end+1) for dim,(start, end) in dimension_slices.items()}
        print(selection_dict)
        ds = ds.sel(selection_dict)
        filename = f"{Piece.var}_{Piece.n_piece}.nc"
        ds.to_netcdf(path=filename, format="NETCDF4")
        print (f"Temporary saved {filename}")
        print(ds.dims)
   
    def file_finished(self):
        # True se su tutte le dim , ends ha raggiunto il limite
        for i in range(len(Piece.dims)):
            if self.end[i] == Piece.shape[i]:
                pass
            else:
                return False
        return True

#####
####
###
##
#

ds_s = xr.open_dataset(URL, engine='pydap')  
for var in VARS:
    Piece.var = var
    Piece.dims = ds_s[var].dims
    Piece.shape = ds_s[var].shape
    print(f"var: {Piece.var}  -- dims:{Piece.dims}  --  shape:{Piece.shape}")
    lenght_tot = reduce(lambda x, y: x * y, Piece.shape, 1)
    if lenght_tot > THRESHOLD:
        # > 2 GB, bisogna partizionarla
        piece = Piece()
        for i in range(len(Piece.dims)):
            _concat_dim = Piece.dims[i]
            d = piece.offset[i]
            
            # Azzera tutte le altre dimensioni già scorse
            for ii in range(i+1,len(Piece.dims)):
                piece.offset[i] = 0
                piece.end[i]
            
            d_continua = True
            d = piece.offset[i]
            while d_continua:
                d_continua = False
                d += 1
                if d <= Piece.shape[i]:
                    d_continua = True
                    piece.end[i] = d
                    if piece.size() > THRESHOLD:
                        Piece.n_piece += 1
                        piece.ask_file()
                        piece.offset[i] = d
                else:
                    for t in range(len(Piece.dims)):
                        print(f"offset[{t}]:{piece.offset[t]}  end[{t}]:{piece.end[t]}")
                    # Raggiunto il limite sulla coordinata
                    d -= 1
                    piece.end[i] = d
        
        ###
        # Merge dei file temporanei .nc
        ###
        print('Merge temp *.nc files....')
        merge_ds = xr.open_mfdataset(f"{Piece.var}_*.nc", concat_dim=_concat_dim, combine='nested')
        merge_ds.to_netcdf(path=f"{Piece.var}.nc", format="NETCDF4")
        
        #Erase temp files
        print('Deleting temp files....')
        for f in glob.glob(f"{Piece.var}_*.nc"):
            os.remove(f)

    else:
        ###
        # Small file
        ###
        ds = xr.open_dataset(URL, engine='pydap')
        ds = ds[var]
        ds.to_netcdf(path=f"{var}.nc", format="NETCDF4")


exit()
