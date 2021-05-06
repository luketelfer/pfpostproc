from glob import glob
from os.path import join
import xarray as xr
from pfpostproc.func import *
from pfpostproc.attrs import *
import pfpostproc.calc as calc
from datetime import datetime,timedelta
xr.set_options(keep_attrs=True)

def get_domain(rundir,runname):
    
    # get files
    fdom = rundir + runname + '.out.00000.nc'
    
    # create dataset
    domain = xr.open_dataset(fdom,drop_variables=['perm_x','perm_y'],chunks={})
    domain = domain.rename({
        'perm_z': 'perm',
        'DZ_Multiplier': 'dz_mult',
        'pressure': 'init_press',
        'saturation': 'init_satur'
    })
    domain = domain.squeeze(drop=True)
    domain['mask'] = domain.mask.where(domain.mask==0,1)
    
    # mask everything
    domain = domain * domain.mask.where(domain.mask==1).isel(z=-1)
    
    # add attrs
    domain.attrs = domain_attrs
    domain.init_press.attrs = init_press_attrs
    domain.init_satur.attrs = init_satur_attrs
    domain.mannings.attrs = mannings_attrs
    domain.perm.attrs = perm_attrs
    domain.porosity.attrs = porosity_attrs
    domain.specific_storage.attrs = specific_storage_attrs
    domain.slopex.attrs = slopex_attrs
    domain.slopey.attrs = slopey_attrs
    domain.dz_mult.attrs = dz_mult_attrs
    domain.mask.attrs = mask_attrs
    
    # add coords
    domain['x'] = np.arange(domain.nx) + 1
    domain['y'] = np.arange(domain.ny) + 1
    domain['z'] = np.arange(domain.nz) + 1
    
    return domain

def get_pf(rundir,runname,domain):
    
    # get files
    fout = sorted(glob(rundir + runname + '.out.[!CLM]*.nc'))[1:]
    findi = rundir + runname + '.indicator.pfb'
    
    # create dataset
    ds = xr.open_mfdataset(fout,chunks={'time':730})
    ds = ds.rename({
        'pressure': 'press',
        'saturation': 'satur',
        'evaptrans_sum': 'evaptranssum',
        'overland_sum': 'overlandsum',
        'time': 't'
    })

    # add indicator file
    ds['indicator'] = xr.DataArray(pfb_arr(findi),dims=['z','y','x']).chunk()
    
    # apply variable dz
    ds['evaptrans'] = ds.evaptrans * domain.dz_mult.values
    ds['evaptranssum'] = ds.evaptranssum * domain.dz_mult.values
    
    # mask everything
    ds = ds * domain.mask.where(domain.mask==1).isel(z=-1).values
    
    # add attrs
    ds.attrs = domain_attrs
    ds.press.attrs = press_attrs
    ds.satur.attrs = satur_attrs
    ds.evaptrans.attrs = evaptrans_attrs
    ds.evaptranssum.attrs = evaptranssum_attrs
    ds.overlandsum.attrs = overlandsum_attrs
    ds.indicator.attrs = indicator_attrs
    
    # add coords
    ds['x'] = np.arange(ds.nx) + 1
    ds['y'] = np.arange(ds.ny) + 1
    ds['z'] = np.arange(ds.nz) + 1
    
    return ds

def get_clm(rundir,runname,domain):
    
    # get files
    fclm = sorted(glob(rundir + runname + '.out.CLM.*.nc'))
    
    # create dataset
    clm = xr.open_mfdataset(fclm,chunks={'time':730})
    t_soil = clm.t_soil.pad({'z':(0,15)}).chunk({'z':1})
    clm = clm.drop_vars('t_soil').merge(t_soil)
    clm = clm.rename({'time':'t'})
    clm = clm.chunk({'t':730,'z':25})
    
    # mask everything
    clm = clm * domain.mask.where(domain.mask==1).isel(z=-1).values
    
    # unit conversions
    clm['qflx_evap_grnd'] = clm['qflx_evap_grnd'] * 3.6 * domain.dx * domain.dy
    clm['qflx_evap_soi'] = clm['qflx_evap_soi'] * 3.6 * domain.dx * domain.dy
    clm['qflx_evap_tot'] = clm['qflx_evap_tot'] * 3.6 * domain.dx * domain.dy
    clm['qflx_evap_veg'] = clm['qflx_evap_veg'] * 3.6 * domain.dx * domain.dy
    clm['qflx_infl'] = clm['qflx_infl'] * 3.6 * domain.dx * domain.dy
    clm['qflx_tran_veg'] = clm['qflx_tran_veg'] * 3.6 * domain.dx * domain.dy
    clm['swe_out'] = clm['swe_out'] * domain.dx
    
    # add attrs
    clm.attrs = domain_attrs
    clm.eflx_lh_tot.attrs = eflx_lh_tot_attrs
    clm.eflx_lwrad_out.attrs = eflx_lwrad_out_attrs
    clm.eflx_sh_tot.attrs = eflx_sh_tot_attrs
    clm.eflx_soil_grnd.attrs = eflx_soil_grnd_attrs
    clm.qflx_evap_grnd.attrs = qflx_evap_grnd_attrs
    clm.qflx_evap_soi.attrs = qflx_evap_soi_attrs
    clm.qflx_evap_tot.attrs = qflx_evap_tot_attrs
    clm.qflx_evap_veg.attrs = qflx_evap_veg_attrs
    clm.qflx_infl.attrs = qflx_infl_attrs
    clm.qflx_tran_veg.attrs = qflx_tran_veg_attrs
    clm.swe_out.attrs = swe_out_attrs
    clm.t_grnd.attrs = t_grnd_attrs
    clm.t_soil.attrs = t_soil_attrs
    
    # add coords
    clm['x'] = np.arange(clm.nx) + 1
    clm['y'] = np.arange(clm.ny) + 1
    clm['z'] = np.arange(clm.nz) + 1
    
    return clm

def get_calc(pf,domain):
    
    # calculations (part I)
    overlandflow = calc.overland_flow(pf,domain).chunk({'t':730})
    init_surfstor = calc.init_surface_storage(domain).chunk()
    init_subsurfstor = calc.init_subsurface_storage(domain).chunk()
    init_gwstor = calc.init_groundwater_storage(init_subsurfstor,domain).chunk()
    init_wtdepth = calc.init_water_table_depth(init_gwstor,domain).chunk()
    surfstor = calc.surface_storage(pf,domain).chunk({'t':730})
    subsurfstor = calc.subsurface_storage(pf,domain).chunk({'t':730})
    gwstor = calc.groundwater_storage(subsurfstor,pf).chunk({'t':730})
    wtdepth = calc.water_table_depth(gwstor,domain).chunk({'t':730})
    
    # create dataset (part)
    ds = xr.Dataset({
        'overlandflow': overlandflow,
        'init_surfstor': init_surfstor,
        'surfstor': surfstor,
        'init_subsurfstor': init_subsurfstor,
        'subsurfstor': subsurfstor,
        'init_gwstor': init_gwstor,
        'gwstor': gwstor,
        'init_wtdepth': init_wtdepth,
        'wtdepth': wtdepth
    })
    
    # calculations (part II)
    ds['tsurf'] = calc.total_surface_storage(ds.surfstor).chunk({'t':730})
    ds['tsubsurf'] = calc.total_subsurface_storage(ds.subsurfstor).chunk({'t':730})
    ds['tstor'] = calc.total_storage(ds.tsurf,ds.tsubsurf).chunk({'t':730})
    ds['dsurf'] = calc.change_in_surface_storage(ds.tsurf,ds.init_surfstor).chunk({'t':730})
    ds['dsubsurf'] = calc.change_in_subsurface_storage(ds.tsubsurf,ds.init_subsurfstor).chunk({'t':730})
    ds['dstor'] = calc.change_in_total_storage(ds.tstor,ds.init_surfstor,ds.init_subsurfstor).chunk({'t':730})
    ds['tevaptranssum'] = calc.total_evaptranssum(pf).chunk({'t':730})
    ds['toverlandsum'] = calc.total_overlandsum(pf).chunk({'t':730})
    ds['tflux'] = calc.pfclm_total_flux(ds.tevaptranssum,ds.toverlandsum).chunk({'t':730})
    ds['expstor'] = calc.expected_storage(ds.tstor,ds.tflux,ds.init_surfstor,ds.init_subsurfstor).chunk({'t':730})
    ds['err'] = calc.water_balance_error(ds.dstor,ds.tflux).chunk({'t':730})
    ds['perc_err'] = calc.water_balance_perc_error(ds.err,ds.expstor).chunk({'t':730})
    
    # add attrs
    ds.attrs = domain_attrs
    
    # add coords
    ds['x'] = np.arange(ds.nx) + 1
    ds['y'] = np.arange(ds.ny) + 1
    ds['z'] = np.arange(ds.nz) + 1
    
    return ds

def get_wrf(domain):
    
    # get files
    fwrf = '/home/ltelfer/scratch/pf_southForkSalmon/input_files/wrf_2006/wrf_forcings.nc'
    
    # create dataset
    wrf = xr.open_dataset(fwrf,chunks={})
    wrf = wrf.rename({
        'DSWR': 'wrf_shortwave',
        'DLWR': 'wrf_longwave',
        'APCP': 'wrf_precip',
        'Temp': 'wrf_temp',
        'UGRD': 'wrf_ewwind',
        'VGRD': 'wrf_nswind',
        'Press': 'wrf_atmospress',
        'SPFH': 'wrf_humid'
    })
    
    # mask everything
    wrf = wrf * domain.mask.isel(z=-1).values
    
    # convert precip units from mm/s to m3
    wrf['wrf_precip'] = wrf.wrf_precip * 3.6 * domain.dx * domain.dy
    
    # add attrs
    wrf.wrf_precip.attrs['units'] = 'm3'
    
    return wrf