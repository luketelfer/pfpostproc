import numpy as np
import xarray as xr
from pfpostproc.attrs import *

def overland_flow(pf,domain):
    # variables
    press = pf.press
    dx = domain.dx
    n = domain.mannings
    sy = domain.slopey
    sx = domain.slopex
    # calculations
    overlandflow = (press.where(press > 0, 0).isel(z = -1) ** (5/3)) * \
                   (dx / n) * \
                   np.sqrt(np.sqrt(sy**2 + sx**2))
    # dataarray
    overlandflow = xr.DataArray(
        overlandflow,
        dims = ['t','y','x'],
        name = 'overlandflow',
        attrs = overlandflow_attrs)
    # return
    return overlandflow

def surface_runoff(overlandflow,domain):
    # variables
    mask = domain.mask
    sy = domain.slopey
    sx = domain.slopex
    # calculations
    inverted_mask = mask.where(mask == 0).where(mask == 1, 0).isel(z = -1)
    surfrun = overlandflow.where(np.isnan(overlandflow), 0)
    surfrun = surfrun.where(np.isnan((surfrun * inverted_mask.shift(x = 1)) + sx.where(sx > 0)), 1) + \
              surfrun.where(np.isnan((surfrun * inverted_mask.shift(x = -1)) + sx.where(sx < 0)), 1) + \
              surfrun.where(np.isnan((surfrun * inverted_mask.shift(y = 1)) + sy.where(sy > 0)), 1) + \
              surfrun.where(np.isnan((surfrun * inverted_mask.shift(y = -1)) + sy.where(sy < 0)), 1)
    surfrun = surfrun * overlandflow
    # dataarray
    surfrun = xr.DataArray(
        surfrun,
        dims = ['t','y','x'],
        name = 'surfrun',
        attrs = surfrun_attrs)
    # return
    return surfrun

def init_surface_storage(domain):
    # variables
    init_press = domain.init_press.isel(z=-1)
    dx = domain.dx
    dy = domain.dy
    mask = domain.mask.isel(z=-1)
    # calculations
    init_surfstor = init_press.where(init_press>0,0)*dx*dy*mask
    # dataarray
    init_surfstor = xr.DataArray(
        init_surfstor,
        dims = ['y','x'],
        name = 'init_surfstor',
        attrs = init_surfstor_attrs)
    # return
    return init_surfstor

def surface_storage(pf,domain):
    # variables
    press = pf.press
    dx = domain.dx
    dy = domain.dy
    mask = domain.mask.isel(z=-1)
    # calculations
    surfstor = press.where(press>0,0).isel(z=-1)*dx*dy*mask
    # dataarray
    surfstor = xr.DataArray(
        surfstor,
        dims = ['t','y','x'],
        name = 'surfstor',
        attrs = surfstor_attrs)
    # return
    return surfstor

def init_subsurface_storage(domain):
    # variables
    init_press = domain.init_press
    init_satur = domain.init_satur
    por = domain.porosity
    spstor = domain.specific_storage
    dx = domain.dx
    dy = domain.dy
    dz = domain.dz * domain.dz_mult
    # calculations
    init_subsurfstor = ((init_satur * por) + (init_press * spstor)) * dx * dy * dz
    # dataarray
    init_subsurfstor = xr.DataArray(
        init_subsurfstor,
        dims = ['z','y','x'],
        name = 'init_subsurfstor',
        attrs = init_subsurfstor_attrs)
    # return
    return init_subsurfstor

def subsurface_storage(pf,domain):
    # variables
    press = pf.press
    satur = pf.satur
    por = domain.porosity
    spstor = domain.specific_storage
    dx = domain.dx
    dy = domain.dy
    dz = domain.dz * domain.dz_mult
    # calculations
    subsurfstor = ((satur * por) + (press * spstor)) * dx * dy * dz
    # dataarray
    subsurfstor = xr.DataArray(
        subsurfstor,
        dims = ['t','z','y','x'],
        name = 'subsurfstor',
        attrs = subsurfstor_attrs)
    # return
    return subsurfstor

def init_groundwater_storage(init_subsurfstor,domain):
    # variables
    init_satur = domain.init_satur
    # calculations
    init_gwstor = init_subsurfstor * init_satur.where(init_satur == 1)
    # dataarray
    init_gwstor = xr.DataArray(
        init_gwstor,
        dims = ['z','y','x'],
        name = 'init_gwstor',
        attrs = init_gwstor_attrs)
    # return
    return init_gwstor

def groundwater_storage(subsurfstor,pf):
    # variables
    satur = pf.satur
    # calculations
    gwstor = subsurfstor * satur.where(satur == 1)
    # dataarray
    gwstor = xr.DataArray(
        gwstor,
        dims = ['t','z','y','x'],
        name = 'gwstor',
        attrs = gwstor_attrs)
    # return
    return gwstor

def init_water_table_depth(init_gwstor,domain):
    # variables
    dz = domain.dz
    nz = domain.nz
    vdz = domain.dz_mult
    mask = domain.mask
    # calculations
    init_wtdepth = (dz * nz) - (init_gwstor.where(np.isnan(init_gwstor), 1) * dz * vdz).sum(dim='z') * mask.isel(z = -1)
    # dataarray
    init_wtdepth = xr.DataArray(
        init_wtdepth,
        dims = ['y','x'],
        name = 'init_wtdepth',
        attrs = init_wtdepth_attrs)
    # return
    return init_wtdepth

def water_table_depth(gwstor,domain):
    # variables
    dz = domain.dz
    nz = domain.nz
    vdz = domain.dz_mult
    mask = domain.mask
    # calculations
    wtdepth = (dz * nz) - (gwstor.where(np.isnan(gwstor), 1) * dz * vdz).sum(dim='z') * mask.isel(z = -1)
    # dataarray
    wtdepth = xr.DataArray(
        wtdepth,
        dims = ['t','y','x'],
        name = 'wtdepth',
        attrs = wtdepth_attrs)
    # return
    return wtdepth

def total_surface_storage(surfstor):
    # calculations
    tsurf = surfstor.sum(dim=['y','x'])
    # dataarray
    tsurf = xr.DataArray(
        tsurf,
        dims = ['t'],
        name = 'tsurf',
        attrs = tsurf_attrs)
    # return
    return tsurf

def total_subsurface_storage(subsurfstor):
    # calculations
    tsubsurf = subsurfstor.sum(dim=['z','y','x'])
    # dataarray
    tsubsurf = xr.DataArray(
        tsubsurf,
        dims = ['t'],
        name = 'tsubsurf',
        attrs = tsubsurf_attrs)
    # return
    return tsubsurf

def total_storage(tsurf,tsubsurf):
    # calculations
    tstor = tsurf + tsubsurf
    # dataarray
    tstor = xr.DataArray(
        tstor,
        dims = ['t'],
        name = 'tstor',
        attrs = tstor_attrs)
    # return
    return tstor

def change_in_surface_storage(tsurf,init_surfstor):
    # variables
    init_tsurf = init_surfstor.sum().expand_dims({'t':[0]})
    tsurf = xr.concat([init_tsurf,tsurf],dim='t')
    # calculations
    dsurf = tsurf.diff(dim='t')
    # dataarray
    dsurf = xr.DataArray(
        dsurf,
        dims = ['t'],
        name = 'dsurf',
        attrs = dsurf_attrs)
    # return
    return dsurf

def change_in_subsurface_storage(tsubsurf,init_subsurfstor):
    # variables
    init_tsubsurf = init_subsurfstor.sum().expand_dims({'t':[0]})
    tsubsurf = xr.concat([init_tsubsurf,tsubsurf],dim='t')
    # calculations
    dsubsurf = tsubsurf.diff(dim='t')
    # dataarray
    dsubsurf = xr.DataArray(
        dsubsurf,
        dims = ['t'],
        name = 'dsubsurf',
        attrs = dsubsurf_attrs)
    # return
    return dsubsurf

def change_in_total_storage(tstor,init_surfstor,init_subsurfstor):
    # variables
    init_tsurf = init_surfstor.sum().expand_dims({'t':[0]})
    init_tsubsurf = init_subsurfstor.sum().expand_dims({'t':[0]})
    init_tstor = init_tsurf + init_tsubsurf
    tstor = xr.concat([init_tstor,tstor],dim='t')
    # calculations
    dstor = tstor.diff(dim='t')
    # dataarray
    dstor = xr.DataArray(
        dstor,
        dims = ['t'],
        name = 'dstor',
        attrs = dstor_attrs)
    # return
    return dstor

def constant_precip(pf,domain,precip,timedump):
    # variables
    template = pf.overlandsum.where(np.isnan(pf.overlandsum),1)
    dx = domain.dx
    dy = domain.dy
    # calculations
    precip = template * precip * dx * dy * timedump
    # dataarray
    precip = xr.DataArray(
        precip,
        dims = ['t','y','x'],
        name = 'precip',
        attrs = precip_attrs)
    # return
    return precip

def total_precip(precip):
    # calculations
    tprecip = precip.sum(dim=['y','x'])
    # dataarray
    tprecip = xr.DataArray(
        tprecip,
        dims = ['t'],
        name = 'tprecip',
        attrs = tprecip_attrs)
    # return
    return tprecip

def total_evaptranssum(pf):
    # variables
    evaptranssum = pf.evaptranssum
    # calculations
    tevaptranssum = evaptranssum.sum(dim=['z','y','x'])
    # dataarray
    tevaptranssum = xr.DataArray(
        tevaptranssum,
        dims = ['t'],
        name = 'tevaptranssum',
        attrs = tevaptranssum_attrs)
    # return
    return tevaptranssum

def total_overlandsum(pf):
    # variables
    overlandsum = pf.overlandsum
    # calculations
    toverlandsum = overlandsum.sum(dim=['y','x'])
    # dataarray
    toverlandsum = xr.DataArray(
        toverlandsum,
        dims = ['t'],
        name = 'toverlandsum',
        attrs = toverlandsum_attrs)
    # return
    return toverlandsum

def pf_total_flux(tprecip,pf):
    # variables
    toverlandsum = pf.toverlandsum
    # calculations
    tflux = tprecip - toverlandsum
    # dataarray
    tflux = xr.DataArray(
        tflux,
        dims = ['t'],
        name = 'tflux',
        attrs = tflux_attrs)
    # return
    return tflux

def pfclm_total_flux(tevaptranssum,toverlandsum):
    # calculations
    tflux = tevaptranssum - toverlandsum
    # dataarray
    tflux = xr.DataArray(
        tflux,
        dims = ['t'],
        name = 'tflux',
        attrs = tflux_attrs)
    # return
    return tflux

def expected_storage(tstor,tflux,init_surfstor,init_subsurfstor):
    # variables
    init_tsurf = init_surfstor.sum().expand_dims({'t':[0]})
    init_tsubsurf = init_subsurfstor.sum().expand_dims({'t':[0]})
    init_tstor = init_tsurf + init_tsubsurf
    tstor = xr.concat([init_tstor,tstor],dim='t')
    # calculations
    expstor = tstor.shift(t=1) + tflux
    # dataarray
    expstor = xr.DataArray(
        expstor,
        dims = ['t'],
        name = 'expstor',
        attrs = expstor_attrs)
    # return
    return expstor

def water_balance_error(dstor,tflux):
    # calculations
    err = dstor - tflux
    # dataarray
    err = xr.DataArray(
        err,
        dims = ['t'],
        name = 'err',
        attrs = err_attrs)
    # return
    return err

def water_balance_perc_error(err,expstor):
    # calculations
    perc_err = (abs(err) / expstor) * 100
    # dataarray
    perc_err = xr.DataArray(
        perc_err,
        dims = ['t'],
        name = 'perc_err',
        attrs = perc_err_attrs)
    # return
    return perc_err