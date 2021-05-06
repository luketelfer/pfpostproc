# pf output attrs
domain_attrs = {
    'dx': 1000,
    'dy': 1000,
    'dz': 40,
    'nx': 24,
    'ny': 64,
    'nz': 25
}
press_attrs = {
    'description': 'pressure head (fluid pressure / weight density)',
    'units':'m'
}
satur_attrs = {
    'description': 'proportion of volume containing water',
    'units':'dimensionless'
}
evaptrans_attrs = {
    'description': 'flux of water (rate) due to precip & et',
    'units':'1/hr'
}
evaptranssum_attrs = {
    'description': 'flux of water (volume) due to precip & et',
    'units':'m^3'
}
overlandsum_attrs = {
    'description': 'flux of water (volume) due to surface runoff',
    'units':'m^3'
}
init_press_attrs = {
    'description': 'initial pressure head',
    'units':'m'
}
init_satur_attrs = {
    'description': 'initial proportion of volume containing water',
    'units':'dimensionless'
}
mannings_attrs = {
    'description': 'coefficient representing surface roughness',
    'units':'hr/m^(1/3)'
}
perm_attrs = {
    'description': '(hydraulic conductivity) ratio of velocity to hydraulic gradient indicating permeability of porous media',
    'units':'m/hr'
}
porosity_attrs = {
    'description': 'proportion of volume occupied by void space',
    'units':'dimensionless'
}
specific_storage_attrs = {
    'description': 'water equivalent proportion of volume flux per unit change in head (relates to water and aquifer compressibility)',
    'units':'1/m'
}
slopex_attrs = {
    'description': 'slopes in the x direction',
    'units':'dimensionless'
}
slopey_attrs = {
    'description': 'slopes in the y direction',
    'units':'dimensionless'
}
dz_mult_attrs = {
    'description': 'thickness of each layer as proportion of domain dz',
    'units':'dimensionless'
}
mask_attrs = {
    'description': 'distinguishes active cells from inactive cells',
    'units':'none'
}
indicator_attrs = {
    'description': 'parameter set for each grid cell',
    '1': 'unburned, sandy loam',
    '2': 'low severity',
    '3': 'moderate severity',
    '4': 'high severity',
    '5': 'granite (saprolite)',
    '6': 'granite (intact)'
}
elev_attrs = {
    'description': 'surface elevation in meters above sea-level',
    'units':'m'
}

# calc attrs
overlandflow_attrs = {
    'description': 'flux of water (volume) out of each grid cell due to surface runoff',
    'equation': '''pf user's manual (p.75, [5.45])''',
    'units':'m^3'
}
surfrun_attrs = {
    'description': '(pftools) flux of water (volume) leaving the domain due to surface runoff',
    'equation': '''pf user's manual (p.75, [5.45])''',
    'units':'m^3'
}
init_surfstor_attrs = {
    'description': 'initial surface storage',
    'equation': '''pf user's manual (p.75, [5.44])''',
    'units':'m^3'
}
surfstor_attrs = {
    'description': 'volume of ponded water above each surface grid cell',
    'equation': '''pf user's manual (p.75, [5.44])''',
    'units':'m^3'
}
init_subsurfstor_attrs = {
    'description': 'initial subsurface storage',
    'equation': '''pf user's manual (p.75, [5.44])''',
    'units':'m^3'
}
subsurfstor_attrs = {
    'description': 'volume of ponded water above each surface grid cell',
    'equation': '''pf user's manual (p.75, [5.44])''',
    'units':'m^3'
}
init_gwstor_attrs = {
    'description': 'initial groundwater storage',
    'units':'m^3'
}
gwstor_attrs = {
    'description': 'volume of water stored in each fully saturated subsurface grid cell',
    'units':'m^3'
}
init_wtdepth_attrs = {
    'description': 'initial water table depth',
    'units':'m'
}
wtdepth_attrs = {
    'description': 'distance between the highest fully saturated grid cell and the domain surface',
    'units':'m'
}
init_tsurf_attrs = {
    'description': 'initial total surface storage',
    'units': 'm^3'
}
tsurf_attrs = {
    'description': 'total surface storage for each timestep',
    'units': 'm^3'
}
init_tsubsurf_attrs = {
    'description': 'initial total subsurface storage',
    'units': 'm^3'
}
tsubsurf_attrs = {
    'description': 'total subsurface storage for each timestep',
    'units': 'm^3'
}
tstor_attrs = {
    'description': 'initial total storage in the domain',
    'units': 'm^3'
}
tstor_attrs = {
    'description': 'total storage in the domain for each timestep (surface + subsurface)',
    'units': 'm^3'
}
dsurf_attrs = {
    'description': 'change in surface storage (t1 - t0)',
    'units': 'm^3'
}
dsubsurf_attrs = {
    'description': 'change in subsurface storage (t1 - t0)',
    'units': 'm^3'
}
dstor_attrs = {
    'description': 'change in total storage (t1 - t0)',
    'units': 'm^3'
}
tevaptranssum_attrs = {
    'description': 'total evapotranspiration/infiltration flux in(+)/out(-) of the domain for each timestep',
    'units': 'm^3'
}
toverlandsum_attrs = {
    'description': 'total runoff out of domain for each timestep',
    'units': 'm^3'
}
tflux_attrs = {
    'description': 'total net flux in(+)/out(-) of domain for each timestep',
    'units': 'm^3'
}
expstor_attrs = {
    'description': 'expected storage for each timestep based on previous storage and current flux',
    'units': 'm^3'
}
err_attrs = {
    'description': 'volumetric difference between the change in total storage and total flux',
    'units': 'm^3'
}
perc_err_attrs = {
    'description': 'water balance error as a percentage of the expected storage',
    'evaluation': 'percerr < 0.0005 for all timesteps means water balance is closing',
    'units': '%'
}


# clm output attrs
eflx_lh_tot_attrs = {
    'description': 'latent heat flux from canopy height to atmosphere',
    'units':'W/m^2'
}
eflx_lwrad_out_attrs = {
    'description': 'outgoing long-wave radiation from ground+canopy',
    'units':'W/m^2'
}
eflx_sh_tot_attrs = {
    'description': 'sensible heat from canopy height to atmosphere',
    'units':'W/m^2'
}
eflx_soil_grnd_attrs = {
    'description': 'ground heat flux',
    'units':'W/m^2'
}
qflx_evap_tot_attrs = {
    'description': 'evapotranspiration from canopy height to atmosphere',
    'units':'m^3'
}
qflx_evap_grnd_attrs = {
    'description': 'ground surface evaporation rate',
    'units':'m^3'
}
qflx_evap_soi_attrs = {
    'description': 'evaporation heat flux from ground',
    'units':'m^3'
}
qflx_evap_veg_attrs = {
    'description': 'evaporation+transpiration from leaves',
    'units':'m^3'
}
qflx_tran_veg_attrs = {
    'description': 'transpiration rate',
    'units':'m^3'
}
qflx_infl_attrs = {
    'description': 'infiltration',
    'units':'m^3'
}
swe_out_attrs = {
    'description': 'snow water equivalent',
    'units':'m^3'
}
t_grnd_attrs = {
    'description': 'ground surface temperature',
    'units':'K'
}
t_soil_attrs = {
    'description': 'soil temperature',
    'units':'K'
}

# usgs gages
shortwave_attrs = {
    'description': 'visible or shart-wave radiation',
    'units': 'W/m^2',
}
longwave_attrs = {
    'description': 'long-wave radiation',
    'units': 'W/m^2',
}
precip_attrs = {
    'description': 'precipitation',
    'units': 'm^3',
}
temp_attrs = {
    'description': 'air temperature',
    'units': 'K',
}
wind_ew_attrs = {
    'description': 'east-west wind speed',
    'units': 'm/s',
}
wind_ns_attrs = {
    'description': 'south-north wind speed',
    'units': 'm/s',
}
atm_press_attrs = {
    'description': 'atmospheric pressure',
    'units': 'Pa',
}
sp_humid_attrs = {
    'description': 'specific humidity',
    'units': 'kg/kg',
}