# Description

## pfpostproc.func

Module includes a couple of functions for adding the indicator pfb and elevation raster to my xarray dataset.

## pfpostproc.workflow

Module includes functions to load and/or assemble xarray datasets and adds attributes.  

1. get_domain

Loads dataset containing domain parameters and initial conditions (pressure, saturation).  This includes the variable dz weights.

2. get_pf

Loads dataset containing model outputs over time.  Note: evaptrans and evaptranssum needed to be multiplied by the variable dz in order to be accurate (only matters when running with CLM because otherwise these variables equal zero).

3. get_clm

Loads dataset containing CLM outputs over time and includes unit conversions that I wanted.  Note: none of these variables were used in my calculations. However, some might be needed for a complete water balance assessment when running with CLM.

Potentially helpful:

* qflx_evap_tol = qflx_evap_soi + qflx_evap_veg
* evaptranssum (from pf outputs) = qflx_infl - qflx_tran_veg (this is the flux between PF and CLM) 

4. get_calc

Calculates additional variables from the model outputs (uses the datasets from `get_domain` and `get_pf`). This was written for a PF-CLM run.  The only difference for PF-only would be to calculate the boundary flux instead of evaptranssum (can do both but evaptranssum = 0) and change the `tflux` calculation function from `calc.pfclm_total_flux` to `calc.pf_total_flux`, using the boundary flux as input instead of evaptranssum.

## pfpostproc.calc

Module includes calculation functions, including water balance.  Many of these calculations were adapted from the fortran/c pftools scripts.  In those cases, I cross-checked my calculations with what the pftools spit out and they matched exactly.  The water balance calculation takes different inputs depending on whether you are running with or without CLM but only looks at the water balance within PF (not CLM) in either case.  The difference is that the boundary flux (in) is explicitly set in the TCL for PF-only, whereas with PF-CLM the evaptranssum variable is used instead (I suppose both could be used but I didn't try doing that).  The problem with PF-only is that the boundary flux does not get output into a variable so it has to be manually set.  It's weird that PF-only outputs evaptranssum but it only contains zeros.  Feels like it ought to contain the boundary flux values.  This gets even more challenging in PF-only when you set up a rain/recession cycle (where it rains for a set amount of time, then stops for a set amount of time).

## pfpostproc.attrs

Metadata dicts for adding attributes to the xarray datasets (contains both domain and variable-specific). Should help explain what the variables are as a reference.

# Water Balance Error, Water Balance Percent Error

From looking at various example scripts, the general procedure is to sum up the storage terms (surface, subsurface) and find the amount that changed from one timestep to the next.  In my case, I get the change in terms of volume of water.  Theoretically, this change should be explained by the net flux volume of water moving in/out of the domain.  So I add up the flux terms ([evaptranssum,overlandsum] for PF-CLM, [boundary flux, overlandsum] for PF-only).  (-) should indicate water leaving the domain and (+) should indicate water entering the domain.  This is tricky with PF-only because the boundary flux set in the TCL script uses (-) for water entering the domain.  

The Water Balance Error is the difference between the change in storage and the net total flux.  This number can be valuable, but since it will never be zero is hard to interpret.  So additional calculations are needed to get a percent error based on the expected storage based on the flux terms.  The expected storage is the total storage at time=t-1 + the total flux at time=t. To get this for the whole time series you need to include the initial storage (which are output separately from the timeseries outputs).