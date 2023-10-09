#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    @author Kylene Cooley
    @brief Data editing and visualization tools to use during data deep
           dives.
"""

import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ooi_data_explorations.common import load_kdata
from ooi_data_explorations.uncabled.process_flort import flort_datalogger
from ooi_data_explorations.uncabled.process_metbk import metbk_datalogger

def nanfill_time_gaps(dataset, freq='3H'):
    """ Use this function to create sections with nans in time series
    at times where no data exists in the record. After preprocessing,
    the time coordinate does not include times where no data exists.
    Adding data variables with nans in these spaces allows user to
    graph time series in line plots with gaps in the time series,
    rather than having a straight line connect data across these gaps.

    Input:
    -------
    :param dataset: xarray Dataset without nans in data variables or
                    coordinates
    :param freq: datetime-like time interval between observations,
                 match to original dataset
                 
    Returns:
    --------
    :return dataset_full: xarray Dataset with same coordinates and
                          variables as dataset but expanded in the time
                          dimension with nans where no data was
                          recorded
    """
    # Set start and end date times from dataset input
    startDT = dataset.time[0]
    endDT = dataset.time[-1]

    # Create regularly spaced time coordinate for the whole time series record
    complete_time = pd.date_range(start=startDT.values, end=endDT.values,
                                  freq=freq)
    variable_array = np.full_like(complete_time, np.nan, dtype=np.float64)

    # Create dataset to combine with original dataset
    nan_ds = xr.Dataset(data_vars=dict(variable=(["time"], variable_array)),
                        coords=dict(time=complete_time))

    # Fill in missing times with nans by matching up coordinates and at
    # least one variable
    dataset_full = xr.combine_by_coords([dataset, nan_ds])

    # Drop working nan variable after datasets are combined
    dataset_full = dataset_full.drop_vars('variable')
    
    return dataset_full


def check_chla_swr(spkir, site, deploy, flort_node):
    """Plot downwelling spectral irradiance in comparison with SWR from
    the surface buoy METBK suite and Chlorophyll-a from the co-located
    FLORT. This calls for the site where the SPKIR is located.
    which will also be the same for METBK and FLORT sensors. While the
    FLORT is also on the NSIF, its data goes through a different DCL such
    that the node for the FLORT is different.
    """
    # Load METBK and FLORT data
    met_node = 'SBD11' # not all subsites have a second METBK (SBD12)
    met_sensor = '06-METBKA000'
    met_method = 'recovered_host'
    met_stream = 'metbk_a_dcl_instrument_recovered'
    metbk = load_kdata(site, met_node, met_sensor, met_method, met_stream,
                       ('*deployment%04d*METBK*.nc' % deploy))
    metbk = metbk_datalogger(metbk)

    flort_sensor = '02-FLORTD000'
    flort_method = 'recovered_host'
    flort_stream = 'flort_sample'
    flort = load_kdata(site, flort_node, flort_sensor, flort_method,
                       flort_stream, ('*deployment%04d*FLORT*.nc' % deploy))
    flort = flort_datalogger(flort)

    # Create subplots
    fig, ax = plt.subplots(3,1, sharex=True, figsize=(15,8))
    for var in spkir.variables:
        if "downwelling_irradiance" in var:
            spkir[var].plot(ax=ax[0], label=spkir[var].radiation_wavelength)
    ax[0].set_ylabel(
        'Downwelling Spectral \n Irradiance \n [uW cm$^{-2}$ nm$^{-1}$]')
    ax[0].legend()

    flort['estimated_chlorophyll'].where(
        flort['estimated_chlorophyll_qc_summary_flag'] != 4).plot(ax=ax[1])

    metbk['shortwave_irradiance'].where(
        metbk['shortwave_irradiance_qc_summary_flag'] != 4).plot(ax=ax[2])
    
    plt.show()
    return metbk, flort, fig, ax