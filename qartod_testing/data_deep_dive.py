#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    @author Kylene Cooley
    @brief facilitates plotting by manipulating OOI data sets
"""

import xarray as xr
import numpy as np
import pandas as pd


def nanfill_time_gaps(dataset, freq='3H'):
    """
        Use this function to create sections with nans in time series at times where no data exists in the record.
        After preprocessing, the time coordinate does not include times where no data exists. Adding data variables with nans 
        in these spaces allows user to graph time series in line plots with gaps in the time series, rather than having a 
        straight line connect data across these gaps.

        :param dataset: xarray Dataset without nans in data variables or coordinates
        :param freq: datetime-like time interval between observations, match to original dataset
        :return dataset_full: xarray Dataset with same coordinates and variables as dataset but expanded in the time dimension with nans where no data was recorded
    """
    # Set start and end date times from dataset input
    startDT = dataset.time[0]
    endDT = dataset.time[-1]

    # Create regularly spaced time coordinate for the whole time series record
    complete_time = pd.date_range(start=startDT.values, end=endDT.values, freq=freq)
    variable_array = np.full_like(complete_time, np.nan, dtype=np.float64)

    # Create dataset to combine with original dataset
    nan_ds = xr.Dataset(data_vars=dict(variable=(["time"], variable_array)), coords=dict(time=complete_time))

    # Fill in missing times with nans by matching up coordinates and at least one variable
    dataset_full = xr.combine_by_coords([dataset, nan_ds])

    # Drop working nan variable after datasets are combined
    dataset_full = dataset_full.drop_vars('variable')
    
    return dataset_full