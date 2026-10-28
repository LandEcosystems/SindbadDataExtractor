#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
In-situ site values from a consolidated csv with one row per site, such as
soilDepth_Shan2026_consolidated.csv, vegFrac_De2026_consolidated.csv and
treeFrac_De2026_consolidated.csv.

The csv at data_path has a Site_ID column and the column sourceVariableName. Site IDs
are matched without case. The value has no time stamp and is placed on every time
step. Sites missing from the csv get NaN.
"""
import xarray as xr
import numpy as np

import extractors.footprint_utils as fu
from extractors.BasexTractor import BasexTractor
import logging
logger = logging.getLogger(__name__)


def extract(dataset, site_info, config):
    bxtr = BasexTractor(dataset=dataset, site_info=site_info, config=config)

    if not bxtr.is_resolution_supported(provided_reso='hourly'):
        return None

    date_ = bxtr.get_date_vec()

    src_data = []

    for tar_name in bxtr.vars_list:

        bxtr.log_var_start(tar_name)

        path = bxtr.vars[tar_name]['data_path']
        val = fu.site_value(path, bxtr.site, bxtr.vars[tar_name]['sourceVariableName'])
        if not np.isfinite(val):
            logger.warning(f"::MISSING:: variable {tar_name} has no data in source {path} for {bxtr.site}. NaN values will be set.")

        data = xr.Dataset({
            tar_name:
            xr.DataArray(data=np.full(len(date_), val), dims=['time'],
                         coords={'time': date_})
        })
        data = bxtr.convert_units(data, tar_name)
        bxtr.log_var_end(data, tar_name, None)
        src_data.append(data)

    src_dataset = bxtr.merge_and_format(src_data)
    return src_dataset


if __name__ == '__main__':
    import inspect, os
    print('---------------------------------------------------')
    print(
        'Provider at: ',
        os.path.dirname(
            os.path.abspath(inspect.getfile(inspect.currentframe()))))
    print('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)
