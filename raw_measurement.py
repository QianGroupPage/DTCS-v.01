import io

import pandas as pd

from lblcrn.xps_data_processing.raw_region import RawRegion
from lblcrn.xps_data_processing.utilities import read_line_blocks


class RawMeasurement:
    def __init__(self, file_path, sequence_number):
        self.regions_list = []
        self.regions = {}
        self.sequence_number = sequence_number
        splitted_lines = read_line_blocks(file_path)
        header_df = pd.read_csv(io.StringIO(''.join(splitted_lines[0])), sep="=")
        num_regions = int(header_df["[Info]"]['Number of Regions'])
        for i in range(num_regions):
            region_header_block = io.StringIO(''.join(splitted_lines[4 * i + 1][1:]))
            region_info_block = io.StringIO(''.join(splitted_lines[4 * i + 2][1:]))
            # TODO: get examples where region_run_mode has something in it, and modify
            # this block accordingly.
            region_run_mode_block = io.StringIO(''.join(splitted_lines[4 * i + 3]))
            region_data_block = io.StringIO(''.join(splitted_lines[4 * i + 4][1:]))
            region = RawRegion(region_header_block, region_info_block, region_run_mode_block, region_data_block)
            self.regions_list.append(region)
            self.regions[region.name] = region

    def combine_regions(self):
        """
        Combine different regions into one common object.
        """
        pass

    def list_region_names(self):
        return [region.name for region in self]

    def calibrate(self, internal_region="VB"):
        """
        :param internal_region: region name or the beginning part of the region to use for alignment.
        :return:
        """
        self.calibrate_by_numerical_value(self.calibration_offset(internal_region=internal_region))

    def has_internal_region(self, internal_region="VB"):
        """
        Return True if the measurement has an internal region that matches with name internal_region.
        """
        for region_name in self.list_region_names():
            if region_name.startswith(internal_region):
                return True
        return False

    def calibration_offset(self, internal_region="VB"):
        region_name_in_use = None
        for region_name in self.list_region_names():
            if region_name.startswith(internal_region):
                if region_name_in_use:
                    raise ValueError(f"Both {region_name_in_use} and {region_name_in_use} start with {internal_region}")
                region_name_in_use = region_name
        return self.regions[region_name_in_use].produce_smoothed_version().find_steepest_section()

    def calibrate_by_numerical_value(self, value):
        """
        Add the index (usually binding energy) of each region with value.
        """
        for region in self:
            region.data.index = region.data.index - value

    def remove_baseline(self, max_iters=50):
        """
        Calculate and subtract Shirley background from each region.
        """
        for region in self:
            region.apply_shirley_background(max_iters=max_iters)

    @property
    def num_regions(self):
        return len(self.regions)

    @property
    def energy_levels(self):
        return list(set([region.energy_level for region in self]))

    def __getitem__(self, i):
        return self.regions_list[i]

