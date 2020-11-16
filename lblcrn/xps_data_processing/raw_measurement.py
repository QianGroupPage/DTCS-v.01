import io
import os

import pandas as pd

from lblcrn.xps_data_processing.raw_region import RawRegion
from lblcrn.xps_data_processing.utilities import read_line_blocks


class RawMeasurement:
    def __init__(self, file_path):
        """
        :param file_path: path of the machine-generated measurement file.
        """
        self.regions_list = []
        self.regions = {}
        self.sequence_number = int(os.path.splitext(file_path)[0].split('_')[-1]) - 1
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
        """
        :return: a list containing the names of all measured regions in the file.
        """
        return [region.name for region in self]

    def calibrate(self, internal_region="VB"):
        """
        Calibrate by an internal region.

        :param internal_region: name or prefix of the name of the region to use for alignment;
        :return: None
        """
        self.calibrate_by_numerical_value(self.calibration_offset(internal_region=internal_region))

    def calibrate_by_numerical_value(self, offset):
        """
        :param offset: the offset value to add to the binding energy column;
        :return: None
        """
        for region in self:
            region.calibrate(offset=offset)

    def has_internal_region(self, internal_region="VB"):
        """
        :param internal_region: name or prefix of the name of the region to use for alignment;
        :return: True if the measurement has an internal region that matches with name internal_region.
        """
        for region_name in self.list_region_names():
            if region_name.startswith(internal_region):
                return True
        return False

    def calibration_offset(self, internal_region="VB"):
        """
        Calibrate all regions in this measurement by the region whose name starts with "VB".

        :param internal_region: name or any prefix of the region by which this calibration is based on;
        :return: numerical value of the offset.
        """
        region_name_in_use = None
        for region_name in self.list_region_names():
            if region_name.startswith(internal_region):
                if region_name_in_use:
                    raise ValueError(f"Both {region_name_in_use} and {region_name_in_use} start with {internal_region}")
                region_name_in_use = region_name
        return self.regions[region_name_in_use].produce_smoothed_version().find_steepest_section()

    def remove_baseline(self, max_iters=50, ignored_species=None):
        """
        Calculate and subtract Shirley background from each region.

        :param max_iters: maximum number of iterations to use in the Shirley algorithm;
        :param ignored_species: if set to None, ignore any species starting with "Survey" or "VB";
                                otherwise, ignore all species in the list.
        :return: None
        """
        for region in self:
            if ignored_species is None or region.name not in ignored_species:
                region.apply_shirley_background(max_iters=max_iters)

    def fit_peaks(self,
                  species=None,
                  ignored_species=None,
                  name_to_peak_fitter=None,
                  automatic_mode=True):
        """
        Fit each region in this measurement into one or several peaks of certain line shapes.

        :param species: the species or the list of species this fitter corresponds to;
        :param ignored_species: the species or the list of species this fitter corresponds to;
                                any species starting with "survey" and "VB" are automatically ignored, unless in the
                                species list;
        :param name_to_peak_fitter: a mapping from region name to a peak_fitter. If a PeakFit object is
                                    provided for a region, it will be used regardless of the automatic_mode parameter;
                                    otherwise, the peak_fitter will depend upon the automatic_mode parameter;
        :param automatic_mode: if set to default value True, the system suggests parameters, and performs the peak
                               fitting;
                               if set to False, return a PeakFit object to allow for manual peak fitting;
        :return a dictionary from region name to its corresponding peak fitter.
                The peak fitter objects come with initial suggested parameters, if name_to_peak_fitter is not provided
                and automatic_mode is set to True;
                Otherwise call each peak_fitter object's add_peak() method to add a suggested peak to the object, call
                this method again with name_to_peak_fitter set to this dictionary of peak fitters.
        """
        ignored_by_default = ["survey", "VB"]
        if species is None:
            species = self.list_region_names()
        species = [s for s in species if s not in ignored_species + ignored_by_default]

        if name_to_peak_fitter is None:
            name_to_peak_fitter = {}

            for region in self:
                if region.name in species:
                    name_to_peak_fitter[region.name] = region.fit_peaks(automatic_mode=automatic_mode)
        else:
            for region in self:
                if region.name in species:
                    if region.name in name_to_peak_fitter:
                        name_to_peak_fitter[region.name] = region.fit_peaks(peak_fitter=
                                                                            name_to_peak_fitter[region.name])
                    else:
                        name_to_peak_fitter[region.name] = region.fit_peaks(automatic_mode=automatic_mode)
        return name_to_peak_fitter

    @property
    def num_regions(self):
        """
        :return: number of regions in the measurement.
        """
        return len(self.regions)

    @property
    def energy_levels(self):
        """
        :return: list of energy levels in the measurement.
        """
        return list(set([region.energy_level for region in self]))

    def __getitem__(self, i):
        return self.regions_list[i]

