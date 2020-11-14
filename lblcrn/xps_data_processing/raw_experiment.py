import os
import re

import matplotlib.pyplot as plt
import pandas as pd
from IPython.display import display

from lblcrn.xps_data_processing.peak_fit import PeakFit
from lblcrn.xps_data_processing.raw_measurement import RawMeasurement
from lblcrn.xps_data_processing.read_digital_book import read_digital_notebook


class RawExperiment:
    def __init__(self, directory_path, automatic_processing=True):
        """
        Process a set of XPS experiment results into a Python class;
        :param directory_path: the path where the experiment is stored;
        :param automatic_processing: when set to default value True, calibrate, calculate Shirley Background, and
                                     perform peak fitting in automatic mode;
                                     otherwise, don't do any of the tasks.
        """
        measurements = []
        for file_path in os.listdir(directory_path):
            full_file_path = f"{directory_path}/{file_path}"

            if self.is_digital_notebook(file_path):
                self.general_comments_line_blocks, self.conditions = \
                    read_digital_notebook(full_file_path)

            if self.is_measurement_file(file_path):
                measurement = RawMeasurement(full_file_path)
                measurements.append(measurement)
        measurements.sort(key=lambda s: s.sequence_number)
        self.measurements = measurements

        for condition in self.conditions:
            included_measurements = []
            for s in condition.measurement_sequence_numbers:
                included_measurements.append(measurements[s - 1])
            condition.set_measurements(included_measurements)

        if automatic_processing:
            self.show_all_conditions()

            print("Start calibration.")
            self.calibrate()
            self.show_calibration_results()

            print("Start Shirley Background calculations.")
            self.remove_baseline()
            self.show_baselines()

            print("Start peak fitting.")
            self.peak_fit()
            self.show_fitting_results()

    def experimental_history(self):
        pass

    def show_all_conditions(self):
        """
        Display all conditions in the experiment as a table.
        :return: None
        """
        column_names_list = self._all_species()

        df = pd.DataFrame(columns=["Condition Identifier"] + column_names_list)
        for i, condition in enumerate(self.conditions):
            df.loc[i] = [condition.id] + [column_name if column_name in condition.list_region_names() else " "
                               for column_name in column_names_list]
        df.set_index("Condition Identifier", inplace=True)

        pd.set_option("display.max_rows", None, "display.max_columns", None)
        display(df)

    def get_measurement(self, measurement_id=None):
        pass

    def calibrate(self, internal_region="VB"):
        """
        Calibrate all measurements in this experiment.
        :param internal_region: region name or the beginning part of the region to use for alignment.
        :return: None
        """
        for i, measurement in enumerate(self.measurements):
            # TODO: quality check.
            if measurement.has_internal_region(internal_region):
                measurement.calibrate(internal_region=internal_region)
            else:
                reference_measurement = self.measurements[self.find_previous_measurement(i)]
                print(f"Region starting with {internal_region} is not present in measurement " +
                      f"{measurement.sequence_number}." +
                      f" Using measurement {reference_measurement.sequence_number} for calibration.")
                measurement.calibrate_by_numerical_value(
                    reference_measurement.calibration_offset(internal_region=internal_region))

    def show_calibration_results(self):
        """
        Visualize the calibration results in a grid.
        """
        cols = len(self._all_species())
        fig, axes = plt.subplots(len(self.conditions), cols)

        for i, condition in enumerate(self.conditions):
            for measurement in condition.measurements:
                for region in measurement:
                    region.show_calibration(ax=axes[i, self._all_species().index(region.name)])

    def remove_baseline(self, max_iters=50):
        """
        Calculate and subtract Shirley background from each measurement.

        :param max_iters: maximum number of iterations for use in Shirley Background calculations.
        :return: None
        """
        for measurement in self.measurements:
            measurement.remove_baseline(max_iters=max_iters)

    def show_baselines(self):
        """
        Visualize all baseline calculations in this RawExperiment.
        :return:
        """
        cols = len(self._all_species())
        fig, axes = plt.subplots(len(self.conditions), cols)

        for i, condition in enumerate(self.conditions):
            for measurement in condition.measurements:
                for region in measurement:
                    region.show_shirley_background_calculations(axes[i, self._all_species().index(region.name)])

    def find_previous_measurement(self, measurement_number, region_name=""):
        """
        :param measurement_number: the measurement for which we're looking for previous energy.
        :param region_name: the name
        :return: measurement number for a previous measurement with the same energy level.
        """
        i = measurement_number - 1
        while i >= 0:
            if self.measurements[i].has_internal_region(region_name):
                current_energy_levels = self.measurements[measurement_number].energy_levels
                measurement_i_energy_levels = self.measurements[i].energy_levels

                if sorted(current_energy_levels) == sorted(measurement_i_energy_levels):
                    return i
            i -= 1
        return -1

    def peak_fit(self,
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
                If name_to_peak_fitter is not provided and automatic_mode is set to True, the peak fitter objects come
                with suggested initial parameters,
                Otherwise call each peak_fitter object's add_peak() method to add a suggested peak to the object, call
                this method again with name_to_peak_fitter set to this returned dictionary of peak fitter objects.
        """
        ignored_by_default = ["survey", "VB"]
        ignored_species.extend(ignored_by_default)

        if species is None:
            species = self._all_species()
        species = [s for s in species if s not in ignored_species]

        name_to_peak_fitter = {s: PeakFit(s, suggest_params=automatic_mode) if s not in name_to_peak_fitter else
                               name_to_peak_fitter[s] for s in species}

        for measurement in self.measurements:
            measurement.fit_peaks(species=measurement.list_region_names(),
                                  ignored_species=ignored_species,
                                  name_to_peak_fitter=name_to_peak_fitter,
                                  automatic_mode=automatic_mode)

        # TODO: Judge which peak_fitter to run by the beginning of the names.
        return name_to_peak_fitter

    def show_fitting_results(self):
        """
        Visualize the results from peak fitting in a grid.
        """
        cols = len(self._all_species())
        fig, axes = plt.subplots(len(self.conditions), cols)

        for i, condition in enumerate(self.conditions):
            for measurement in condition.measurements:
                for region in measurement:
                    region.peak_fitter.plot(axes[i, self._all_species().index(region.name)])

    @staticmethod
    def is_measurement_file(filename):
        """
        :param filename: any filename;
        :return: True if the file ends with "_" following 4 digits.
        """
        return re.match(r'\w*_\d\d\d\d.txt$', filename)

    @staticmethod
    def is_digital_notebook(filename):
        """
        Enforce the following format for the digital notebook:
        'digital{any separator}notebook' must be part of the file name. Word cases don't matter.

        :param filename: any filename;
        :return: True if filename follows the format rule.
        """
        allowed_separators = [" ", "_"]
        for separator in allowed_separators:
            filename.replace(separator, " ")

        splitted_lower_filename = filename.lower().split()
        for i, w in enumerate(splitted_lower_filename):
            if w == "digital" and i < len(splitted_lower_filename) and splitted_lower_filename[i + 1] == "notebook":
                return True
        return False

    def _all_species(self):
        """
        :return: a list with all species in the experiment.
        """
        region_names_list = []
        for condition in self.conditions:
            for name in condition.list_region_names():
                if name not in region_names_list:
                    region_names_list.append(name)
        return region_names_list

    def __getitem__(self, name):
        return self.measurements[name]
