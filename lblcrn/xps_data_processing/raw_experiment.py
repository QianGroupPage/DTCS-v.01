import os
import re

import pandas as pd
from IPython.display import display

from lblcrn.xps_data_processing.raw_measurement import RawMeasurement
from lblcrn.xps_data_processing.read_digital_book import read_digital_notebook


class RawExperiment:
    def __init__(self, directory_path):
        measurements = []
        for file_path in os.listdir(directory_path):
            full_file_path = f"{directory_path}/{file_path}"

            if self.is_digital_notebook(file_path):
                self.general_comments_line_blocks, self.conditions = \
                    read_digital_notebook(full_file_path)

            if self.is_measurement_file(file_path):
                sequence_number = int(os.path.splitext(file_path)[0][-4:]) - 1
                measurement = RawMeasurement(full_file_path, sequence_number)
                measurements.append(measurement)
        measurements.sort(key=lambda s: s.sequence_number)
        self.measurements = measurements

        for condition in self.conditions:
            included_measurements = []
            for s in condition.measurement_sequence_numbers:
                included_measurements.append(measurements[s - 1])
            condition.set_measurements(included_measurements)

    def experimental_history(self):
        pass

    def show_all_conditions(self):
        column_names_list = []
        for condition in self.conditions:
            for name in condition.list_region_names():
                if name not in column_names_list:
                    column_names_list.append(name)

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
        for i, measurement in enumerate(self.measurements):
            # TODO: quality check.
            if measurement.has_internal_region(internal_region):
                measurement.calibrate(internal_region=internal_region)
            else:
                reference_measurement = self.measurements[self.find_previous_measurement(i)]
                print(f"Region starting with {internal_region} is not present in measurement {measurement.sequence_number}." +
                      f" Using measurement {reference_measurement.sequence_number} for calibration.")
                measurement.calibrate_by_numerical_value(
                    reference_measurement.calibration_offset(internal_region=internal_region))

    def remove_baseline(self, max_iters=50):
        """
        Calculate and subtract Shirley background from each measurement.
        """
        for measurement in self.measurements:
            measurement.remove_baseline(max_iters=max_iters)

    def find_previous_measurement(self, measurement_number, region_name=""):
        """
        Find a previous measurement at the same energy level.
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


    @staticmethod
    def is_measurement_file(filename):
        """
        A measurement file must end with "_" following 4 digits.
        """
        return re.match(r'\w*_\d\d\d\d.txt$', filename)

    @staticmethod
    def is_digital_notebook(filename):
        """
        Enforce the following format for the digital notebook:
        'digital{any separator}notebook' must be part of the file name. Word cases don't matter.
        """
        allowed_separators = [" ", "_"]
        for separator in allowed_separators:
            filename.replace(separator, " ")

        splitted_lower_filename = filename.lower().split()
        for i, w in enumerate(splitted_lower_filename):
            if w == "digital" and i < len(splitted_lower_filename) and splitted_lower_filename[i + 1] == "notebook":
                return True
        return False

    def __getitem__(self, name):
        return self.measurements[name]
