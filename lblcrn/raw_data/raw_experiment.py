import os
import re
import pandas as pd
from lblcrn.raw_data.raw_sample import RawSample
from IPython.display import display


class RawExperiment:
    def __init__(self, directory_path):
        samples = []
        for file_path in os.listdir(directory_path):
            if self.is_sample_file(file_path):
                sequence_number = int(os.path.splitext(file_path)[0][-4:])
                sample = RawSample(f"{directory_path}/{file_path}", sequence_number)
                samples.append(sample)
        samples.sort(key=lambda s: s.sequence_number)
        self.samples = samples

    def experimental_history(self):
        pass

    def show_all_samples(self):
        column_names_list = []
        for sample in self.samples:
            for name in sample.list_region_names():
                if name not in column_names_list:
                    column_names_list.append(name)

        df = pd.DataFrame(columns=["Sample Identifier"] + column_names_list)
        for i, sample in enumerate(self.samples):
            df.loc[i] = [i] + [column_name if column_name in sample.list_region_names() else " "
                               for column_name in column_names_list]
        df.set_index("Sample Identifier", inplace=True)

        pd.set_option("display.max_rows", None, "display.max_columns", None)
        display(df)

    def get_sample(self, sample_id=None):
        pass

    def calibrate(self, internal_region="VB"):
        for sample in self.samples:
            # TODO: handle the case when internal region is not there.
            # TODO: quality check.
            if sample.has_internal_region(internal_region):
                sample.calibrate(internal_region=internal_region)
            else:
                print(f"Region {internal_region} is not present in sample {sample.sequence_number}")

    @staticmethod
    def is_sample_file(filename):
        """
        A sample file must end with "_" following 4 digits.
        """
        return re.match(r'\w*_\d\d\d\d.txt$', filename)

    def __getitem__(self, name):
        return self.samples[name]
