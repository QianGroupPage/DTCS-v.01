import pandas as pd
import io
from raw_region import RawRegion


class RawSample:
    def __init__(self, file_path, sequence_number):
        self.regions_list = []
        self.regions = {}
        self.sequence_number = sequence_number
        with open(file_path, "r+") as f:
            lines = f.readlines()

            # Split into lists of lines by newline characters
            size = len(lines)
            idx_list = [idx for idx, val in
                        enumerate(lines) if val == "\n"]

            splitted_lines = [lines[i: j] for i, j in
                              zip([0] + [idx + 1 for idx in idx_list], idx_list +
                                  ([size] if idx_list[-1] != size else [])) if i != j]

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

    @property
    def num_regions(self):
        return len(self.regions)

    def __getitem__(self, i):
        return self.regions_list[i]

