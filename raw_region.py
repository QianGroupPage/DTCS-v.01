import numpy as np
import pandas as pd
from scipy.signal import argrelextrema, savgol_filter

from lblcrn.xps_data_processing.baseline import shirley_background


class RawRegion:
    def __init__(self, header_block=None, info_block=None, run_mode_block=None, data_block=None):
        if header_block and info_block and run_mode_block and data_block:
            region_header_df = pd.read_csv(header_block, sep="=", header=None, index_col=0)
            region_info_df = pd.read_csv(info_block, sep="=", header=None, index_col=0)
            # region_run_mode_df = pd.read_csv(run_mode_block, sep="=")
            region_data_df = pd.read_csv(data_block, delim_whitespace=True, header=None)

            self.info_df = region_info_df
            self.header_df = region_header_df

            # TODO: figure out whether to use the scale column in header. It apparently is supposed to completely match
            # first column in data_df.
            region_header_df[1]["Dimension 1 scale"] = \
                [float(word) for word in region_header_df[1]["Dimension 1 scale"].split()]
            # region_data_df[region_header_df[1]["Dimension 1 name"]] = region_header_df[1]["Dimension 1 scale"]
            # region_data_df.set_index(region_header_df[1]["Dimension 1 name"], inplace=True)

            region_data_df.set_index(0, inplace=True)
            region_data_df.rename_axis(region_header_df[1]["Dimension 1 name"], inplace=True)
            region_data_df.rename(columns={1: "data"}, inplace=True)

            self.name = region_header_df[1]["Region Name"]
            self.data = region_data_df
            self.parent = None

            # Placeholders for signal processing results
            self.shirley_background = None

    def plot(self, data=None, **kwargs):
        """
        A wrapper around Pandas' plotting function, except the following functionalities:
            1. The axis is always inversed;
            2. The title is assigned as current region's name.
        """
        if data is None:
            data = self.data
        if "xlim" in kwargs:
            kwargs["xlim"] = kwargs["xlim"].sorted()
        else:
            kwargs["xlim"] = data.index.max(), data.index.min()
        if "title" not in kwargs:
            kwargs["title"] = self.name
        data.plot(**kwargs)

    def produce_smoothed_version(self):
        """
        Savgol filter: the standard tool for smoothening noisy signals.
        """
        new_version = RawRegion()
        new_version.name = self.name
        new_version_data = self.data.copy()

        # print(self.data.size)
        # print(len(savgol_filter(self.data["data"].to_numpy(), 61, 3, mode="nearest")))
        new_version_data["data"] = savgol_filter(self.data["data"].to_numpy(), 61, 3, mode="nearest")

        new_version_data["1st derivative"] = savgol_filter(self.data["data"].to_numpy(), 61, 3, deriv=1, mode="nearest")
        new_version_data["2nd derivative"] = savgol_filter(self.data["data"].to_numpy(), 61, 3, deriv=2, mode="nearest")

        new_version.data = new_version_data
        return new_version

    def find_steepest_section(self, column="data", extrema_func="extreme"):
        if extrema_func == "extreme":
            extremas = sorted(argrelextrema(self.data[column].to_numpy(), np.greater)[0].tolist() + (
                argrelextrema(self.data[column].to_numpy(), np.less)[0].tolist()))
        elif extrema_func == "zero":
            # The following should take an array
            extremas = argrelextrema(self.data[column].to_numpy(),
                                     lambda x1, x2: np.less(np.abs(x1), np.abs(x2)))[0].tolist()
        sections = [(self.data.index[a], self.data.index[b]) for a, b in zip(extremas[:-1], extremas[1:])]
        best_a, best_b = sections[0]
        max_deriv = np.abs((self.data["data"][best_a] - self.data["data"][best_b]) / (best_a - best_b))
        for a, b in sections:
            deriv = np.abs((self.data["data"][a] - self.data["data"][b]) / (a - b))
            if deriv > max_deriv:
                best_a, best_b = a, b
        return (best_a + best_b) / 2

    def apply_shirley_background(self, max_iters=50):
        """
        Calculate Shirley Background and subtract it from self.data.
        """
        if self.shirley_background is not None:
            raise Exception("Shirley background has been calculated.")
        self.shirley_background = shirley_background(self.data.copy(), max_iters=max_iters)
        self.data = self.data - self.shirley_background

    def show_shirley_background_calculations(self):
        """
        Show the original data and the shirley background calculated based on it.
        """
        calculations_df = self.data.copy()
        if self.shirley_background is None:
            raise Exception("Shirley Background has not been calculated. "
                            + "Use apply_shirley_background method first.")
        calculations_df += self.shirley_background
        calculations_df["Shirley Background"] = self.shirley_background
        calculations_df.rename(columns={'data': 'Original Signal'}, inplace=True)
        self.plot(data=calculations_df)

    @property
    def energy_level(self):
        """
        Find the excitation energy in this region.
        :return:
        """
        return int(self.info_df[1]["Excitation Energy"])
