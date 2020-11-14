import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import argrelextrema, savgol_filter

from lblcrn.xps_data_processing.baseline import shirley_background
from lblcrn.xps_data_processing.peak_fit import PeakFit


class RawRegion:
    def __init__(self, header_block=None, info_block=None, run_mode_block=None, data_block=None):
        """
        :param header_block: string buffer for the header;
        :param info_block: string buffer for the section containing block information;
        :param run_mode_block: string buffer for the section containing run mode;
        :param data_block: string buffer for the section containing raw data.
        """
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
            self.calibration_offset = None
            self.shirley_background = None

            # The PeakFitter for this object and a Curve object to store and visualize the result from the peak_fitter.
            self.peak_fitter = None
            self.peak_fit_result = None

    def plot(self, data=None, **kwargs):
        # TODO: see if kwargs explanation is good.
        """
        A wrapper around Pandas' plotting function, except the following functionalities:
            1. The axis is always inversed;
            2. The title is assigned as current region's name.

        :param data: the data to visualize; if set to None, visualize data in this region itself;
        :param kwargs: any arguments to send into the plotting method of a Pandas dataframe;
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
        Apply Savgol filter to the region. Savgol filter is the standard tool for smoothening noisy signals.
        :return: a copy of the region object with "data" column smoothed, and "1st derivative" and "2nd derivative"
        columns attached to the dataframe.
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

    def fit_peaks(self,
                  peak_fitter=None,
                  automatic_mode=True):
        """
        Fit current region into one or several peaks of certain line shapes.

        :param peak_fitter: a PeakFit object to fit the region with. Default value is None. If a PeakFit object is
                            provided, it will be used regardless of the automatic_mode parameter;
        :param automatic_mode: if set to default value True, the system suggests parameters, and performs the peak
                               fitting;
                               if set to False, return a PeakFit object to allow for manual peak fitting;
        :return: a PeakFitter object used to fit the data stored in this region and organize and visualize the fitting
                 results.
        """
        if peak_fitter:
            self.peak_fitter = peak_fitter
            self.peak_fit_result = peak_fitter.fit(curve=self.data)
        elif automatic_mode:
            peak_fitter = PeakFit(self.name, suggest_params=True, curve=self.data)
            self.peak_fitter = peak_fitter
            self.peak_fit_result = peak_fitter.fit(curve=self.data)
        else:
            peak_fitter = PeakFit(self.name, suggest_params=False, curve=self.data)
            self.peak_fitter = peak_fitter
        return peak_fitter

    def find_steepest_section(self, column="data", extrema_func="extreme"):
        """
        :param column: the column in self.data for which we find the steepest section;
        :param extrema_func: if set to "extreme", take the largest and the smallest data points as extrema;
                             if set to "zero", take the largest and the smallest values closest to 0 as extrema;
        :return: the point where the self.data is the steepest.
        """
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

    def calibrate(self, offset=0):
        """
        Calibrate this region by a numerical value. Practically, this method deducts the binding energy column by an
        offset.

        :param offset: the calibration value.
        :return: None.
        """
        if self.calibration_offset is None:
            self.calibration_offset = offset
            self.data.index = self.data.index - offset
        else:
            print(f"This region has been calibrated with offset={offset}. Consider calling undo_calibration() method" +
                  f" first")

    def undo_calibration(self):
        """
        Undo the current calibration, whose offset value has been stored in self.calibration_offset.
        """
        if self.calibration_offset is None:
            print("Calibration has not been performed. Please calibrate this region first.")
        else:
            self.data.index = self.data.index + self.calibration_offset
            self.calibration_offset = None

    def show_calibration(self, ax=None):
        """
        Visualize the original region before calibration was applied and the calibration offset value.

        :param ax: the Matplotlib axis object to use for this plot.
                   When set to default value None, use Matplotlib's current axis.
        """
        if self.calibration_offset is None:
            print("Calibration has not been performed. Please calibrate this region first.")
        else:
            original_data = self.data.copy()
            original_data.index = original_data.index + self.calibration_offset
            if ax is None:
                ax = plt.gca()
            self.plot(data=original_data, ax=ax)
            ax.axvline(x=self.calibration_offset, ymin=0, ymax=1, linestyle='dashed', lw=1,
                       label=f"offset={round(self.calibration_offset, 2)}")

    def apply_shirley_background(self, max_iters=50):
        """
        Calculate the Shirley Background and subtract it from self.data.

        :param max_iters: maximum number of iterations to use in the Shirley algorithm;
        :return: None
        """
        if self.shirley_background is not None:
            raise Exception("Shirley background has been calculated.")
        self.shirley_background = shirley_background(self.data.copy(), max_iters=max_iters)
        self.data = self.data - self.shirley_background

    def show_shirley_background_calculations(self, ax=None):
        """
        Plot the original data and the shirley background calculated based on it.

        :param ax: the Matplotlib axis object to use for this plot.
                   When set to default value None, use Matplotlib's current axis.
        :return: None
        """
        calculations_df = self.data.copy()
        if self.shirley_background is None:
            raise Exception("Shirley Background has not been calculated. "
                            + "Use apply_shirley_background method first.")
        calculations_df += self.shirley_background
        calculations_df["Shirley Background"] = self.shirley_background
        calculations_df.rename(columns={'data': 'Original Signal'}, inplace=True)
        self.plot(data=calculations_df, ax=ax)

    @property
    def energy_level(self):
        """
        :return: the excitation energy in this region.
        """
        return int(self.info_df[1]["Excitation Energy"])
