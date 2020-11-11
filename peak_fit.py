import pprint
from collections import OrderedDict
from inspect import signature

import numpy as np
from scipy.optimize import curve_fit

from lblcrn.common.curve import Curve
from lblcrn.common.util import unweave
from lblcrn.xps_data_processing.fitting_suggestions import \
    suggest_fitting_params
from lblcrn.xps_data_processing.line_shapes import line_shapes


class PeakFit:
    def __init__(self,
                 species,
                 line_shape="glp",
                 suggest_params=True,
                 curve=None,
                 ):
        """
        :param species: the name of the species to use in fitting;
        :param suggest_params: suggest parameters if set to True;
        :param line_shape: line_shape to use in the fitting. Default value is Gaussian-Lorentzian Product line shape;
        :param curve: a dataframe whose only column is the curve to fit.
        """
        self.curve = curve
        self.line_shape = line_shape
        # Store mappings from line_shape argument names to a list of [initial guess, lower bound, upper bound].
        # One mapping should exist for each peak. Usually there should be up to two peaks.
        self.kwargs_list = []

        if suggest_params:
            suggestions = suggest_fitting_params(species)

            if suggestions is None or suggestions.empty:
                print(f"No suggestions found for species {species}. Please supply your initial guesses.")

            if not line_shape:
                line_shape = suggestions.iloc[0]["line_shape"]

            line_shape_signature = list(signature(line_shapes[line_shape]["function"]).parameters.keys())

            line_shape_params_size = len(line_shape_signature) - 2 - 2
            suggested_line_shape_params = suggestions["line_shape_params"]

            if line_shape_params_size != len(suggested_line_shape_params[0]):
                print(f"Suggestions database has {len(suggested_line_shape_params[0])} " +
                      f"suggested parameters whereas line shape {line_shape} requires {line_shape_params_size} " +
                      f"number of parameters")

            for peak_number in range(suggestions["be"].size):
                kwargs = OrderedDict()
                center_guess = suggestions["be"][peak_number]
                center_error = suggestions["be_error"][peak_number]
                kwargs["center"] = [center_guess, center_guess - center_error, center_guess + center_error]

                fwhm_guess = suggestions["fwhm"][peak_number]
                fwhm_error = suggestions["fwhm_error"][peak_number]
                kwargs["fwhm"] = [fwhm_guess, fwhm_guess - fwhm_error, fwhm_guess + fwhm_error]

                i = 0
                for param_name in line_shape_signature:
                    if param_name not in ["x", "height", "center", "fwhm"]:
                        kwargs[param_name] = [suggested_line_shape_params[peak_number][i]]
                        i += 1

                self.add_peak(**kwargs)

        # A Curve object to organize and visualize the fitting results.
        self.fitting_result = None

    def add_peak(self, peak_index=-1, **kwargs):
        """
        Add or update an individual peak in the fitting process
        :param peak_index: index of the peak whose initial argument to add to or to overwrite;
                            default value -1 means that this peak is totally new.
        :param kwargs: keyword arguments for passing to the line_shape. All values have to be list.
                       If the list has only one element, the element is interpreted as a strong guess, with an error
                       at 1% of the value;
                       If the list has two numbers, both numbers are used as two sides of a hard constraint on the
                       parameter. Initial guesses are taken as the mean of the two numbers;
                       If the list has three numbers, the first number is taken as the initial guess; both subsequent
                       numbers are used as hard constraints;
        :return: None
        """
        resulting_kwargs = OrderedDict()
        if "height" not in kwargs:
            height_upper_bounds = self.curve.iloc[:, 0].max()
            height_initial_guess = height_upper_bounds
            resulting_kwargs["height"] = [height_initial_guess, 0, height_upper_bounds]

        if peak_index != -1:
            resulting_kwargs.update(self.kwargs_list[peak_index])
        else:
            resulting_kwargs.update(kwargs)
        for k, args in kwargs.items():
            if len(args) == 1:
                args = [args[0], args[0] * 0.8, args[0] * 1.2]
            elif len(args) == 2:
                args = [(args[0] + args[1]) / 2] + args
            resulting_kwargs[k] = args

        if peak_index != -1:
            self.kwargs_list[peak_index] = resulting_kwargs
        else:
            self.kwargs_list.append(resulting_kwargs)

    def fit(self, curve=None):
        """
        Fit the self.curve or input dataframe curve according to fitting parameters in self.kwargs_list according to
        the line shape function identifier self.line_shape.
        :param curve: the curve to use the parameters in this object to fit for;
                      default value is None,
        :return: if curve is provided, a Curve object to visualize fitting results;
                 otherwise, store the resulting Curve object in self.fitting_result and visualize the fitting results.
        """
        f = self._produce_line_shape_function(self.line_shape)
        initial_guesses = np.array([[param_list[0] for param_list in kwargs.values()]
                                    for kwargs in self.kwargs_list]).flatten()

        lower_bound = np.array([[param_list[1] for param_list in kwargs.values()]
                                for kwargs in self.kwargs_list]).flatten()

        upper_bound = np.array([[param_list[2] for param_list in kwargs.values()]
                                for kwargs in self.kwargs_list]).flatten()

        print(f"Fitting {len(self.kwargs_list)} peaks")
        self.show_fitting_params()

        curve_to_fit = self.curve if curve is None else curve
        popt, pcov = curve_fit(f, np.flip(curve_to_fit.index.to_numpy()), np.flip(curve_to_fit.iloc[:, 0].to_numpy()),
                               p0=initial_guesses, bounds=[lower_bound, upper_bound])
        final_params = unweave(popt, line_shapes[self.line_shape]["num_params"] - 1)
        print(f"Fitted parameters {list(zip(*[l.tolist() for l in final_params]))}")

        fitting_result = Curve(f,
                               curve_to_fit.index.min(),
                               curve_to_fit.index.max(),
                               len(curve_to_fit.index),
                               curve_to_fit,
                               *final_params)

        if curve is None:
            self.fitting_result = fitting_result
            self.plot()
        else:
            return fitting_result

    def plot(self):
        """
        Plot the result of the fitting, including individual peaks, combined curve, against the original curve.
        """
        if self.fitting_result:
            return self.fitting_result.plot()
        else:
            print("Fitting not completed.")

    def show_fitting_params(self):
        """
        Show the current fitting parameters.
        """
        pp = pprint.PrettyPrinter(indent=4)

        # Round numerical values in kwargs_list
        rounded_kwargs_list = []
        for kwargs in self.kwargs_list:
            rounded_kwargs = OrderedDict()
            for k, v in kwargs.items():
                rounded_kwargs[k] = [round(p, 3) for p in v]
            rounded_kwargs_list.append(rounded_kwargs)
        pp.pprint(dict(enumerate(rounded_kwargs_list)))

    @staticmethod
    def _produce_line_shape_function(type="glp"):
        """
        :param type: the type of the fitting function;
        :return: a function used for peak fitting.
        """

        func = line_shapes[type]["function"]
        num_params = line_shapes[type]["num_params"] - 1

        def f(xs, *params):
            """
            Function of the assumed peak. This is the model; we optimize the params for this
            function towards certain goals.

            :param xs: numpy array containing x values;
            :param params: a list of additional parameters; every num_params parameters forms the parameters
                           for one individual peak;
            :return: a numpy array produced
            """
            params = unweave(params, num_params)
            ys = np.zeros_like(xs)
            for param_set in zip(*params):
                ys += func(xs, *param_set)
            return ys
        return f
