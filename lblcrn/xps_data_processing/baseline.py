import numpy as np


def shirley_background(signal, max_iters=50):
    """
    Based on an Igor macro developed by James Mudd at Warick.
    https://www.jamesmudd.com/downloads, accessed Oct 3, 2020.

    :param signal: the input curve, with an index that represents binding energy in decreasing order with perfectly even
    steps;
    :param max_iters: maximum number of iterations to use in the Shirley algorithm;
    :return: the background as a dataframe, where 0-th column is the background.
    """
    # Assuming binding energy decreasing: stop_index < start_index, delta_index < 0
    start_index, stop_index = min(range(0, len(signal.index)), key=lambda i: signal.index[i]), \
                              max((0, len(signal.index) - 1), key=lambda i: signal.index[i])
    delta_index = abs(signal.index[1] - signal.index[0])
    signal_copy = signal.copy() - signal.min()

    background = signal.copy()
    background.iloc[:, 0] = np.zeros_like(signal.index)
    k = 0.001

    iters = 0
    if signal_copy.iloc[start_index, 0] > signal_copy.iloc[stop_index, 0]:
        while iters < max_iters:
            bg_temp = background.copy()
            # Loop through the energies i is in energy
            for curve_iloc_i in range(stop_index, start_index + 1):
                integral = 0
                for curve_iloc_j in range(curve_iloc_i, stop_index - 1, -1):
                    # Shirley Background
                    integral += (signal_copy.iloc[curve_iloc_j, 0] - bg_temp.iloc[curve_iloc_j, 0]) * delta_index
                # CHANGE: WE ALWAYS KNOW THE EXACT VALUE
                background.iloc[curve_iloc_i, 0] = k * integral

            # Adjust k to converge fitting
            background_offset = background.iloc[start_index, 0] - signal_copy.iloc[start_index, 0]
            k = k - (background_offset / background.iloc[start_index, 0])*k*0.25 # 0.25 is to ensure convergence
            background_offset_old = background_offset

            print("start_value", signal_copy.iloc[start_index, 0], "end_value", signal_copy.iloc[stop_index, 0])

            if abs(background_offset) <= 0.000001 * background.iloc[start_index, 0]:
                break
            iters += 1
    else:
        # Switch end points and therefore the direction of the integral.
        start_index, stop_index = stop_index, start_index

        while iters < max_iters:
            bg_temp = background.copy()
            # Loop through the energies i is in energy
            for curve_iloc_i in inclusive_range(stop_index, start_index):
                integral = 0
                for curve_iloc_j in inclusive_range(curve_iloc_i, stop_index):
                    # Shirley Background
                    integral += (signal_copy.iloc[curve_iloc_j, 0] - bg_temp.iloc[curve_iloc_j, 0]) * delta_index
                # CHANGE: WE ALWAYS KNOW THE EXACT VALUE
                background.iloc[curve_iloc_i, 0] = k * integral

            # Adjust k to converge fitting
            background_offset = background.iloc[start_index, 0] - signal_copy.iloc[start_index, 0]
            k = k - (background_offset / background.iloc[start_index, 0]) * k * 0.25  # 0.25 is to ensure convergence
            background_offset_old = background_offset

            if abs(background_offset) <= 0.000001 * background.iloc[start_index, 0]:
                break
            iters += 1

    if iters == max_iters:
        print(f"Shirley algorithm did not converge at allowed maximum number of iterations max_iters={max_iters}.")

    background += signal.min()
    return background


def inclusive_range(x, y):
    """
    :param x: any integer;
    :param y: any integer;
    :return: [x, x - 1, ..., y] if x >= y or [x, x + 1, ..., y] if x < y
    """
    if x < y:
        return range(x, y + 1)
    else:
        return range(y, x + 1)[::-1]
