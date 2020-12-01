from lblcrn.common.line_shapes.phonopy_ir import CalculateIRIntensity, ParseOUTCAR
import re
import numpy as np


def ir_intensity(outcar_file):
    """
    Return the IR intensity based on an outcarfile.

    :param outcarfile: the path for an outcarfile from which the ir intensity is computed;
    :return: the ir intensity.
    """
    # intensity_input_dict = ParseOUTCAR(outcarfile, ['born_charges', 'phonon_modes'])
    # return CalculateIRIntensity(eigendisplacement=intensity_input_dict['phonon_modes'],
    #                             becTensors=intensity_input_dict["born_charges"])
    bec_tensor = parse_born_charges(outcar_file)
    print(bec_tensor)
    # return CalculateIRIntensity(eigendisplacement=intensity_input_dict['phonon_modes'],
    #                             becTensors=bec_tensor)


def parse_born_charges(outcar_file):
    # A list of rows of tensors; every 3 row should correspond to an atom.
    tensor_rows = []
    with open(outcar_file, "r") as file:
        for line in file:
            if re.search("BORN EFFECTIVE CHARGE FOR ION", line):
                floats = [float(s) for s in line.split()[-3:]]
                tensor_rows.append(floats)
    bec_tensor = []
    for i in range(0, len(tensor_rows), 3):
        bec_tensor.append(np.array(tensor_rows[i: i+3], dtype=np.float64))
    return bec_tensor



