from lblcrn.common.line_shapes.phonopy_ir import CalculateIRIntensity, ParseOUTCAR


def ir_intensity(outcarfile):
    """
    Return the IR intensity based on an outcarfile.

    :param outcarfile: the path for an outcarfile from which the ir intensity is computed;
    :return: the ir intensity.
    """
    # intensity_input_dict = ParseOUTCAR(outcarfile, ['born_charges', 'phonon_modes'])
    # return CalculateIRIntensity(eigendisplacement=intensity_input_dict['phonon_modes'],
    #                             becTensors=intensity_input_dict["born_charges"])
    intensity_input_dict = ParseOUTCAR(outcarfile, ['phonon_modes'])
    return CalculateIRIntensity(eigendisplacement=intensity_input_dict['phonon_modes'])

