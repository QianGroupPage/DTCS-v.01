# dtcs.SpectroscoPy/Interfaces/Phonopy.py


# ---------
# Docstring
# ---------

""" Routines for interfacing with the Phono(3)py packages. """


# -------
# Imports
# -------

import math;

import h5py;

import yaml;

# This code tries to use the C version of the YAML Loader and falls back to the standard Loader if this is not available.
# This code was taken from https://pyyaml.org/wiki/PyYAMLDocumentation.

try:
    from yaml import CLoader as Loader, CDumper as Dumper;
except ImportError:
    from yaml import Loader, Dumper;

import numpy as np;

# Needed for parsing Phonopy BORN files.

from phonopy.file_IO import parse_BORN;
from phonopy.structure.atoms import Atoms

from dtcs.sim.ir.SpectroscoPy.Constants import ZeroTolerance;


# ----------
# YAML Files
# ----------

def ReadPhonopyYAML(filePath = r"phonopy.yaml"):
    """
    Read a Phonopy phonopy.yaml file and extract the structure and atomic masses.

    Return value:
        A (structure, atomic_masses) tuple containing the extracted data.
        The structure is returned as a tuple of (lattice_vectors, atomic_symbols, atom_positions_frac) data.
    """

    # Load and parse YAML file.

    inputYAML = None;

    with open(filePath, 'r') as inputReader:
        inputYAML = yaml.load(inputReader, Loader = Loader);

    # Get primitive cell used in the calculation.

    primitiveCell = inputYAML['primitive_cell'];

    # Get structure.

    latticeVectors = [
        np.array(vector, dtype = np.float64)
            for vector in primitiveCell['lattice']
        ];

    atomicSymbols = [
        atom['symbol'] for atom in primitiveCell['points']
        ];

    atomPositions = [
        np.array(atom['coordinates'], dtype = np.float64)
            for atom in primitiveCell['points']
        ];

    structure = (latticeVectors, atomicSymbols, atomPositions);

    # Get atomic masses.

    atomicMasses = [
        atom['mass'] for atom in primitiveCell['points']
        ];

    # Return data.

    return (structure, atomicMasses);

def ReadMeshYAML(filePath = r"mesh.yaml"):
    """
    Read a Phonopy mesh.yaml file and extract the Gamma-point eigenvectors.

    Return value:
        A tuple of (phonon_modes, structure, atomic_masses) containing the Gamma-point frequencies/eigenvectors, and the structure and a list of atomic masses if available.
        The phonon modes are returned as a tuple of (frequencies, eigenvectors) lists with 3N elements; the eigenvectors are Nx3 NumPy arrays.
        The structure is returned as a tuple of (lattice_vectors, atomic_symbols, atom_positions_frac).

    Notes:
        If the input file does not contain frequencies and eigenvectors at q_Gamma = (0, 0, 0), errors are raised.
        YAML files exported by older versions of Phonopy do not contain the structure and atomic masses - in this case, structure and atomic_masses are set to None.
    """

    # Load and parse YAML file.

    inputYAML = None;

    with open(filePath, 'r') as inputReader:
        inputYAML = yaml.load(inputReader, Loader = Loader);

    # Get structure and atomic masses if available.

    structure, atomicMasses = None, None;

    latticeVectors, atomPositions, atomicSymbols = None, None, None;

    if 'lattice' in inputYAML:
        latticeVectors = [
            np.array(vector, dtype = np.float64)
                for vector in inputYAML['lattice']
            ];

    # Get atom positions, atomic symbols and atomic masses.

    if 'atoms' or 'points' in inputYAML:
        # Workaround for differences between older and newer versions of Phonopy...

        atomsKey = 'atoms';

        if atomsKey not in inputYAML:
            atomsKey = 'points';

        atomicSymbols = [
            atom['symbol'] for atom in inputYAML[atomsKey]
            ];

        # ... and again!

        positionKey = 'position';

        if positionKey not in inputYAML[atomsKey][0]:
            positionKey = 'coordinates';

        atomPositions = [
            np.array(atom[positionKey], dtype = np.float64)
                for atom in inputYAML[atomsKey]
            ];

        atomicMasses = [
            atom['mass'] for atom in inputYAML[atomsKey]
            ];

    # If structural information was read, create a (lattice_vectors, atomic_symbols, atom_positions) tuple.

    if latticeVectors != None and atomPositions != None and atomicSymbols != None:
        structure = (latticeVectors, atomicSymbols, atomPositions);

    # Get Gamma-point frequencies and eigenvectors.

    phononModes = None;

    frequencies, eigenvectors = None, None;

    for qPoint in inputYAML['phonon']:
        qx, qy, qz = qPoint['q-position'];

        if math.fabs(qx) < ZeroTolerance and math.fabs(qy) < ZeroTolerance and math.fabs(qz) < ZeroTolerance:
            frequencies, eigenvectors = [], [];

            for mode in qPoint['band']:
                frequencies.append(mode['frequency']);

                if 'eigenvector' not in mode:
                    raise Exception("Error: Eigenvectors not found in YAML file \"{0}\".".format(filePath));

                eigenvector = [
                     [atomDisp[i][0] for i in range(0, 3)]
                        for atomDisp in mode['eigenvector']
                    ];

                eigenvectors.append(
                    np.array(eigenvector, dtype = np.float64)
                    );

            break;

    # Check Gamma-point frequencies and eigenvectors were found.

    if frequencies == None:
        raise Exception("Error: Data for q_Gamma = (0, 0, 0) not found in YAML file \"{0}\".".format(filePath));

    phononModes = (frequencies, eigenvectors);

    # Return Gamma-point frequencies and eigenvectors, structure and atomic masses.

    return (phononModes, structure, atomicMasses);

def ReadIrRepsYAML(filePath = r"irreps.yaml"):
    """
    Read a Phonopy irreps.yaml file and extract the ir. reps. and associated band indices.

    The ir. reps. must have been calculated at the Gamma point.

    Return value:
        A list of (irrep_symbol, band_indices) tuples mapping ir. reps. to band indices.
    """

    inputYAML = None;

    with open(filePath, 'r') as inputReader:
        inputYAML = yaml.load(inputReader);

    # Check the ir. reps. are given for the Gamma point modes.

    qx, qy, qz = inputYAML['q-position'];

    if math.fabs(qx) > ZeroTolerance or math.fabs(qy) > ZeroTolerance or math.fabs(qz) > ZeroTolerance:
        raise Exception("Error: Ir. reps. are required for the Gamma point modes.");

    irRepLabels, bandIndices = [], [];

    for mode in inputYAML['normal_modes']:
        irRepLabels.append(mode['ir_label']);

        # Band indices in the irreps.yaml files start from one -> adjust to zero-based.

        bandIndices.append(
            [index - 1 for index in mode['band_indices']]
            );

    return [item for item in zip(irRepLabels, bandIndices)];


# ----------
# HDF5 Files
# ----------

def ReadMeshHDF5(filePath):
    """
    Read a Phonopy mesh.hdf5 file and return the Gamma-point eigenvectors.

    Return value:
        A tuple of (frequencies, eigenvectors) lists containing the 3N Gamma-point frequencies and eigenvectors.
        The eigenvectors are Nx3 NumPy arrays.

    Notes:
        Errors are raised if the file does not contain an entry for q_Gamma = (0, 0, 0) or does not contain eigenvectors.
    """

    frequencies, eigenvectors = None, None;

    meshHDF5 = h5py.File(filePath, 'r');

    try:
        # Locate the Gamma point in the q-point sampling mesh.

        qIndex = None;

        for i, (qx, qy, qz) in enumerate(meshHDF5['qpoint']):
            if math.fabs(qx) < ZeroTolerance and math.fabs(qy) < ZeroTolerance and math.fabs(qz) < ZeroTolerance:
                qIndex = i;
                break;

        if qIndex == None:
            raise Exception("Error: q_Gamma = (0, 0, 0) not found in the q-point mesh specified in \"{0}\".".format(filePath));

        # Extract a list of Gamma-point frequencies.

        frequencies = [
            frequency for frequency in meshHDF5['frequency'][qIndex]
            ];

        if 'eigenvector' not in meshHDF5:
            raise Exception("Error: Eigenvectors not found in \"{0}\".".format(filePath));

        # Sanity check.

        assert len(frequencies) % 3 == 0;

        nAtoms = len(frequencies) // 3;

        # Extract the Gamma-point eigenvectors.

        eigenvectors = [];

        for i in range(0, len(frequencies)):
            eigenvector = [];

            # Indexing taken from the Mesh.write_yaml() routine in Phonopy.

            for j in range(0, nAtoms):
                eigenvector.append(
                    [meshHDF5['eigenvector'][qIndex, j * 3 + k, i].real for k in range(0, 3)]
                    );

            eigenvectors.append(
                np.array(eigenvector, dtype = np.float64)
                );

    finally:
        # Make sure the input file gets closed.

        meshHDF5.close();

    # Return the Gamma-point frequencies and eigenvectors.

    return (frequencies, eigenvectors);

def ReadPhono3pyHDF5(filePath, linewidthTemperature):
    """ Read a Phono3py HDF5 file and return the Gamma-point mode linewidths at the specified temperature. """

    linewidths = None;

    phono3pyHDF5 = h5py.File(filePath, 'r');

    try:
        # Locate the temperature index corresponding to linewidthTemperature.

        tIndex = None;

        for i, t in enumerate(phono3pyHDF5['temperature']):
            if t == linewidthTemperature:
                tIndex = i;
                break;

        if tIndex == None:
            raise Exception("Error: Linewidth temperature {0} not found under the 'temperature' key.".format(linewidthTemperature));

        # If the file has a 'qpoint' key, locate the q-point index corresponding to the Gamma point.

        qIndex = None;

        if 'qpoint' in phono3pyHDF5.keys():
            for i, (qx, qy, qz) in enumerate(phono3pyHDF5['qpoint']):
                if qx == 0.0 and qy == 0.0 and qz == 0.0:
                    qIndex = i;
                    break;

            if qIndex == None:
                raise Exception("Error: q = (0, 0, ) not found under the 'qpoint' key.");

        # Extract the linewidths.

        if qIndex != None:
            # If the 'qpoint' key is present in the file, the linewidths dataset has shape (num_temperatures, num_q_points, num_bands).

            linewidths = [
                linewidth for linewidth in phono3pyHDF5['gamma'][tIndex][qIndex]
                ];
        else:
            # If not, the dataset has shape (num_temperatures, num_bands).

            linewidths = [
                linewidth for linewidth in phono3pyHDF5['gamma'][tIndex]
                ];
    finally:
        # Make sure we close the input file.

        phono3pyHDF5.close();

    return linewidths;


# ----------
# Misc Files
# ----------

def ReadBORN(structure, filePath = "BORN"):
    """
    Read a Phonopy BORN file, and expand the charges for the supplied structure, and return a list of Born effective-charge tensors for each atom in the structure.

    Arguments:
        structure -- a tuple of (lattice_vectors, atomic_symbols, atom_positions) specifying the structure to be used to expand the Born charges to the full set of atoms.

    Keyword arguments:
        filePath -- path to the Phonopy BORN-format file to read (default: BORN).

    Return value:
        A list of 3x3 Born effective-charge tensors for each atom in the supplied structure.

    Notes:
        The atom positions supplied with structure are assumed to be in fractional coordinates.
    """

    # Convert the supplied structure to a Phonopy Atoms object.

    latticeVectors, atomicSymbols, atomPositions = structure;

    cell = Atoms(
        symbols = atomicSymbols, cell = latticeVectors, scaled_positions = atomPositions
        );

    # Read the BORN file.

    bornData = parse_BORN(cell, filename = filePath);

    # Return the Born effective-charge tensors from the dictionary returned by parse_BORN.

    return [becTensor for becTensor in bornData['born']];
