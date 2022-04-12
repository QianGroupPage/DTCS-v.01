from dtcs.sim.ir import EigenvectorsToEigendisplacements

from dtcs.sim.ir import ReadBORN

from dtcs.sim.ir.SpectroscoPy.CLI.Parser import PostProcessArgs
from dtcs.sim.ir.SpectroscoPy.CLI.Phonopy import Phonopy_LoadData_Core, Phonopy_LoadData_Optional
from dtcs.sim.ir import RunMode_IR


def ir_sim_wrapper(
        MeshHDF5='mesh.hdf5',
        MeshYAML='mesh.yaml',
        PhonopyYAML='phonopy.yaml',
        CellFile='POSCAR',
        BORNFile='BORN',
        LinewidthHDF5=None,
        LinewidthTemperature=None,
        IrRepsYAML='irreps.yaml',
        FrequencyUnits='inv_cm',
        IrReps=False,
        Linewidth=16.5,
        SpectrumRange=None,
        SpectrumResolution=None,
        InstrumentBroadening=None,
        InstrumentBroadeningShape='gaussian',
        OutputPrefix=None,
        DataFormat='dat',
):

    args = lambda: None

    # Emulate Phonopy_UpdateParser(parser, 'ir');
    args.MeshHDF5 = MeshHDF5  # 'mesh.hdf5'
    args.MeshYAML = MeshYAML  # 'mesh.yaml'
    args.PhonopyYAML = PhonopyYAML  # 'phonopy.yaml'
    args.CellFile = CellFile  # 'POSCAR'
    args.BORNFile = BORNFile  # 'BORN'
    args.LinewidthHDF5 = LinewidthHDF5  # None
    args.LinewidthTemperature = LinewidthTemperature  # None
    args.IrRepsYAML = IrRepsYAML  # 'irreps.yaml'

    # Emulate UpdateParser(parser, 'ir', supportedFeatures = ['ir_reps'])
    args.FrequencyUnits = FrequencyUnits  # 'inv_cm'
    #  Peak table
    args.IrReps = IrReps  # False
    #  Spectrum simulation
    args.Linewidth = Linewidth  # 16.5
    args.SpectrumRange = SpectrumRange  # None
    args.SpectrumResolution = SpectrumResolution  # None
    args.InstrumentBroadening = InstrumentBroadening  # None
    args.InstrumentBroadeningShape = InstrumentBroadeningShape  # 'gaussian'
    #  Data output
    args.OutputPrefix = OutputPrefix  # None
    args.DataFormat = DataFormat  # 'dat'


    PostProcessArgs(args, 'ir');

    # Read input data.

    inputData = Phonopy_LoadData_Core(
        args, extractList = ['structure', 'atomic_masses', 'phonon_modes']
    )

    structure = inputData['structure']
    atomicMasses = inputData['atomic_masses']

    frequencies, eigenvectors = inputData['phonon_modes']

    # Convert eigenvectors to eigendisplacements.

    eigendisplacements = EigenvectorsToEigendisplacements(
        eigenvectors, atomicMasses
    )

    # Read Born effective-charge tensors.

    becTensors = ReadBORN(structure, filePath = args.BORNFile)

    # Read ir. rep. data and/or linewidths, if required.

    linewidths, irRepData = Phonopy_LoadData_Optional(args)

    RunMode_IR(
        frequencies=frequencies,

        frequencyUnits='thz',

        eigendisplacements=eigendisplacements,

        becTensors=becTensors,
        args=args,

        linewidths=linewidths,
        irRepData=irRepData,
    )