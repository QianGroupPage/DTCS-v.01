from __future__ import print_function

import itertools

import lblcrn.sim.surface_crn.surface_crns.models.grids as grids
import lblcrn.sim.surface_crn.surface_crns.readers.manifest_readers as manifest_readers
import lblcrn.sim.surface_crn.surface_crns.simulators.queue_simulator as simulators
from lblcrn.sim.surface_crn.surface_crns.profiling import TimeProfiler


class GridSimTimeProfiler():
    """
    Class for profiling grid-like surface CRN simulations. Can profile run times
    and output at specified locations for one or more initial states.
    """

    def __init__(self, manifest_filename=None, options=None):
        self.DEBUG = False

        if manifest_filename and options:
            raise Exception("GridSimTimeProfiler can be initialized with " +
                            "either a manifest file or a dictionary of " +
                            "options, but not both.")
        if manifest_filename:
            self.options = manifest_readers.read_manifest(manifest_filename)
        elif options:
            self.options = options
        else:
            raise Exception("GridSimTimeProfiler must be initialized with " +
                            "either a manifest file name or a dictionary " +
                            "of options read from a manifest file.")

        # Read out miscellaneous options from manifest
        if 'debug' in self.options:
            debug_flag = self.options['debug'].lower()
            if debug_flag in ['true', 'on', 'yes']:
                self.DEBUG = True
            elif debug_flag in ['false', 'off', 'no']:
                self.DEBUG = False
            else:
                self.DEBUG = bool(int(self.options['debug']))
        else:
            DEBUG = debug
        if 'rng_seed' in self.options:
            self.RAND_SEED = int(self.options['rng_seed'])
        else:
            self.RAND_SEED = 0
        if 'max_duration' in self.options:
            self.SIMULATION_DURATION = float(self.options['max_duration'])
        else:
            self.SIMULATION_DURATION = 100

        # Load transition rules
        if 'transition_rules' in self.options:
            self.transition_rules = self.options['transition_rules']
            # Check validity of rules and add colors if necessary.
            for rule in self.transition_rules:
                # Check that reactions are unimolecular or bimolecular
                if len(rule.inputs) != len(rule.outputs):
                    raise Exception("Invalid transition rule " + str(rule) +
                                    "\nTransition rules for a surface CRN must" +
                                    "have the same number of inputs and " +
                                    "outputs.")
                if len(rule.inputs) > 2:
                    raise Exception("Invalid transition rule " + str(rule) +
                                    "\nOnly transition rules between one or " +
                                    "two species are allowed in this " +
                                    "implementation.")
            if self.DEBUG:
                print("Transition rules:\n" + str(self.transition_rules))
        else:
            raise Exception("Transition rules required in manifest.")

        # Set up grid and define initial state
        if 'init_state' in self.options:
            self.initial_state = self.options['init_state']
            if self.DEBUG:
                print("Initial state:")
                print(str(self.initial_state))
        else:
            raise Exception("Initial state required in manifest")

        self.grid = grids.SquareGrid(self.initial_state.shape[0],
                                     self.initial_state.shape[1])
        self.grid.set_global_state(self.initial_state)

        # Set up simulation timer
        self.simulation = simulators.QueueSimulator(surface=self.grid,
                                                    transition_rules=self.transition_rules,
                                                    seed=self.RAND_SEED,
                                                    simulation_duration=self.SIMULATION_DURATION)

    def timing_test(self, num_runs, stop_criteria, output_file):
        """
        Run the simulation num_runs times, stopping according to user-defined
        stop criteria, and printing simulation times to completion

        Parameters:
            -num_runs: the number of times to run the simulation
            -stop_criteria: a function taking a simulation and returning True if
                                the simulation is over according to some
                                criteria and False if the simulation is not.
            -output_file: the name of a file to print the output to. Output is a
                            TSV file (readable by gnuplot) with times for each
                            run.
        """
        timer = TimeProfiler(self.simulation)
        timer.run_simulations(num_runs, stop_criteria, output_file)

    def correctness_test(self, num_runs, start_states, start_state_labels,
                         stop_criteria, output_function, output_labels,
                         output_file):
        """
        Run the simulation multiple times on different inputs, printing the
        times to completion and the states of any number of output positions at
        the end of each simulation.

        Parameters:
            -num_runs: the number of times to run the simulation from each
                        initial condition.
            -start_states: a list or generator of initial states, in the form
                            of a numpy array (or list of lists). The simulation
                            will be run for each starting state num_runs times.
            -start_state_labels: list or generator of names for the states in
                                    start_states.
            -stop_criteria: a function taking a simulation and returning True if
                                the simulation is over according to some
                                criteria and False otherwise.
            -output_function: A function taking a simulation and returning a
                                list of outputs to be reported.
            -output_labels: list or generator of names for teh outputs of
                                output_function
            -output_file: the name of a file to print the output to. Output is a
                            TSV file (readable by gnuplot) with times for each
                            run.
        """
        with open(output_file, 'w') as outstream:
            outstream.write("Start State\tRun Number\tTime")
            states_and_labels = itertools.izip(start_states, start_state_labels)
            for label in output_labels:
                outstream.write("\t" + label)
            x = 0
            for init_state, state_label in states_and_labels:
                print("x: " + str(x))
                x += 1
                grid = grids.SquareGrid(init_state.shape[0],
                                        init_state.shape[1])
                grid.set_global_state(init_state)
                self.simulation = simulators.QueueSimulator(
                    surface=grid,
                    transition_rules=self.transition_rules,
                    seed=self.RAND_SEED,
                    simulation_duration=self.SIMULATION_DURATION)
                timer = TimeProfiler(self.simulation)
                run_num = 0
                for sim_end in timer.simulation_results(num_runs, stop_criteria):
                    run_num += 1
                    outstream.write("\n")
                    outstream.write(state_label)
                    outstream.write("\t" + str(run_num))
                    outstream.write("\t" + str(sim_end.time))
                    for output in output_function(sim_end):
                        outstream.write("\t" + output)
