class TimeProfiler():
    '''
    Class for timing the in-simulation run times of surface CRN simulations.
    '''

    def __init__(self, simulation):
        self.simulation = simulation

    def run_simulations(self, num_runs, stop_criteria, output_file):
        '''
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
        '''
        with open(output_file, 'w') as outstream:
            outstream.write("Run Times")
            for sim in self.simulation_results(num_runs, stop_criteria):
                outstream.write("\n" + str(self.simulation.time))

    def simulation_results(self, num_runs, stop_criteria):
        '''
        Generator yielding the end results of repeated simulations.

        Parameters:
            -num_runs: the number of times to run the simulation.
            -stop_criteria: a function taking a simulation and returning True if 
                                the simulation is over according to some 
                                criteria and False otherwisey. 
        '''
        for i in range(num_runs):
            self.simulation.reset()
            while not (self.simulation.done() or
                       stop_criteria(self.simulation)):
                self.simulation.process_next_reaction()
            yield self.simulation