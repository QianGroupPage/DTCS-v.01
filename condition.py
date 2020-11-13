class Condition:
    def __init__(self, line_block):
        """
        :param line_block: a list of lines describing a set of experiments for under a condition.
                           The first line in the line-block describes the condition for the block;
                           Each of the following line describes an individual set of experimental measurements.
        """
        condition_info = self.read_condition_line(line_block[0])
        self.pressure, self.temperature = condition_info[:2]
        self.comments = " ".join(condition_info[2:])

        self.measurement_sequence_numbers = []
        self.measurement_info_dict = {}
        for line in line_block[1:]:
            sequence_number, species, comments = self.read_measurement_entry_line(line)
            self.measurement_sequence_numbers.append(sequence_number)
            self.measurement_info_dict[sequence_number] = [species, comments]

        self.measurements = []

    def set_measurements(self, measurements):
        """
        This is called internally when the list of measurements for this condition has been processed.
        :param measurements: a list of RawMeasurement objects.
        :return: None
        """
        self.measurements = measurements

    def list_region_names(self):
        """
        :return: a list enumerating all the region names included in all measurements in this condition.
        """
        region_names = []
        for measurement in self.measurements:
            region_names.extend(measurement.list_region_names())
        return region_names

    @property
    def id(self):
        """
        :return: a unique identifier for this condition, temperature followed by pressure.
        """
        return " ".join([self.pressure, self.temperature])

    @staticmethod
    def read_condition_line(line):
        """
        :param line: a line consisting of P, T, comments, each separated by any of ",", ":", or space.
        :return: a list of strings describing measurement conditions.
        """
        line = line.rstrip()

        allowed_separators = [",", ". ", ":"]
        for separator in allowed_separators:
            line = line.replace(separator, ",")
        return [line_section.strip() for line_section in line.split(",")]

    @staticmethod
    def read_measurement_entry_line(line):
        """
        :param line: a line that describes the measurement entry, starting with its sequence number.
        :return: the sequence number of the measurement, a list of species names in it, a list of comment strings
                 for the measurement.
        """
        sequence_number = line.split()[0]
        line_left = "".join(line.split()[1:])

        comments = []
        for c in line_left.split(":")[1:]:
            comments.extend(c.split(","))
        species = line_left.split(":")[0].split(',')
        return int(sequence_number), species, comments
