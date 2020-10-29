class Condition:
    def __init__(self, line_block):
        """
        Read information from a block in the notebook provided by the experimentalist.
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
        Set the measurements list of this condition.
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
        return " ".join([self.pressure, self.temperature])

    @staticmethod
    def read_condition_line(line):
        """
        Condition line consists of P, T, comments, each separated by any of ",", ":", or space.
        """
        line = line.rstrip()

        allowed_separators = [",", ". ", ":"]
        for separator in allowed_separators:
            line = line.replace(separator, ",")
        return [line_section.strip() for line_section in line.split(",")]

    @staticmethod
    def read_measurement_entry_line(line):
        """
        A line for a measurement entry consists of sequence number for the entry
        :param line:
        :return:
        """
        sequence_number = line.split()[0]
        line_left = "".join(line.split()[1:])
        comments = line_left.split(":")[1:]
        species = line_left.split(":")[0].split()
        return int(sequence_number), species, comments
