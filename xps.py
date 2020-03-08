first_line = '[region 2]'
start_info = '[info'
start_data = '[data'
info_keys = {'region name', 'center energy'}

def read_data(filename):
    f = open(filename)
    lines = []
    for l in f:
        lines.append(l)
    f.close()

    it = iter(lines)
    line = next(it, None)

    data = []
    # Skip the survey scan
    while line is not None and line.strip().lower() != first_line:
        line = next(it, None)

    while line is not None:
        data.append(parse_region(it))
        line = next(it, None)

    return data

def parse_region(it):
    info = {}
    x, y = [], []
    line = next(it, None)

    # Skip to info section
    while line is not None:
        line = line.strip().lower()
        if start_info in line:
            break

        line = next(it, None)

    while line is not None:
        line = line.strip().lower()
        # Parse a region
        if start_data in line:
            x, y = parse_data(it)
            break

        if line[:line.find('=')] in info_keys:
            info[line[:line.find('=')]] = line[line.find('=')+1:]

        line = next(it, None)

    return XPSData(info, x, y)

def parse_data(it):
    x, y = [], []
    line = next(it, None)

    while line is not None:
        line = line.strip()
        if not line:
           break 

        nums = line.split()

        x.append(float(nums[0]))
        y.append(float(nums[1]))
        line = next(it, None)

    return x, y

class XPSData:
    def __init__(self, info, binding_energy, intensity):
        self.info = info
        self.binding_energy = binding_energy
        self.intensity = intensity
